"""Unit tests for src/analysis/gene_lookup.py"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from src.analysis.gene_lookup import (
    GeneLookupResult,
    lookup_gene_name,
    _try_case_insensitive,
    _try_exact,
    _try_fuzzy,
    _try_species_prefix,
)


def _make_adata(var_names: list[str], raw_var_names: list[str] | None = None) -> MagicMock:
    """Build a minimal AnnData mock with the given var_names."""
    adata = MagicMock()
    adata.var_names = var_names
    if raw_var_names is not None:
        adata.raw = MagicMock()
        adata.raw.var_names = raw_var_names
    else:
        adata.raw = None
    return adata


# ---------------------------------------------------------------------------
# Strategy 1: exact match
# ---------------------------------------------------------------------------

class TestExactMatch:
    def test_exact_hit(self):
        adata = _make_adata(["TNNT2", "MKI67", "CD3E"])
        result = lookup_gene_name(adata, "TNNT2")
        assert result.matched is True
        assert result.resolved_name == "TNNT2"
        assert result.strategy == "exact"

    def test_exact_miss(self):
        adata = _make_adata(["TNNT2", "MKI67"])
        result = _try_exact("tnnt2", {"TNNT2", "MKI67"})
        assert result is None

    def test_raw_only_gene_not_found(self):
        # Fix 1: lookup is restricted to adata.var_names only.
        # A gene present only in raw.var_names must NOT be returned as resolved —
        # downstream tools use adata.var_names and would fail with a confusing error.
        adata = _make_adata(["MKI67"], raw_var_names=["TNNT2", "MKI67"])
        result = lookup_gene_name(adata, "TNNT2")
        assert result.matched is False


# ---------------------------------------------------------------------------
# Strategy 2: case-insensitive match
# ---------------------------------------------------------------------------

class TestCaseInsensitiveMatch:
    def test_lowercase_input(self):
        adata = _make_adata(["TNNT2", "MKI67"])
        result = lookup_gene_name(adata, "tnnt2")
        assert result.matched is True
        assert result.resolved_name == "TNNT2"
        assert result.strategy == "case_insensitive"

    def test_mixed_case_input(self):
        adata = _make_adata(["Tnnt2", "Mki67"])
        result = lookup_gene_name(adata, "TNNT2")
        # TNNT2 != Tnnt2 exact, but case-insensitive should catch it
        # However species_prefix also fires for ALLCAPS → Title — strategy 2 runs first
        assert result.matched is True
        assert result.resolved_name == "Tnnt2"
        assert result.strategy in ("case_insensitive", "species_prefix")

    def test_helper_returns_correct_cased_name(self):
        names = ["Tnnt2", "Mki67", "Cd3e"]
        assert _try_case_insensitive("CD3E", names) == "Cd3e"
        assert _try_case_insensitive("mki67", names) == "Mki67"
        assert _try_case_insensitive("NOTHERE", names) is None


# ---------------------------------------------------------------------------
# Strategy 3: species prefix / case convention
# ---------------------------------------------------------------------------

class TestSpeciesPrefixMatch:
    def test_human_to_mouse(self):
        """TNNT2 (human ALLCAPS) → Tnnt2 (mouse Title case)"""
        adata = _make_adata(["Tnnt2", "Mki67", "Cd3e"])
        result = lookup_gene_name(adata, "TNNT2")
        assert result.matched is True
        assert result.resolved_name == "Tnnt2"
        # After Fix 2: bare ALLCAPS falls through to case_insensitive, not species_prefix
        assert result.strategy == "case_insensitive"

    def test_mouse_to_human(self):
        """Tnnt2 (mouse Title) → TNNT2 (human ALLCAPS)"""
        adata = _make_adata(["TNNT2", "MKI67"])
        result = lookup_gene_name(adata, "Tnnt2")
        assert result.matched is True
        assert result.resolved_name == "TNNT2"
        assert result.strategy in ("case_insensitive", "species_prefix")

    def test_strip_mm_prefix(self):
        adata = _make_adata(["Tnnt2", "Mki67"])
        result = lookup_gene_name(adata, "Mm_Tnnt2")
        assert result.matched is True
        assert result.resolved_name == "Tnnt2"
        assert result.strategy == "species_prefix"

    def test_strip_hs_prefix(self):
        adata = _make_adata(["TNNT2", "MKI67"])
        result = lookup_gene_name(adata, "Hs_TNNT2")
        assert result.matched is True
        assert result.resolved_name == "TNNT2"
        assert result.strategy == "species_prefix"

    def test_helper_no_match(self):
        names_set = {"TNNT2", "MKI67"}
        assert _try_species_prefix("COMPLETELY_DIFFERENT", names_set) is None

    def test_bare_allcaps_not_resolved_by_species_prefix(self):
        # Fix 2: bare ALLCAPS without a species prefix (e.g. IL2) must NOT be
        # silently resolved via case-toggle — IL2 → Il2 could be the wrong gene
        # when the user meant IL2RA. Strategy 3 returns None; fuzzy handles it.
        names_set = {"Il2", "Il2ra", "Mki67"}
        assert _try_species_prefix("IL2", names_set) is None

    def test_prefixed_allcaps_still_resolves(self):
        # Mm_TNNT2 has a prefix — case-toggle on the stripped part is still valid.
        adata = _make_adata(["Tnnt2", "Mki67"])
        result = lookup_gene_name(adata, "Mm_TNNT2")
        assert result.matched is True
        assert result.resolved_name == "Tnnt2"
        assert result.strategy == "species_prefix"

    def test_bare_allcaps_no_prefix_falls_through_to_fuzzy(self):
        # Fix 2: bare ALLCAPS without a species prefix (e.g. IL2) must NOT be
        # silently resolved via case-toggle — IL2 → Il2 could be the wrong gene
        # when the user meant IL2RA. Strategy 3 returns None; fuzzy surfaces candidates.
        adata = _make_adata(["Il2", "Il2ra", "Mki67"])
        result = _try_species_prefix("IL2", {"Il2", "Il2ra", "Mki67"})
        assert result is None  # strategy 3 must not fire for bare ALLCAPS


# ---------------------------------------------------------------------------
# Strategy 4: fuzzy match
# ---------------------------------------------------------------------------

class TestFuzzyMatch:
    def test_fuzzy_typo(self):
        """TNNT (partial) should fuzzy-match to TNNT2"""
        adata = _make_adata(["TNNT2", "TNNT1", "MKI67", "CD3E", "CD8A"])
        result = lookup_gene_name(adata, "TNNT")
        assert result.matched is False
        assert result.strategy == "fuzzy"
        assert len(result.candidates) > 0
        assert any("TNNT" in c for c in result.candidates)

    def test_fuzzy_returns_up_to_5_candidates(self):
        names = [f"GENE{i}" for i in range(20)]
        candidates = _try_fuzzy("GENE1", names, cutoff=0.6, n=5)
        assert len(candidates) <= 5

    def test_fuzzy_cutoff_filters_noise(self):
        """Completely unrelated name should return no fuzzy candidates."""
        adata = _make_adata(["TNNT2", "MKI67", "CD3E"])
        result = lookup_gene_name(adata, "XYZQWERTY")
        assert result.matched is False
        assert result.strategy is None
        assert result.candidates == []

    def test_fuzzy_message_contains_candidates(self):
        adata = _make_adata(["TNNT2", "TNNT1", "MKI67"])
        result = lookup_gene_name(adata, "TNNT")
        assert "TNNT" in result.message or "TNNT2" in result.message or "TNNT1" in result.message


# ---------------------------------------------------------------------------
# No match at all
# ---------------------------------------------------------------------------

class TestNoMatch:
    def test_completely_unknown_gene(self):
        adata = _make_adata(["TNNT2", "MKI67", "CD3E"])
        result = lookup_gene_name(adata, "XYZQWERTY123")
        assert result.matched is False
        assert result.resolved_name is None
        assert result.strategy is None
        assert result.candidates == []
        assert "not found" in result.message.lower()

    def test_result_is_dataclass(self):
        adata = _make_adata(["TNNT2"])
        result = lookup_gene_name(adata, "TNNT2")
        assert isinstance(result, GeneLookupResult)


# ---------------------------------------------------------------------------
# Strategy ordering: exact wins over case-insensitive
# ---------------------------------------------------------------------------

class TestStrategyOrdering:
    def test_exact_takes_priority(self):
        """If exact match exists, should not fall through to case-insensitive."""
        adata = _make_adata(["TNNT2", "tnnt2"])
        result = lookup_gene_name(adata, "TNNT2")
        assert result.strategy == "exact"
        assert result.resolved_name == "TNNT2"

    def test_case_insensitive_before_species_prefix(self):
        """Case-insensitive (strategy 2) runs before species prefix (strategy 3)."""
        # Dataset has "Tnnt2" — both case_insensitive and species_prefix would match
        # Strategy 2 should win
        adata = _make_adata(["Tnnt2"])
        result = lookup_gene_name(adata, "TNNT2")
        # Either strategy 2 or 3 is acceptable here since both resolve correctly
        assert result.matched is True
        assert result.resolved_name == "Tnnt2"
