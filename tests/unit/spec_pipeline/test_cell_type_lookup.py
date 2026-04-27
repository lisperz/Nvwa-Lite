"""Unit tests for lookup_cell_type_name() — exhaustive branch coverage.

Covers every matching strategy in src/domain/resolver/cell_type_lookup.py:
  - Exact match (literal string, case-sensitive)
  - Normalized match (lowercase + whitespace/dash/underscore collapse)
  - Initials match (uppercase letters only — NK preservation rule)
  - Initials match (first-letter of each word — EC first-letter rule)
  - Unambiguous / ambiguous / no-match outcomes
  - Missing cell_type column
  - Edge inputs: whitespace-only, empty string
"""

from __future__ import annotations

import pytest

from src.domain.resolver.cell_type_lookup import CellTypeLookupResult, lookup_cell_type_name


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_STANDARD_LABELS = [
    "T cell",
    "B cell",
    "NK cell",
    "endothelial cell",
    "epithelial cell",
    "Dendritic Cell",
]
"""Labels used in most tests; inlined into adata.obs['cell_type']."""


def _set_cell_type(adata, labels: list[str]) -> None:
    """Assign a repeating sequence of *labels* to adata.obs['cell_type']."""
    import numpy as np
    import pandas as pd

    n = adata.n_obs
    values = [labels[i % len(labels)] for i in range(n)]
    adata.obs["cell_type"] = pd.Categorical(values)


# ---------------------------------------------------------------------------
# 1. Exact match (strategy == "exact")
# ---------------------------------------------------------------------------


class TestExactMatch:
    def test_exact_match_returns_matched_true(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T cell")
        assert result.matched is True
        assert result.resolved_name == "T cell"
        assert result.strategy == "exact"
        assert result.input_name == "T cell"

    def test_exact_match_preserves_input_casing(self, adata):
        # "Dendritic Cell" is stored with that casing — must match exactly
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "Dendritic Cell")
        assert result.matched is True
        assert result.strategy == "exact"
        assert result.resolved_name == "Dendritic Cell"

    def test_exact_match_wrong_case_does_not_use_exact(self, adata):
        # "t cell" is NOT in the list, so exact branch should not fire
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "t cell")
        # Falls through to normalized — still expected to match but NOT strategy=exact
        assert result.strategy != "exact"

    def test_exact_match_nk_cell(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "NK cell")
        assert result.matched is True
        assert result.strategy == "exact"
        assert result.resolved_name == "NK cell"


# ---------------------------------------------------------------------------
# 2. Normalized match (strategy == "normalized")
# ---------------------------------------------------------------------------


class TestNormalizedMatch:
    def test_lowercase_resolves_to_canonical(self, adata):
        # "T cell" stored; "t cell" should normalize-match
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "t cell")
        assert result.matched is True
        assert result.strategy == "normalized"
        assert result.resolved_name == "T cell"

    def test_trailing_whitespace_stripped(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "B cell  ")
        assert result.matched is True
        assert result.strategy == "normalized"
        assert result.resolved_name == "B cell"

    def test_leading_whitespace_stripped(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "  B cell")
        assert result.matched is True
        assert result.strategy == "normalized"
        assert result.resolved_name == "B cell"

    def test_underscore_treated_as_space(self, adata):
        # _normalize converts "_" → "_", so "T_cell" normalizes to "t_cell"
        # and "T cell" also normalizes to "t_cell" — should match
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T_cell")
        assert result.matched is True
        assert result.strategy == "normalized"
        assert result.resolved_name == "T cell"

    def test_hyphen_treated_as_separator(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T-cell")
        assert result.matched is True
        assert result.strategy == "normalized"
        assert result.resolved_name == "T cell"

    def test_mixed_case_with_separator_variants(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "b-Cell")
        assert result.matched is True
        assert result.strategy == "normalized"
        assert result.resolved_name == "B cell"

    def test_normalized_result_contains_original_input(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "t_cell")
        assert result.input_name == "t_cell"
        assert result.resolved_name == "T cell"

    def test_normalized_candidates_empty_on_match(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "t cell")
        assert result.candidates == []


# ---------------------------------------------------------------------------
# 3. Initials matching — uppercase letters only (NK preservation rule)
#    _initials_uppercase extracts all uppercase letters from the name.
#    "NK cell" → uppercase letters: N, K → "NK"
# ---------------------------------------------------------------------------


class TestInitialsUppercase:
    def test_nk_resolves_unambiguously(self, adata):
        # "NK cell" is the only name whose uppercase letters spell "NK"
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "NK")
        assert result.matched is True
        assert result.strategy == "initials_uppercase"
        assert result.resolved_name == "NK cell"

    def test_nk_input_name_preserved(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "NK")
        assert result.input_name == "NK"

    def test_nk_candidates_empty_on_unambiguous_match(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "NK")
        assert result.candidates == []

    def test_initials_uppercase_ambiguous_returns_candidates(self, adata):
        # Add two names that both have uppercase letters spelling "DC"
        labels = ["Dendritic Cell", "Dust Cell", "T cell"]
        # "Dendritic Cell" → DC; "Dust Cell" → DC
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "DC")
        assert result.matched is False
        assert result.strategy == "initials_uppercase"
        assert len(result.candidates) == 2
        assert set(result.candidates) == {"Dendritic Cell", "Dust Cell"}

    def test_initials_uppercase_no_match_falls_through(self, adata):
        # "XY" has no uppercase-initials match in standard labels; falls to
        # first-letter check and ultimately no-match
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "XY")
        assert result.matched is False
        assert result.strategy is None

    def test_initials_check_requires_all_uppercase(self, adata):
        # "nk" is NOT all-uppercase, so initials branch should not fire
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "nk")
        # Should NOT be strategy=initials_uppercase; may match via normalized or not at all
        assert result.strategy != "initials_uppercase"

    def test_initials_check_requires_len_le_3(self, adata):
        # 4-character all-uppercase string skips the initials branch entirely
        labels = ["Natural Killer Lymphocyte Cell", "T cell"]
        # uppercase initials of first label = "NKLC" — but len("NKLC") > 3
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "NKLC")
        # initials branch requires len(raw) <= 3, so no initials strategy
        assert result.strategy != "initials_uppercase"
        assert result.strategy != "initials_firstletter"
        assert result.matched is False

    def test_single_uppercase_letter_can_match(self, adata):
        # len("B") == 1, isupper() == True — initials branch fires
        # "B cell" has uppercase "B" only → uppercase initials = "B"
        labels = ["B cell", "T cell"]
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "B")
        assert result.matched is True
        assert result.strategy == "initials_uppercase"
        assert result.resolved_name == "B cell"


# ---------------------------------------------------------------------------
# 4. Initials matching — first-letter of each word (EC / first-letter rule)
#    _initials_first_letters splits on whitespace/underscore/hyphen and
#    takes the first letter (uppercased) of each token.
#    "endothelial cell" → E, C → "EC"
#    "epithelial cell"  → E, C → "EC"  — ambiguous
# ---------------------------------------------------------------------------


class TestInitialsFirstLetter:
    def test_ec_matches_endothelial_and_epithelial_ambiguous(self, adata):
        labels = ["endothelial cell", "epithelial cell", "T cell"]
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "EC")
        assert result.matched is False
        assert result.strategy == "initials_firstletter"
        assert len(result.candidates) == 2
        assert set(result.candidates) == {"endothelial cell", "epithelial cell"}

    def test_ec_ambiguous_candidates_not_empty(self, adata):
        labels = ["endothelial cell", "epithelial cell"]
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "EC")
        assert result.candidates  # non-empty

    def test_first_letter_unambiguous_returns_matched(self, adata):
        # "TC" → T cell only (endothelial has uppercase "E" but _initials_uppercase
        # for "TC" fires first — we need "TC" to match via first-letter only).
        # "T cell" → first letters T, C → "TC"
        # No other label starts with T-something-C
        labels = ["T cell", "B cell", "NK cell"]
        _set_cell_type(adata, labels)
        # "TC" not all-uppercase for initials_uppercase match (it is all uppercase)
        # _initials_uppercase("T cell") = "T" (only 'T' is uppercase)
        # so uppercase branch fires first, finds "T" not "TC" — no uppercase match
        # then first-letter branch: "T cell" → "TC" — unambiguous
        result = lookup_cell_type_name(adata, "TC")
        assert result.matched is True
        assert result.strategy == "initials_firstletter"
        assert result.resolved_name == "T cell"

    def test_first_letter_initials_input_name_preserved(self, adata):
        labels = ["endothelial cell", "T cell"]
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "EC")
        assert result.input_name == "EC"

    def test_first_letter_falls_through_after_uppercase_no_match(self, adata):
        # "BC" — uppercase initials of "B cell" = "B" (not "BC");
        # first-letter of "B cell" = "BC" — should match
        labels = ["B cell", "endothelial cell"]
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "BC")
        assert result.matched is True
        assert result.strategy == "initials_firstletter"
        assert result.resolved_name == "B cell"

    def test_hyphen_separated_name_in_first_letter_initials(self, adata):
        # Stored label uses hyphen separator
        labels = ["natural-killer cell", "T cell"]
        _set_cell_type(adata, labels)
        # first letters of "natural-killer cell" → N, K, C → "NKC"
        result = lookup_cell_type_name(adata, "NKC")
        assert result.matched is True
        assert result.strategy == "initials_firstletter"
        assert result.resolved_name == "natural-killer cell"


# ---------------------------------------------------------------------------
# 5. No match outcome
# ---------------------------------------------------------------------------


class TestNoMatch:
    def test_no_match_returns_matched_false(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "mast cell")
        assert result.matched is False
        assert result.resolved_name is None
        assert result.strategy is None

    def test_no_match_candidates_empty(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "mast cell")
        assert result.candidates == []

    def test_no_match_input_name_preserved(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "mast cell")
        assert result.input_name == "mast cell"

    def test_no_match_message_nonempty(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "mast cell")
        assert result.message  # not empty

    def test_partially_similar_name_does_not_match(self, adata):
        # "T" alone: len=1, isupper=True — fires initials branch
        # "T cell" → _initials_uppercase = "T" → matches "T cell"
        # This is actually a legitimate upstream match, so test a truly absent name
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "plasma cell")
        assert result.matched is False
        assert result.strategy is None


# ---------------------------------------------------------------------------
# 6. Missing cell_type column
# ---------------------------------------------------------------------------


class TestMissingColumn:
    def test_no_cell_type_column_returns_matched_false(self, adata):
        # Do not set any cell_type column — use raw adata from fixture
        # Ensure the column is absent
        if "cell_type" in adata.obs.columns:
            del adata.obs["cell_type"]
        result = lookup_cell_type_name(adata, "T cell")
        assert result.matched is False
        assert result.resolved_name is None
        assert result.strategy is None

    def test_no_cell_type_column_message_mentions_column(self, adata):
        if "cell_type" in adata.obs.columns:
            del adata.obs["cell_type"]
        result = lookup_cell_type_name(adata, "T cell")
        assert "cell_type" in result.message

    def test_no_cell_type_column_input_name_preserved(self, adata):
        if "cell_type" in adata.obs.columns:
            del adata.obs["cell_type"]
        result = lookup_cell_type_name(adata, "T cell")
        assert result.input_name == "T cell"

    def test_custom_obs_col_missing_returns_matched_false(self, adata):
        # obs_col="leiden" is not present → graceful failure
        result = lookup_cell_type_name(adata, "0", obs_col="leiden")
        assert result.matched is False
        assert result.strategy is None

    def test_custom_obs_col_present_matches(self, adata):
        import pandas as pd

        adata.obs["leiden"] = pd.Categorical(["0"] * adata.n_obs)
        result = lookup_cell_type_name(adata, "0", obs_col="leiden")
        assert result.matched is True
        assert result.strategy == "exact"


# ---------------------------------------------------------------------------
# 7. Edge inputs — whitespace-only and empty string
# ---------------------------------------------------------------------------


class TestEdgeInputs:
    def test_whitespace_only_raw_no_match(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "   ")
        # "   " is not in all_names (exact fails)
        # _normalize("   ") = "" after strip → "" not in norm_map unless there's an empty-string label
        # len("   ") = 3 but NOT isupper() → initials branch skipped
        assert result.matched is False

    def test_empty_string_raw_no_match(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "")
        # len("") = 0, not isupper() → initials branch skipped
        # "" not in all_names; _normalize("") = "" not in norm_map
        assert result.matched is False

    def test_empty_string_input_name_preserved(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "")
        assert result.input_name == ""

    def test_whitespace_only_input_name_preserved(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "   ")
        assert result.input_name == "   "

    def test_single_space_no_match(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, " ")
        assert result.matched is False


# ---------------------------------------------------------------------------
# 8. Result type and field contract
# ---------------------------------------------------------------------------


class TestResultType:
    def test_returns_cell_type_lookup_result_instance(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T cell")
        assert isinstance(result, CellTypeLookupResult)

    def test_candidates_default_is_list(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T cell")
        assert isinstance(result.candidates, list)

    def test_message_default_is_str(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T cell")
        assert isinstance(result.message, str)

    def test_matched_is_bool(self, adata):
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T cell")
        assert isinstance(result.matched, bool)


# ---------------------------------------------------------------------------
# 9. Strategy priority — exact fires before normalized
# ---------------------------------------------------------------------------


class TestStrategyPriority:
    def test_exact_beats_normalized_when_both_would_match(self, adata):
        # "T cell" stored; raw="T cell" would match both exact and normalized
        # exact must fire first (implementation: exact check is lines 37-41)
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "T cell")
        assert result.strategy == "exact"

    def test_normalized_fires_before_initials(self, adata):
        # "nk cell" normalizes to "nk_cell"; "NK cell" normalizes to "nk_cell"
        # normalized match should fire before the initials branch
        _set_cell_type(adata, _STANDARD_LABELS)
        result = lookup_cell_type_name(adata, "nk cell")
        assert result.strategy == "normalized"
        assert result.resolved_name == "NK cell"

    def test_initials_uppercase_fires_before_first_letter(self, adata):
        # Add a label where uppercase initials match unambiguously, before
        # first-letter would also match
        labels = ["Natural Killer", "T cell"]
        # _initials_uppercase("Natural Killer") = "NK"
        # _initials_first_letters("Natural Killer") = "NK" — same result
        # uppercase branch fires first; if unambiguous it returns before first-letter
        _set_cell_type(adata, labels)
        result = lookup_cell_type_name(adata, "NK")
        assert result.matched is True
        assert result.strategy == "initials_uppercase"
        assert result.resolved_name == "Natural Killer"


# ---------------------------------------------------------------------------
# 10. Behavioral edge cases — pin current behavior
# ---------------------------------------------------------------------------


def test_whitespace_only_query_rejected_even_with_whitespace_label(adata):
    """Whitespace-only raw query is rejected by the degenerate-input guard.

    Before the guard, "   " in all_names matched a "   " cell_type label
    and returned matched=True — semantically wrong. The guard at the top of
    `lookup_cell_type_name` now strips raw and returns matched=False when
    empty, regardless of what labels are stored in obs.
    """
    import pandas as pd

    adata.obs["cell_type"] = pd.Categorical(["   "] * adata.n_obs)
    result = lookup_cell_type_name(adata, "   ")
    assert result.matched is False
    assert result.candidates == []
