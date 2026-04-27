"""Exhaustive-branch unit tests for lookup_condition_name().

Covers every branch in the domain/resolver/condition_lookup.py contract:
  - exact match in a single declared column
  - match in one of multiple declared columns
  - case-insensitive (normalized) match
  - raw not found in any declared column -> matched=False
  - ambiguous: same label in two different declared columns -> candidates=[...]
  - ambiguous: normalized match in two different declared columns -> candidates=[...]
  - no nvwa_meta key at all -> graceful fallback (defaults to ["condition"])
  - nvwa_meta present but no condition_cols key -> graceful fallback
  - condition_cols is an empty list -> not found
  - declared column does not exist in adata.obs (stale declaration) -> no crash

The shared `adata` fixture comes from conftest.py (per-test copy, mutable).
"""

from __future__ import annotations

import pytest

from src.domain.resolver.condition_lookup import (
    ConditionLookupResult,
    lookup_condition_name,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _set_condition_cols(adata, cols: list[str]) -> None:
    """Write nvwa_meta.condition_cols onto adata.uns in place."""
    if "nvwa_meta" not in adata.uns:
        adata.uns["nvwa_meta"] = {}
    adata.uns["nvwa_meta"]["condition_cols"] = cols


# ---------------------------------------------------------------------------
# 1. Exact match — single declared column
# ---------------------------------------------------------------------------

class TestExactMatchSingleColumn:
    """Raw value found verbatim in the one declared column."""

    def test_returns_matched_true(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is True

    def test_resolved_name_equals_raw(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.resolved_name == "WT"

    def test_obs_column_set(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.obs_column == "treatment"

    def test_strategy_is_exact(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.strategy == "exact"

    def test_input_name_preserved(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.input_name == "WT"

    def test_candidates_empty(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.candidates == []

    def test_message_non_empty(self, adata):
        adata.obs["treatment"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "WT")
        assert result.message != ""

    def test_multiple_unique_values_in_column(self, adata):
        """Column has several values; only one matches -> still exact."""
        import pandas as pd
        vals = (["WT"] * (adata.n_obs // 2)) + (["KO"] * (adata.n_obs - adata.n_obs // 2))
        adata.obs["treatment"] = pd.Categorical(vals)
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "KO")
        assert result.matched is True
        assert result.resolved_name == "KO"
        assert result.obs_column == "treatment"


# ---------------------------------------------------------------------------
# 2. Exact match — one of multiple declared columns
# ---------------------------------------------------------------------------

class TestExactMatchMultipleColumns:
    """Raw value present in one column; absent from sibling columns."""

    def test_matched_true(self, adata):
        import pandas as pd
        adata.obs["genotype"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["timepoint"] = pd.Categorical(["D7"] * adata.n_obs)
        _set_condition_cols(adata, ["genotype", "timepoint"])
        result = lookup_condition_name(adata, "D7")
        assert result.matched is True

    def test_correct_column_identified(self, adata):
        import pandas as pd
        adata.obs["genotype"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["timepoint"] = pd.Categorical(["D7"] * adata.n_obs)
        _set_condition_cols(adata, ["genotype", "timepoint"])
        result = lookup_condition_name(adata, "D7")
        assert result.obs_column == "timepoint"

    def test_resolved_name_correct(self, adata):
        import pandas as pd
        adata.obs["genotype"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["timepoint"] = pd.Categorical(["D7"] * adata.n_obs)
        _set_condition_cols(adata, ["genotype", "timepoint"])
        result = lookup_condition_name(adata, "D7")
        assert result.resolved_name == "D7"

    def test_strategy_is_exact(self, adata):
        import pandas as pd
        adata.obs["genotype"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["timepoint"] = pd.Categorical(["D7"] * adata.n_obs)
        _set_condition_cols(adata, ["genotype", "timepoint"])
        result = lookup_condition_name(adata, "D7")
        assert result.strategy == "exact"

    def test_three_columns_match_in_last(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["X"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["Y"] * adata.n_obs)
        adata.obs["col_c"] = pd.Categorical(["TARGET"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b", "col_c"])
        result = lookup_condition_name(adata, "TARGET")
        assert result.matched is True
        assert result.obs_column == "col_c"


# ---------------------------------------------------------------------------
# 3. Case-insensitive (normalized) match
# ---------------------------------------------------------------------------

class TestNormalizedMatch:
    """Raw value differs only by case / whitespace / separator from stored value."""

    def test_lowercase_raw_matches_uppercase_stored(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "wt")
        assert result.matched is True

    def test_strategy_is_normalized(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "wt")
        assert result.strategy == "normalized"

    def test_resolved_name_is_stored_form(self, adata):
        """resolved_name should be the canonical stored form, not the raw input."""
        adata.obs["condition"] = ["Wild Type"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "wild_type")
        assert result.resolved_name == "Wild Type"

    def test_hyphen_in_raw_matches_space_in_stored(self, adata):
        """Normalization collapses hyphens and spaces to underscore."""
        adata.obs["condition"] = ["Wild Type"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "wild-type")
        assert result.matched is True

    def test_obs_column_set_on_normalized_match(self, adata):
        adata.obs["treatment"] = ["Knock Out"] * adata.n_obs
        _set_condition_cols(adata, ["treatment"])
        result = lookup_condition_name(adata, "knock_out")
        assert result.obs_column == "treatment"

    def test_normalized_match_in_one_of_multiple_columns(self, adata):
        import pandas as pd
        adata.obs["genotype"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["phase"] = pd.Categorical(["G1 Phase"] * adata.n_obs)
        _set_condition_cols(adata, ["genotype", "phase"])
        result = lookup_condition_name(adata, "g1_phase")
        assert result.matched is True
        assert result.obs_column == "phase"


# ---------------------------------------------------------------------------
# 4. Not found in any declared column
# ---------------------------------------------------------------------------

class TestNotFound:
    """Raw value absent from all declared columns."""

    def test_matched_false(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert result.matched is False

    def test_resolved_name_is_none(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert result.resolved_name is None

    def test_obs_column_is_none(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert result.obs_column is None

    def test_strategy_is_none(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert result.strategy is None

    def test_candidates_empty(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert result.candidates == []

    def test_message_mentions_raw(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert "NONEXISTENT" in result.message

    def test_input_name_preserved(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "NONEXISTENT")
        assert result.input_name == "NONEXISTENT"

    def test_not_found_in_multiple_columns(self, adata):
        import pandas as pd
        adata.obs["genotype"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["timepoint"] = pd.Categorical(["D7"] * adata.n_obs)
        _set_condition_cols(adata, ["genotype", "timepoint"])
        result = lookup_condition_name(adata, "MISSING")
        assert result.matched is False


# ---------------------------------------------------------------------------
# 5. Ambiguous — exact match in two or more declared columns
# ---------------------------------------------------------------------------

class TestAmbiguousExactMatch:
    """Same label appears in two different declared columns -> ambiguous."""

    def test_matched_false(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is False

    def test_candidates_populated(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert len(result.candidates) == 2

    def test_candidates_format_column_colon_value(self, adata):
        """Each candidate should be 'column:value'."""
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert "col_a:WT" in result.candidates
        assert "col_b:WT" in result.candidates

    def test_resolved_name_is_none(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert result.resolved_name is None

    def test_obs_column_is_none(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert result.obs_column is None

    def test_strategy_is_exact(self, adata):
        """Ambiguous exact -> strategy should be 'exact' (not normalized)."""
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert result.strategy == "exact"

    def test_three_columns_all_ambiguous(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_c"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b", "col_c"])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is False
        assert len(result.candidates) == 3


# ---------------------------------------------------------------------------
# 6. Ambiguous — normalized match in two or more declared columns
# ---------------------------------------------------------------------------

class TestAmbiguousNormalizedMatch:
    """Different stored forms that normalize identically across two columns."""

    def test_matched_false(self, adata):
        import pandas as pd
        # "Wild Type" and "wild_type" both normalize to "wild_type"
        adata.obs["col_a"] = pd.Categorical(["Wild Type"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["wild_type"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WILD-TYPE")
        assert result.matched is False

    def test_candidates_populated(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["Wild Type"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["wild_type"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WILD-TYPE")
        assert len(result.candidates) == 2

    def test_candidates_format(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["Wild Type"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["wild_type"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WILD-TYPE")
        assert "col_a:Wild Type" in result.candidates
        assert "col_b:wild_type" in result.candidates

    def test_strategy_is_normalized(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["Wild Type"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["wild_type"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WILD-TYPE")
        assert result.strategy == "normalized"


# ---------------------------------------------------------------------------
# 7. No nvwa_meta key at all -> graceful fallback to ["condition"]
# ---------------------------------------------------------------------------

class TestNoNvwaMeta:
    """adata.uns has no nvwa_meta -> _get_condition_cols falls back to ['condition']."""

    def test_no_nvwa_meta_and_condition_col_present(self, adata):
        """Falls back to 'condition'; exact match succeeds if that col exists."""
        adata.uns.pop("nvwa_meta", None)
        adata.obs["condition"] = ["WT"] * adata.n_obs
        result = lookup_condition_name(adata, "WT")
        assert result.matched is True
        assert result.obs_column == "condition"

    def test_no_nvwa_meta_and_condition_col_absent(self, adata):
        """Falls back to 'condition'; stale col -> not found, no crash."""
        adata.uns.pop("nvwa_meta", None)
        # Ensure 'condition' col is not in obs
        if "condition" in adata.obs.columns:
            del adata.obs["condition"]
        result = lookup_condition_name(adata, "WT")
        assert result.matched is False

    def test_no_nvwa_meta_does_not_raise(self, adata):
        adata.uns.pop("nvwa_meta", None)
        try:
            lookup_condition_name(adata, "anything")
        except Exception as exc:
            pytest.fail(f"Should not raise; got {exc}")

    def test_no_nvwa_meta_returns_condition_lookup_result(self, adata):
        adata.uns.pop("nvwa_meta", None)
        result = lookup_condition_name(adata, "WT")
        assert isinstance(result, ConditionLookupResult)


# ---------------------------------------------------------------------------
# 8. nvwa_meta present but no condition_cols key -> graceful fallback
# ---------------------------------------------------------------------------

class TestNvwaMetaNoConditionCols:
    """nvwa_meta exists but has no condition_cols key -> falls back to ['condition']."""

    def test_fallback_to_condition_col(self, adata):
        adata.uns["nvwa_meta"] = {"species": "human"}
        adata.obs["condition"] = ["KO"] * adata.n_obs
        result = lookup_condition_name(adata, "KO")
        assert result.matched is True
        assert result.obs_column == "condition"

    def test_does_not_raise(self, adata):
        adata.uns["nvwa_meta"] = {}
        try:
            lookup_condition_name(adata, "anything")
        except Exception as exc:
            pytest.fail(f"Should not raise; got {exc}")

    def test_returns_correct_type(self, adata):
        adata.uns["nvwa_meta"] = {}
        result = lookup_condition_name(adata, "anything")
        assert isinstance(result, ConditionLookupResult)


# ---------------------------------------------------------------------------
# 9. condition_cols is an empty list -> not found
# ---------------------------------------------------------------------------

class TestEmptyConditionCols:
    """condition_cols = [] means no columns to scan; every lookup should miss."""

    def test_matched_false(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, [])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is False

    def test_does_not_raise(self, adata):
        _set_condition_cols(adata, [])
        try:
            lookup_condition_name(adata, "WT")
        except Exception as exc:
            pytest.fail(f"Should not raise; got {exc}")

    def test_candidates_empty(self, adata):
        _set_condition_cols(adata, [])
        result = lookup_condition_name(adata, "WT")
        assert result.candidates == []

    def test_input_name_preserved(self, adata):
        _set_condition_cols(adata, [])
        result = lookup_condition_name(adata, "WT")
        assert result.input_name == "WT"


# ---------------------------------------------------------------------------
# 10. Declared column does not exist in adata.obs (stale declaration)
# ---------------------------------------------------------------------------

class TestStaleColumnDeclaration:
    """Column listed in condition_cols but absent from adata.obs -> no crash."""

    def test_stale_column_does_not_raise(self, adata):
        _set_condition_cols(adata, ["nonexistent_col"])
        try:
            lookup_condition_name(adata, "WT")
        except Exception as exc:
            pytest.fail(f"Should not raise on stale column; got {exc}")

    def test_stale_column_returns_not_found(self, adata):
        _set_condition_cols(adata, ["nonexistent_col"])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is False

    def test_stale_col_plus_valid_col_still_matches(self, adata):
        """Stale entry is skipped; valid sibling column still resolves."""
        adata.obs["real_col"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["nonexistent_col", "real_col"])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is True
        assert result.obs_column == "real_col"

    def test_all_stale_returns_not_found(self, adata):
        """All declared columns are stale -> not found, no crash."""
        _set_condition_cols(adata, ["ghost_a", "ghost_b"])
        result = lookup_condition_name(adata, "WT")
        assert result.matched is False
        assert result.candidates == []

    def test_all_stale_returns_result_type(self, adata):
        _set_condition_cols(adata, ["ghost_a", "ghost_b"])
        result = lookup_condition_name(adata, "anything")
        assert isinstance(result, ConditionLookupResult)


# ---------------------------------------------------------------------------
# 11. Return type and field completeness
# ---------------------------------------------------------------------------

class TestReturnType:
    """Every code path returns a well-formed ConditionLookupResult."""

    def test_exact_match_result_type(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "WT")
        assert isinstance(result, ConditionLookupResult)

    def test_not_found_result_type(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "MISSING")
        assert isinstance(result, ConditionLookupResult)

    def test_ambiguous_result_type(self, adata):
        import pandas as pd
        adata.obs["col_a"] = pd.Categorical(["WT"] * adata.n_obs)
        adata.obs["col_b"] = pd.Categorical(["WT"] * adata.n_obs)
        _set_condition_cols(adata, ["col_a", "col_b"])
        result = lookup_condition_name(adata, "WT")
        assert isinstance(result, ConditionLookupResult)

    def test_all_fields_present_on_match(self, adata):
        adata.obs["condition"] = ["WT"] * adata.n_obs
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "WT")
        assert hasattr(result, "matched")
        assert hasattr(result, "input_name")
        assert hasattr(result, "resolved_name")
        assert hasattr(result, "obs_column")
        assert hasattr(result, "strategy")
        assert hasattr(result, "candidates")
        assert hasattr(result, "message")

    def test_all_fields_present_on_miss(self, adata):
        _set_condition_cols(adata, ["condition"])
        result = lookup_condition_name(adata, "MISSING")
        assert hasattr(result, "matched")
        assert hasattr(result, "input_name")
        assert hasattr(result, "resolved_name")
        assert hasattr(result, "obs_column")
        assert hasattr(result, "strategy")
        assert hasattr(result, "candidates")
        assert hasattr(result, "message")
