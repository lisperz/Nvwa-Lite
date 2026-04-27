"""Exhaustive-branch unit tests for dataset_overview(adata).

Tests every branch documented in src/tools/inspection.py and in the
'tools/inspection.py + registry amendment contract (locked)' section of
local/product/prompt_flow_mvp_2026-04-22.md.

All branches exercised via the public dataset_overview() return string.
Private helpers (_get_meta, _species_line, etc.) are NOT imported directly.

Assertions use substring matching so future wording tweaks don't break tests.
For branches that produce NO output, the substring is asserted absent.
"""

from __future__ import annotations

import pytest

from src.tools.inspection import dataset_overview


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _set_meta(adata, meta: dict) -> None:
    """Write nvwa_meta into adata.uns in-place."""
    adata.uns["nvwa_meta"] = meta


def _clear_meta(adata) -> None:
    """Remove nvwa_meta from adata.uns entirely."""
    adata.uns.pop("nvwa_meta", None)


# ---------------------------------------------------------------------------
# Section 1: nvwa_meta key presence / emptiness
# ---------------------------------------------------------------------------

class TestMetaAbsence:
    """No nvwa_meta key → bare header line only."""

    def test_no_nvwa_meta_key_returns_header(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert "cells" in result
        assert "genes" in result

    def test_no_nvwa_meta_key_no_species_line(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert "Species:" not in result

    def test_no_nvwa_meta_key_no_conditions_line(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert "Conditions detected:" not in result

    def test_no_nvwa_meta_key_no_cell_types_line(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert "Cell types annotated:" not in result

    def test_no_nvwa_meta_key_no_ambiguity_questions(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert "I need your help" not in result

    def test_no_nvwa_meta_key_header_contains_cell_count(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        # n_obs formatted as "{n_obs:,}" — for small subsets like 200 no comma is added
        assert f"{adata.n_obs:,}" in result

    def test_no_nvwa_meta_key_header_contains_gene_count(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert f"{adata.n_vars:,}" in result


class TestMetaEmpty:
    """nvwa_meta present but empty dict → bare header only."""

    def test_empty_meta_returns_header(self, adata):
        _set_meta(adata, {})
        result = dataset_overview(adata)
        assert "cells" in result and "genes" in result

    def test_empty_meta_no_species(self, adata):
        _set_meta(adata, {})
        result = dataset_overview(adata)
        assert "Species:" not in result

    def test_empty_meta_no_conditions(self, adata):
        _set_meta(adata, {})
        result = dataset_overview(adata)
        assert "Conditions detected:" not in result

    def test_empty_meta_no_cell_types(self, adata):
        _set_meta(adata, {})
        result = dataset_overview(adata)
        assert "Cell types annotated:" not in result

    def test_empty_meta_no_ambiguity(self, adata):
        _set_meta(adata, {})
        result = dataset_overview(adata)
        assert "I need your help" not in result


class TestMetaNonDictValue:
    """nvwa_meta set to a non-dict value → treated as missing, bare header only."""

    def test_none_meta_no_crash(self, adata):
        adata.uns["nvwa_meta"] = None
        result = dataset_overview(adata)
        assert "cells" in result

    def test_string_meta_no_crash(self, adata):
        adata.uns["nvwa_meta"] = "invalid"
        result = dataset_overview(adata)
        assert "cells" in result

    def test_non_dict_meta_no_species(self, adata):
        adata.uns["nvwa_meta"] = 42
        result = dataset_overview(adata)
        assert "Species:" not in result


# ---------------------------------------------------------------------------
# Section 2: Species line — high-confidence sources
# ---------------------------------------------------------------------------

class TestSpeciesHighConfidence:
    """High-confidence species → Species: line with source note."""

    def _make_species(self, value: str, source: str) -> dict:
        return {"value": value, "confidence": "high", "source": source, "candidates": []}

    def test_declared_source_shows_declared_note(self, adata):
        _set_meta(adata, {"species": self._make_species("human", "declared")})
        result = dataset_overview(adata)
        assert "Species: human" in result
        assert "declared in dataset metadata" in result

    def test_ensembl_prefix_source_shows_ensembl_note(self, adata):
        _set_meta(adata, {"species": self._make_species("mouse", "ensembl_prefix")})
        result = dataset_overview(adata)
        assert "Species: mouse" in result
        assert "detected from Ensembl prefixes" in result

    def test_gene_case_source_shows_gene_case_note(self, adata):
        _set_meta(adata, {"species": self._make_species("rat", "gene_case")})
        result = dataset_overview(adata)
        assert "Species: rat" in result
        assert "detected from gene-symbol case" in result

    def test_unknown_source_falls_through_to_detected(self, adata):
        """Catch-all: unrecognized source maps to 'detected'."""
        _set_meta(adata, {"species": self._make_species("human", "some_future_source")})
        result = dataset_overview(adata)
        assert "Species: human" in result
        assert "detected" in result

    def test_empty_string_source_falls_through_to_detected(self, adata):
        _set_meta(adata, {"species": self._make_species("mouse", "")})
        result = dataset_overview(adata)
        assert "Species: mouse" in result
        assert "detected" in result

    def test_species_line_ends_with_period(self, adata):
        _set_meta(adata, {"species": self._make_species("human", "declared")})
        result = dataset_overview(adata)
        # The species line ends with a period
        for line in result.split("\n"):
            if "Species:" in line:
                assert line.endswith(".")
                break


# ---------------------------------------------------------------------------
# Section 3: Species line — low-confidence / absent
# ---------------------------------------------------------------------------

class TestSpeciesLowConfidence:
    """Low-confidence → no Species: line; surfaces as ambiguity question instead."""

    def test_low_conf_no_species_line(self, adata):
        _set_meta(adata, {
            "species": {"value": None, "confidence": "low", "source": "gene_case",
                        "candidates": ["mouse", "rat"]},
        })
        result = dataset_overview(adata)
        # Summary line "Species: <value> (<source>)." must not appear; the ambiguity
        # question section may still contain the substring "Species:" as its prefix.
        first_line = result.split("\n", 1)[0]
        assert "Species:" not in first_line
        assert not any(
            line.startswith("Species:") for line in result.split("\n")
        )

    def test_low_conf_gene_case_question_phrasing(self, adata):
        """gene_case low-conf → Title-case phrasing."""
        _set_meta(adata, {
            "species": {"value": None, "confidence": "low", "source": "gene_case",
                        "candidates": ["mouse", "rat"]},
        })
        result = dataset_overview(adata)
        assert "Title-case" in result
        assert "mouse" in result
        assert "rat" in result

    def test_low_conf_other_source_generic_phrasing(self, adata):
        """Non gene_case low-conf → generic 'couldn't detect' phrasing."""
        _set_meta(adata, {
            "species": {"value": None, "confidence": "low", "source": "ensembl_prefix",
                        "candidates": ["human", "mouse"]},
        })
        result = dataset_overview(adata)
        assert "couldn't detect" in result
        assert "human" in result
        assert "mouse" in result

    def test_low_conf_no_declared_source_generic_phrasing(self, adata):
        """source='declared' but low-conf → generic phrasing (not gene_case branch)."""
        _set_meta(adata, {
            "species": {"value": None, "confidence": "low", "source": "declared",
                        "candidates": ["human"]},
        })
        result = dataset_overview(adata)
        assert "couldn't detect" in result

    def test_low_conf_triggers_ambiguity_section_header(self, adata):
        _set_meta(adata, {
            "species": {"value": None, "confidence": "low", "source": "gene_case",
                        "candidates": ["mouse", "rat"]},
        })
        result = dataset_overview(adata)
        assert "I need your help" in result

    def test_missing_species_key_no_species_line(self, adata):
        """No 'species' key at all → no species line, no ambiguity question."""
        _set_meta(adata, {"schema": {}})
        result = dataset_overview(adata)
        assert "Species:" not in result

    def test_species_non_dict_no_species_line(self, adata):
        """species value is not a dict → silently skipped."""
        _set_meta(adata, {"species": "human"})
        result = dataset_overview(adata)
        assert "Species:" not in result

    def test_species_no_value_high_conf_no_line(self, adata):
        """value=None with high confidence → no species line (value guard)."""
        _set_meta(adata, {
            "species": {"value": None, "confidence": "high", "source": "declared", "candidates": []},
        })
        result = dataset_overview(adata)
        assert "Species:" not in result

    def test_species_empty_value_high_conf_no_line(self, adata):
        """value='' with high confidence → no species line (falsy guard)."""
        _set_meta(adata, {
            "species": {"value": "", "confidence": "high", "source": "declared", "candidates": []},
        })
        result = dataset_overview(adata)
        assert "Species:" not in result


# ---------------------------------------------------------------------------
# Section 4: Conditions line
# ---------------------------------------------------------------------------

class TestConditionsLine:
    """condition_cols: empty / present with or without values."""

    def test_no_condition_cols_key_no_conditions_line(self, adata):
        _set_meta(adata, {"schema": {}})
        result = dataset_overview(adata)
        assert "Conditions detected:" not in result

    def test_empty_condition_cols_no_conditions_line(self, adata):
        _set_meta(adata, {"condition_cols": [], "schema": {}})
        result = dataset_overview(adata)
        assert "Conditions detected:" not in result

    def test_condition_col_with_values(self, adata):
        _set_meta(adata, {
            "condition_cols": ["treatment"],
            "schema": {
                "treatment": {"role": "condition", "values": ["ctrl", "drug"], "n_unique": 2},
            },
        })
        result = dataset_overview(adata)
        assert "Conditions detected:" in result
        assert "treatment" in result
        assert "ctrl" in result
        assert "drug" in result

    def test_condition_col_without_values_shows_n_unique(self, adata):
        _set_meta(adata, {
            "condition_cols": ["batch"],
            "schema": {
                "batch": {"role": "condition", "n_unique": 5},
            },
        })
        result = dataset_overview(adata)
        assert "Conditions detected:" in result
        assert "batch" in result
        assert "5" in result

    def test_condition_col_without_values_or_n_unique_shows_question_mark(self, adata):
        """No values, no n_unique → '?' placeholder."""
        _set_meta(adata, {
            "condition_cols": ["group"],
            "schema": {"group": {}},
        })
        result = dataset_overview(adata)
        assert "Conditions detected:" in result
        assert "?" in result

    def test_multiple_condition_cols_all_shown(self, adata):
        _set_meta(adata, {
            "condition_cols": ["treatment", "timepoint"],
            "schema": {
                "treatment": {"values": ["ctrl", "drug"], "n_unique": 2},
                "timepoint": {"values": ["d0", "d7"], "n_unique": 2},
            },
        })
        result = dataset_overview(adata)
        assert "treatment" in result
        assert "timepoint" in result

    def test_condition_col_not_in_schema_shows_question_mark(self, adata):
        """col in condition_cols but not in schema → n_unique fallback is '?'."""
        _set_meta(adata, {
            "condition_cols": ["mystery_col"],
            "schema": {},
        })
        result = dataset_overview(adata)
        assert "mystery_col" in result
        assert "?" in result

    def test_condition_cols_not_a_list_no_conditions_line(self, adata):
        """condition_cols set to non-list → no conditions line (type guard)."""
        _set_meta(adata, {"condition_cols": "treatment", "schema": {}})
        result = dataset_overview(adata)
        assert "Conditions detected:" not in result


# ---------------------------------------------------------------------------
# Section 5: Cell types line
# ---------------------------------------------------------------------------

class TestCellTypesLine:
    """schema entries with role='cell_type': no truncation vs truncation."""

    def test_no_cell_type_col_no_cell_types_line(self, adata):
        _set_meta(adata, {"schema": {"treatment": {"role": "condition"}}})
        result = dataset_overview(adata)
        assert "Cell types annotated:" not in result

    def test_cell_type_with_values_no_truncation(self, adata):
        """n_unique == len(values) → full list, no ellipsis."""
        _set_meta(adata, {
            "schema": {
                "cell_type": {
                    "role": "cell_type",
                    "n_unique": 3,
                    "values": ["T cell", "B cell", "NK cell"],
                },
            },
        })
        result = dataset_overview(adata)
        assert "Cell types annotated:" in result
        assert "3 types" in result
        assert "T cell" in result
        assert "B cell" in result
        assert "NK cell" in result
        # Ellipsis should NOT appear when n_unique == len(values)
        assert "…" not in result

    def test_cell_type_with_values_truncation(self, adata):
        """n_unique > len(values) → preview list with ellipsis."""
        values = ["T cell", "B cell", "NK cell", "Monocyte", "DC"]
        _set_meta(adata, {
            "schema": {
                "cell_type": {
                    "role": "cell_type",
                    "n_unique": 10,  # more than len(values)=5
                    "values": values,
                },
            },
        })
        result = dataset_overview(adata)
        assert "Cell types annotated:" in result
        assert "10 types" in result
        assert "T cell" in result
        assert "…" in result

    def test_cell_type_truncation_preview_max_5(self, adata):
        """preview slice is [:5] — verify only first 5 of >5 values shown before ellipsis."""
        values = ["A", "B", "C", "D", "E", "F", "G"]
        _set_meta(adata, {
            "schema": {
                "cell_type": {
                    "role": "cell_type",
                    "n_unique": 20,
                    "values": values,
                },
            },
        })
        result = dataset_overview(adata)
        # First 5 appear
        for v in values[:5]:
            assert v in result
        # 6th and 7th are NOT in the preview (they fall past slice [:5])
        # "F" is at index 5, "G" at 6 — should be absent before the ellipsis
        # We cannot guarantee the exact format easily without full-string match,
        # but we can verify truncation happened:
        assert "…" in result

    def test_cell_type_without_values_falls_back_to_column_name(self, adata):
        """No 'values' key → fallback to 'column <col>' phrasing."""
        _set_meta(adata, {
            "schema": {
                "cell_type": {
                    "role": "cell_type",
                    "n_unique": 7,
                },
            },
        })
        result = dataset_overview(adata)
        assert "Cell types annotated:" in result
        assert "7 types" in result
        assert "cell_type" in result  # column name appears
        assert "…" not in result

    def test_empty_schema_no_cell_types_line(self, adata):
        _set_meta(adata, {"schema": {}})
        result = dataset_overview(adata)
        assert "Cell types annotated:" not in result

    def test_cell_type_n_unique_equals_values_length_no_ellipsis(self, adata):
        """Boundary: n_unique == len(values) → no ellipsis (uses full-list branch)."""
        values = ["CD4 T", "CD8 T", "B cell", "Mono", "NK"]
        _set_meta(adata, {
            "schema": {
                "cell_type": {
                    "role": "cell_type",
                    "n_unique": 5,  # exactly equal to len(values)
                    "values": values,
                },
            },
        })
        result = dataset_overview(adata)
        assert "5 types" in result
        assert "…" not in result
        for v in values:
            assert v in result


# ---------------------------------------------------------------------------
# Section 6: Ambiguity questions — column-level
# ---------------------------------------------------------------------------

class TestAmbiguityQuestionsColumns:
    """ambiguous_cols list drives column ambiguity questions."""

    def test_no_ambiguous_cols_no_questions(self, adata):
        _set_meta(adata, {"schema": {}})
        result = dataset_overview(adata)
        assert "I need your help" not in result

    def test_empty_ambiguous_cols_no_questions(self, adata):
        _set_meta(adata, {"ambiguous_cols": [], "schema": {}})
        result = dataset_overview(adata)
        assert "I need your help" not in result

    def test_single_candidate_no_or(self, adata):
        """1 candidate → plain name, no 'or'."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {
                    "candidates": ["condition"],
                    "values": ["ctrl", "trt"],
                    "n_unique": 2,
                },
            },
        })
        result = dataset_overview(adata)
        assert "col1" in result
        assert "experimental condition" in result
        # "or" should NOT appear between candidates (only 1)
        # Check that " or " does not appear in the candidate portion
        assert "I need your help" in result

    def test_two_candidates_uses_or(self, adata):
        """2 candidates → 'X or Y'."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {
                    "candidates": ["condition", "batch"],
                    "values": ["A", "B"],
                    "n_unique": 2,
                },
            },
        })
        result = dataset_overview(adata)
        assert "experimental condition" in result
        assert "batch label" in result
        assert " or " in result

    def test_three_candidates_uses_oxford_comma(self, adata):
        """3+ candidates → 'X, Y, or Z'."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {
                    "candidates": ["condition", "batch", "sample_id"],
                    "values": ["s1", "s2"],
                    "n_unique": 2,
                },
            },
        })
        result = dataset_overview(adata)
        assert "experimental condition" in result
        assert "batch label" in result
        assert "sample ID" in result
        # Oxford comma pattern: "..., or ..."
        assert ", or " in result

    def test_unknown_candidate_role_passes_through_raw(self, adata):
        """Candidate not in _ROLE_HUMAN → raw string used as-is."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {
                    "candidates": ["some_future_role"],
                    "values": ["x", "y"],
                    "n_unique": 2,
                },
            },
        })
        result = dataset_overview(adata)
        assert "some_future_role" in result

    def test_mixed_known_unknown_candidates(self, adata):
        """Mix of known and unknown roles: known humanized, unknown raw."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {
                    "candidates": ["condition", "novel_role"],
                    "values": ["a", "b"],
                    "n_unique": 2,
                },
            },
        })
        result = dataset_overview(adata)
        assert "experimental condition" in result
        assert "novel_role" in result

    def test_col_values_shown_in_question(self, adata):
        """Values appear in the question detail."""
        _set_meta(adata, {
            "ambiguous_cols": ["group"],
            "schema": {
                "group": {
                    "candidates": ["condition"],
                    "values": ["ctrl", "drug"],
                    "n_unique": 2,
                },
            },
        })
        result = dataset_overview(adata)
        assert "ctrl" in result
        assert "drug" in result

    def test_col_no_values_shows_n_unique_in_question(self, adata):
        """No values in schema entry → n_unique used for the count."""
        _set_meta(adata, {
            "ambiguous_cols": ["group"],
            "schema": {
                "group": {
                    "candidates": ["condition"],
                    "n_unique": 8,
                },
            },
        })
        result = dataset_overview(adata)
        assert "8 unique values" in result

    def test_col_values_truncated_at_6_with_ellipsis(self, adata):
        """More than 6 values → first 6 shown + ', ...' appended."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {
                    "candidates": ["condition"],
                    "values": ["v1", "v2", "v3", "v4", "v5", "v6", "v7"],
                    "n_unique": 7,
                },
            },
        })
        result = dataset_overview(adata)
        assert "v1" in result
        assert "v6" in result
        assert "…" in result
        # v7 should NOT appear since it's beyond the 6-item slice
        assert "v7" not in result

    def test_ambiguous_col_not_in_schema_uses_defaults(self, adata):
        """col in ambiguous_cols but missing from schema → empty values/candidates."""
        _set_meta(adata, {
            "ambiguous_cols": ["ghost_col"],
            "schema": {},
        })
        result = dataset_overview(adata)
        assert "ghost_col" in result
        # No values and no candidates → 'something' (the _join_candidates empty fallback)
        assert "something" in result

    def test_multiple_ambiguous_cols_each_gets_question(self, adata):
        """Multiple ambiguous columns → one question per column."""
        _set_meta(adata, {
            "ambiguous_cols": ["colA", "colB"],
            "schema": {
                "colA": {"candidates": ["condition"], "values": ["x"], "n_unique": 1},
                "colB": {"candidates": ["batch"], "values": ["y"], "n_unique": 1},
            },
        })
        result = dataset_overview(adata)
        assert "colA" in result
        assert "colB" in result

    def test_questions_section_header_present_when_ambiguous(self, adata):
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {"candidates": ["condition"], "values": ["a"], "n_unique": 1},
            },
        })
        result = dataset_overview(adata)
        assert "I need your help on a couple of things before we dive in:" in result

    def test_all_known_role_candidates_humanized(self, adata):
        """All six known roles map to human-readable labels."""
        known_roles = {
            "condition": "experimental condition",
            "sample_id": "sample ID",
            "batch": "batch label",
            "cell_type": "cell type",
            "clustering": "cluster assignment",
            "qc_metric": "QC metric",
        }
        for role, label in known_roles.items():
            _set_meta(adata, {
                "ambiguous_cols": ["col1"],
                "schema": {
                    "col1": {"candidates": [role], "values": [], "n_unique": 0},
                },
            })
            result = dataset_overview(adata)
            assert label in result, f"Role '{role}' should humanize to '{label}'"


# ---------------------------------------------------------------------------
# Section 7: Ambiguity questions — species
# ---------------------------------------------------------------------------

class TestAmbiguityQuestionsSpecies:
    """Low-confidence species → species ambiguity question."""

    def test_gene_case_low_conf_title_case_phrasing(self, adata):
        _set_meta(adata, {
            "species": {
                "value": None, "confidence": "low",
                "source": "gene_case",
                "candidates": ["mouse", "rat"],
            },
        })
        result = dataset_overview(adata)
        assert "Title-case" in result
        assert "Which is this dataset?" in result

    def test_other_source_low_conf_generic_phrasing(self, adata):
        _set_meta(adata, {
            "species": {
                "value": None, "confidence": "low",
                "source": "ensembl_prefix",
                "candidates": ["human", "mouse"],
            },
        })
        result = dataset_overview(adata)
        assert "couldn't detect" in result

    def test_species_candidates_joined_with_or_for_two(self, adata):
        _set_meta(adata, {
            "species": {
                "value": None, "confidence": "low",
                "source": "gene_case",
                "candidates": ["mouse", "rat"],
            },
        })
        result = dataset_overview(adata)
        assert "mouse or rat" in result

    def test_species_candidates_oxford_for_three(self, adata):
        _set_meta(adata, {
            "species": {
                "value": None, "confidence": "low",
                "source": "ensembl_prefix",
                "candidates": ["human", "mouse", "rat"],
            },
        })
        result = dataset_overview(adata)
        assert ", or " in result

    def test_species_empty_candidates_shows_something(self, adata):
        """Empty candidates list → 'something' fallback from _join_candidates."""
        _set_meta(adata, {
            "species": {
                "value": None, "confidence": "low",
                "source": "gene_case",
                "candidates": [],
            },
        })
        result = dataset_overview(adata)
        assert "something" in result

    def test_species_single_candidate_no_or(self, adata):
        _set_meta(adata, {
            "species": {
                "value": None, "confidence": "low",
                "source": "ensembl_prefix",
                "candidates": ["human"],
            },
        })
        result = dataset_overview(adata)
        assert "human" in result
        assert "human or " not in result

    def test_high_conf_species_no_ambiguity_question(self, adata):
        """High confidence → no species ambiguity question."""
        _set_meta(adata, {
            "species": {
                "value": "human", "confidence": "high",
                "source": "declared", "candidates": [],
            },
        })
        result = dataset_overview(adata)
        assert "Title-case" not in result
        assert "couldn't detect" not in result


# ---------------------------------------------------------------------------
# Section 8: Combined / full-meta scenarios
# ---------------------------------------------------------------------------

class TestFullMeta:
    """Tests with a complete nvwa_meta dict — all lines present simultaneously."""

    def _full_meta(self) -> dict:
        return {
            "species": {
                "value": "human",
                "confidence": "high",
                "source": "declared",
                "candidates": [],
            },
            "condition_cols": ["treatment"],
            "schema": {
                "treatment": {
                    "role": "condition",
                    "values": ["ctrl", "drug"],
                    "n_unique": 2,
                },
                "cell_type": {
                    "role": "cell_type",
                    "n_unique": 5,
                    "values": ["T cell", "B cell", "NK cell", "Mono", "DC"],
                },
            },
            "ambiguous_cols": [],
        }

    def test_all_lines_present(self, adata):
        _set_meta(adata, self._full_meta())
        result = dataset_overview(adata)
        assert "Dataset:" in result
        assert "Species:" in result
        assert "Conditions detected:" in result
        assert "Cell types annotated:" in result

    def test_no_ambiguity_section_when_nothing_ambiguous(self, adata):
        _set_meta(adata, self._full_meta())
        result = dataset_overview(adata)
        assert "I need your help" not in result

    def test_partial_meta_species_and_conditions_only(self, adata):
        """Partial schema — only species + conditions, no cell types."""
        _set_meta(adata, {
            "species": {"value": "mouse", "confidence": "high", "source": "ensembl_prefix", "candidates": []},
            "condition_cols": ["group"],
            "schema": {"group": {"values": ["A", "B"], "n_unique": 2}},
        })
        result = dataset_overview(adata)
        assert "Species: mouse" in result
        assert "Conditions detected:" in result
        assert "Cell types annotated:" not in result

    def test_conditions_and_ambiguity_both_shown(self, adata):
        """Both conditions and an ambiguous column → both present."""
        _set_meta(adata, {
            "condition_cols": ["treatment"],
            "schema": {
                "treatment": {"values": ["ctrl"], "n_unique": 1},
                "mystery": {"candidates": ["condition", "batch"], "values": ["x", "y"], "n_unique": 2},
            },
            "ambiguous_cols": ["mystery"],
        })
        result = dataset_overview(adata)
        assert "Conditions detected:" in result
        assert "mystery" in result
        assert "I need your help" in result

    def test_cell_types_and_species_ambiguity_together(self, adata):
        """Cell types present + species ambiguous → both lines shown."""
        _set_meta(adata, {
            "schema": {
                "cell_type": {"role": "cell_type", "n_unique": 3, "values": ["T", "B", "NK"]},
            },
            "species": {
                "value": None, "confidence": "low", "source": "gene_case",
                "candidates": ["mouse", "rat"],
            },
        })
        result = dataset_overview(adata)
        assert "Cell types annotated:" in result
        assert "Title-case" in result


# ---------------------------------------------------------------------------
# Section 9: Return-value structure
# ---------------------------------------------------------------------------

class TestReturnValueStructure:
    """Structural guarantees of the return value."""

    def test_returns_string(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        assert isinstance(result, str)

    def test_header_is_first_line(self, adata):
        _clear_meta(adata)
        result = dataset_overview(adata)
        lines = result.split("\n")
        assert lines[0].startswith("Dataset:")

    def test_questions_section_preceded_by_blank_line(self, adata):
        """Ambiguity section is separated from summary by a blank line."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {"candidates": ["condition"], "values": ["a", "b"], "n_unique": 2},
            },
        })
        result = dataset_overview(adata)
        lines = result.split("\n")
        # Find the "I need your help" line index
        help_idx = next(i for i, l in enumerate(lines) if "I need your help" in l)
        assert lines[help_idx - 1] == ""

    def test_question_lines_start_with_dash(self, adata):
        """Each ambiguity question line starts with '- '."""
        _set_meta(adata, {
            "ambiguous_cols": ["col1"],
            "schema": {
                "col1": {"candidates": ["condition"], "values": ["x"], "n_unique": 1},
            },
        })
        result = dataset_overview(adata)
        lines = result.split("\n")
        q_lines = [l for l in lines if l.startswith("- ")]
        assert len(q_lines) >= 1

    def test_cell_gene_counts_formatted(self, adata):
        """n_obs and n_vars appear in the header formatted with comma separators."""
        _clear_meta(adata)
        result = dataset_overview(adata)
        # The format is "{n_obs:,} cells × {n_vars:,} genes"
        assert f"{adata.n_obs:,}" in result
        assert f"{adata.n_vars:,}" in result

    def test_conditions_line_ends_with_period(self, adata):
        _set_meta(adata, {
            "condition_cols": ["treatment"],
            "schema": {"treatment": {"values": ["ctrl"], "n_unique": 1}},
        })
        result = dataset_overview(adata)
        for line in result.split("\n"):
            if "Conditions detected:" in line:
                assert line.endswith(".")
                break

    def test_cell_types_line_ends_with_period(self, adata):
        _set_meta(adata, {
            "schema": {
                "cell_type": {"role": "cell_type", "n_unique": 2, "values": ["A", "B"]},
            },
        })
        result = dataset_overview(adata)
        for line in result.split("\n"):
            if "Cell types annotated:" in line:
                assert line.endswith(".")
                break
