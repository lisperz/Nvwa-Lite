"""Unit tests for src/domain/resolver/column_classifier.py — classify() public API.

Covers every branch documented in the
  domain/resolver/column_classifier.py contract (locked)
section of local/product/prompt_flow_mvp_2026-04-22.md.

Run:
    python -m pytest tests/unit/spec_pipeline/test_column_classifier.py -v
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.domain.resolver.column_classifier import Role, classify

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_LOW = 20   # _LOW_CARDINALITY_THRESHOLD (from source)
_PREVIEW_CAP = 20  # _MAX_VALUES_PREVIEW (from source)


def _nvwa(adata):
    """Shorthand to get adata.uns['nvwa_meta']."""
    return adata.uns["nvwa_meta"]


def _schema(adata, col):
    """Shorthand to get the schema entry dict for a single column."""
    return _nvwa(adata)["schema"][col]


# ---------------------------------------------------------------------------
# 1. All 8 Role values produced — one representative column each
# ---------------------------------------------------------------------------

class TestAllRolesProduced:
    """Each of the 8 Role values must appear at least once after classify()."""

    def test_role_condition_via_known_name(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "treatment")["role"] == Role.CONDITION.value

    def test_role_cell_type_via_known_name(self, adata):
        adata.obs["cell_type"] = pd.Categorical(["T cell", "B cell"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "cell_type")["role"] == Role.CELL_TYPE.value

    def test_role_clustering_via_known_name(self, adata):
        adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "leiden")["role"] == Role.CLUSTERING.value

    def test_role_qc_metric_via_known_pattern(self, adata):
        adata.obs["nCount_RNA"] = np.random.default_rng(0).integers(500, 5000, adata.n_obs).astype(float)
        classify(adata)
        assert _schema(adata, "nCount_RNA")["role"] == Role.QC_METRIC.value

    def test_role_sample_id_via_known_name(self, adata):
        adata.obs["orig.ident"] = pd.Categorical(["donor_A", "donor_B"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "orig.ident")["role"] == Role.SAMPLE_ID.value

    def test_role_batch_via_known_name(self, adata):
        adata.obs["batch"] = pd.Categorical(["batch1", "batch2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "batch")["role"] == Role.BATCH.value

    def test_role_ambiguous_via_hint_name(self, adata):
        adata.obs["sample"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "sample")["role"] == Role.AMBIGUOUS.value

    def test_role_other_via_high_cardinality_string(self, adata):
        # 21+ unique string values, unknown name → OTHER
        vals = [f"cell_{i}" for i in range(adata.n_obs)]
        adata.obs["cell_barcode"] = vals
        classify(adata)
        assert _schema(adata, "cell_barcode")["role"] == Role.OTHER.value


# ---------------------------------------------------------------------------
# 2. Known-name table hits
# ---------------------------------------------------------------------------

class TestKnownNameTable:
    """Verify every key entry in the _KNOWN_NAMES contract."""

    # --- SAMPLE_ID ---
    def test_orig_ident(self, adata):
        adata.obs["orig.ident"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "orig.ident")["role"] == Role.SAMPLE_ID.value

    def test_donor(self, adata):
        adata.obs["donor"] = pd.Categorical(["d1", "d2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "donor")["role"] == Role.SAMPLE_ID.value

    def test_donor_id(self, adata):
        adata.obs["donor_id"] = pd.Categorical(["d1", "d2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "donor_id")["role"] == Role.SAMPLE_ID.value

    def test_subject(self, adata):
        adata.obs["subject"] = pd.Categorical(["sub1", "sub2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "subject")["role"] == Role.SAMPLE_ID.value

    def test_subject_id(self, adata):
        adata.obs["subject_id"] = pd.Categorical(["sub1", "sub2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "subject_id")["role"] == Role.SAMPLE_ID.value

    def test_sample_id_exact(self, adata):
        adata.obs["sample_id"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "sample_id")["role"] == Role.SAMPLE_ID.value

    # --- CLUSTERING ---
    def test_seurat_clusters(self, adata):
        adata.obs["seurat_clusters"] = pd.Categorical([str(i % 8) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "seurat_clusters")["role"] == Role.CLUSTERING.value

    def test_leiden(self, adata):
        adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "leiden")["role"] == Role.CLUSTERING.value

    def test_louvain(self, adata):
        adata.obs["louvain"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "louvain")["role"] == Role.CLUSTERING.value

    def test_cluster(self, adata):
        adata.obs["cluster"] = pd.Categorical([str(i % 6) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "cluster")["role"] == Role.CLUSTERING.value

    def test_clusters(self, adata):
        adata.obs["clusters"] = pd.Categorical([str(i % 6) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "clusters")["role"] == Role.CLUSTERING.value

    # --- CONDITION ---
    def test_phase(self, adata):
        adata.obs["phase"] = pd.Categorical((["G1", "S", "G2M"] * (adata.n_obs // 3 + 1))[:adata.n_obs])
        classify(adata)
        assert _schema(adata, "phase")["role"] == Role.CONDITION.value

    def test_genotype(self, adata):
        adata.obs["genotype"] = pd.Categorical(["WT", "KO"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "genotype")["role"] == Role.CONDITION.value

    def test_treatment_exact(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "treatment")["role"] == Role.CONDITION.value

    def test_condition_exact(self, adata):
        adata.obs["condition"] = pd.Categorical(["healthy", "disease"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "condition")["role"] == Role.CONDITION.value

    def test_timepoint(self, adata):
        adata.obs["timepoint"] = pd.Categorical((["d0", "d7", "d14"] * (adata.n_obs // 3 + 1))[:adata.n_obs])
        classify(adata)
        assert _schema(adata, "timepoint")["role"] == Role.CONDITION.value

    def test_sex(self, adata):
        adata.obs["sex"] = pd.Categorical(["M", "F"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "sex")["role"] == Role.CONDITION.value

    # --- CELL_TYPE ---
    def test_celltype(self, adata):
        adata.obs["celltype"] = pd.Categorical(["T cell", "B cell"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "celltype")["role"] == Role.CELL_TYPE.value

    def test_cell_types(self, adata):
        adata.obs["cell_types"] = pd.Categorical(["T cell", "B cell"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "cell_types")["role"] == Role.CELL_TYPE.value

    # --- BATCH ---
    def test_batch_exact(self, adata):
        adata.obs["batch"] = pd.Categorical(["b1", "b2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "batch")["role"] == Role.BATCH.value


# ---------------------------------------------------------------------------
# 3. Case-insensitive matching
# ---------------------------------------------------------------------------

class TestCaseInsensitive:
    """Known-name lookup is case-insensitive per the contract."""

    def test_orig_ident_title_case(self, adata):
        adata.obs["Orig.Ident"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "Orig.Ident")["role"] == Role.SAMPLE_ID.value

    def test_seurat_clusters_upper(self, adata):
        adata.obs["SEURAT_CLUSTERS"] = pd.Categorical([str(i % 8) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "SEURAT_CLUSTERS")["role"] == Role.CLUSTERING.value

    def test_leiden_upper(self, adata):
        adata.obs["LEIDEN"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "LEIDEN")["role"] == Role.CLUSTERING.value

    def test_phase_upper(self, adata):
        adata.obs["Phase"] = pd.Categorical((["G1", "S", "G2M"] * (adata.n_obs // 3 + 1))[:adata.n_obs])
        classify(adata)
        assert _schema(adata, "Phase")["role"] == Role.CONDITION.value

    def test_cell_type_mixed_case(self, adata):
        adata.obs["Cell_Type"] = pd.Categorical(["T cell", "B cell"] * (adata.n_obs // 2))
        classify(adata)
        assert _schema(adata, "Cell_Type")["role"] == Role.CELL_TYPE.value


# ---------------------------------------------------------------------------
# 4. _snn_res. substring pattern (Seurat family)
# ---------------------------------------------------------------------------

class TestSnnResPattern:
    """Any column whose lower-cased name contains '_snn_res.' → CLUSTERING."""

    def test_rna_snn_res(self, adata):
        adata.obs["RNA_snn_res.0.5"] = pd.Categorical([str(i % 10) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "RNA_snn_res.0.5")["role"] == Role.CLUSTERING.value

    def test_integrated_snn_res(self, adata):
        adata.obs["integrated_snn_res.0.8"] = pd.Categorical([str(i % 12) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "integrated_snn_res.0.8")["role"] == Role.CLUSTERING.value

    def test_sct_snn_res(self, adata):
        adata.obs["SCT_snn_res.1"] = pd.Categorical([str(i % 15) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "SCT_snn_res.1")["role"] == Role.CLUSTERING.value

    def test_snn_res_high_confidence(self, adata):
        adata.obs["RNA_snn_res.0.5"] = pd.Categorical([str(i % 10) for i in range(adata.n_obs)])
        classify(adata)
        assert _schema(adata, "RNA_snn_res.0.5")["confidence"] == "high"


# ---------------------------------------------------------------------------
# 5. QC pattern substring matching
# ---------------------------------------------------------------------------

class TestQcPatterns:
    """Columns matching QC substrings → QC_METRIC at high confidence."""

    def test_ncount_prefix(self, adata):
        adata.obs["nCount_RNA"] = np.arange(adata.n_obs, dtype=float)
        classify(adata)
        entry = _schema(adata, "nCount_RNA")
        assert entry["role"] == Role.QC_METRIC.value
        assert entry["confidence"] == "high"

    def test_nfeature_prefix(self, adata):
        adata.obs["nFeature_RNA"] = np.arange(adata.n_obs, dtype=float)
        classify(adata)
        assert _schema(adata, "nFeature_RNA")["role"] == Role.QC_METRIC.value

    def test_percent_mt(self, adata):
        adata.obs["percent.mt"] = np.random.default_rng(1).uniform(0, 30, adata.n_obs)
        classify(adata)
        assert _schema(adata, "percent.mt")["role"] == Role.QC_METRIC.value

    def test_pct_counts_prefix(self, adata):
        adata.obs["pct_counts_mt"] = np.random.default_rng(2).uniform(0, 30, adata.n_obs)
        classify(adata)
        assert _schema(adata, "pct_counts_mt")["role"] == Role.QC_METRIC.value

    def test_total_counts(self, adata):
        adata.obs["total_counts"] = np.arange(adata.n_obs, dtype=float)
        classify(adata)
        assert _schema(adata, "total_counts")["role"] == Role.QC_METRIC.value

    def test_n_genes(self, adata):
        adata.obs["n_genes"] = np.arange(adata.n_obs, dtype=float)
        classify(adata)
        assert _schema(adata, "n_genes")["role"] == Role.QC_METRIC.value

    def test_n_counts(self, adata):
        adata.obs["n_counts"] = np.arange(adata.n_obs, dtype=float)
        classify(adata)
        assert _schema(adata, "n_counts")["role"] == Role.QC_METRIC.value


# ---------------------------------------------------------------------------
# 6. Ambiguous hint names
# ---------------------------------------------------------------------------

class TestAmbiguousHints:
    """Names in {sample, group, ident} → AMBIGUOUS with candidate list."""

    def test_sample_hint(self, adata):
        adata.obs["sample"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        entry = _schema(adata, "sample")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert entry["confidence"] == "low"
        assert Role.CONDITION.value in entry["candidates"]
        assert Role.SAMPLE_ID.value in entry["candidates"]
        assert Role.BATCH.value in entry["candidates"]

    def test_group_hint(self, adata):
        adata.obs["group"] = pd.Categorical(["ctrl", "treat"] * (adata.n_obs // 2))
        classify(adata)
        entry = _schema(adata, "group")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert entry["candidates"] is not None and len(entry["candidates"]) > 0

    def test_ident_hint(self, adata):
        adata.obs["ident"] = pd.Categorical((["A", "B", "C"] * (adata.n_obs // 3 + 1))[:adata.n_obs])
        classify(adata)
        entry = _schema(adata, "ident")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert entry["candidates"] is not None and len(entry["candidates"]) > 0

    def test_ambiguous_hint_appears_in_ambiguous_cols(self, adata):
        adata.obs["sample"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        assert "sample" in _nvwa(adata)["ambiguous_cols"]

    def test_ambiguous_hint_not_in_condition_cols(self, adata):
        adata.obs["group"] = pd.Categorical(["ctrl", "treat"] * (adata.n_obs // 2))
        classify(adata)
        assert "group" not in _nvwa(adata)["condition_cols"]


# ---------------------------------------------------------------------------
# 7. Cardinality threshold = 20
# ---------------------------------------------------------------------------

class TestCardinalityThreshold:
    """19 unique values → low-cardinality path; 21 → high-cardinality path."""

    def test_numeric_low_cardinality_19_is_ambiguous(self, adata):
        # 19 distinct integer values, unknown name → numeric low-cardinality → AMBIGUOUS
        vals = [i % 19 for i in range(adata.n_obs)]
        adata.obs["my_score"] = np.array(vals, dtype=float)
        classify(adata)
        entry = _schema(adata, "my_score")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert entry["n_unique"] == 19

    def test_numeric_high_cardinality_21_is_qc(self, adata):
        # 21+ distinct float values, unknown name → numeric high-cardinality → QC_METRIC
        rng = np.random.default_rng(99)
        vals = rng.uniform(0, 1, adata.n_obs)
        # Guarantee 21+ unique values (continuous floats almost certainly unique)
        adata.obs["my_score2"] = vals
        classify(adata)
        entry = _schema(adata, "my_score2")
        assert entry["role"] == Role.QC_METRIC.value

    def test_string_low_cardinality_19_is_ambiguous(self, adata):
        # 19 unique string labels, unknown name → low-cardinality string → AMBIGUOUS
        vals = [f"cat_{i % 19}" for i in range(adata.n_obs)]
        adata.obs["my_label"] = pd.Categorical(vals)
        classify(adata)
        entry = _schema(adata, "my_label")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert entry["n_unique"] == 19

    def test_string_high_cardinality_21_is_other(self, adata):
        # 21+ unique string labels, unknown name → fallthrough → OTHER
        vals = [f"cell_{i}" for i in range(adata.n_obs)]  # each unique
        adata.obs["barcodes"] = vals
        classify(adata)
        entry = _schema(adata, "barcodes")
        assert entry["role"] == Role.OTHER.value

    def test_boundary_exactly_20_unique_is_low_cardinality(self, adata):
        # Exactly 20 unique values → <= 20 → low-cardinality AMBIGUOUS (unknown name)
        vals = [f"cat_{i % 20}" for i in range(adata.n_obs)]
        adata.obs["border_col"] = pd.Categorical(vals)
        classify(adata)
        entry = _schema(adata, "border_col")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert entry["n_unique"] == 20

    def test_boundary_exactly_21_unique_is_high_cardinality(self, adata):
        # 21 unique values → > 20 → high-cardinality → OTHER (unknown string name)
        vals = [f"cat_{i % 21}" for i in range(adata.n_obs)]
        adata.obs["border_col2"] = pd.Categorical(vals)
        classify(adata)
        entry = _schema(adata, "border_col2")
        assert entry["role"] == Role.OTHER.value
        assert entry["n_unique"] == 21

    def test_known_name_low_cardinality_populates_values(self, adata):
        # Known name + low cardinality → values field populated
        adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        entry = _schema(adata, "leiden")
        assert entry["values"] is not None

    def test_known_name_high_cardinality_no_values(self, adata):
        # Known name + > 20 unique values → values field is None
        adata.obs["cell_type"] = [f"cell_{i}" for i in range(adata.n_obs)]
        classify(adata)
        entry = _schema(adata, "cell_type")
        assert entry["values"] is None


# ---------------------------------------------------------------------------
# 8. Preview cap at 20 values
# ---------------------------------------------------------------------------

class TestPreviewCap:
    """_values_preview returns at most 20 entries."""

    def test_preview_capped_at_20(self, adata):
        # 19 unique categories → low-cardinality → values populated, capped at 20
        vals = [f"cat_{i % 19}" for i in range(adata.n_obs)]
        adata.obs["sample"] = pd.Categorical(vals)  # "sample" is an ambiguous hint
        classify(adata)
        entry = _schema(adata, "sample")
        assert entry["values"] is not None
        assert len(entry["values"]) <= _PREVIEW_CAP

    def test_preview_capped_when_exactly_20_unique(self, adata):
        # Known low-card name with 20 unique → values has exactly 20
        vals = [f"cat_{i % 20}" for i in range(adata.n_obs)]
        adata.obs["leiden"] = pd.Categorical(vals)
        classify(adata)
        entry = _schema(adata, "leiden")
        assert entry["values"] is not None
        assert len(entry["values"]) <= _PREVIEW_CAP

    def test_preview_values_are_strings(self, adata):
        # All entries in values must be plain strings
        vals = [i % 5 for i in range(adata.n_obs)]
        adata.obs["sample"] = pd.Categorical(vals)
        classify(adata)
        entry = _schema(adata, "sample")
        assert all(isinstance(v, str) for v in entry["values"])


# ---------------------------------------------------------------------------
# 9. Merge semantics — preserves other nvwa_meta keys
# ---------------------------------------------------------------------------

class TestMergeSemantics:
    """classify() overwrites schema/condition_cols/ambiguous_cols only."""

    def test_preserves_species_key(self, adata):
        adata.uns["nvwa_meta"] = {"species": {"value": "human", "confidence": "high"}}
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        meta = _nvwa(adata)
        assert "species" in meta
        assert meta["species"]["value"] == "human"

    def test_preserves_user_resolutions(self, adata):
        adata.uns["nvwa_meta"] = {"user_resolutions": {"sample": "condition"}}
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        meta = _nvwa(adata)
        assert meta["user_resolutions"]["sample"] == "condition"

    def test_overwrites_schema(self, adata):
        adata.uns["nvwa_meta"] = {"schema": {"stale_col": {"role": "other"}}}
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert "stale_col" not in _nvwa(adata)["schema"]
        assert "treatment" in _nvwa(adata)["schema"]

    def test_overwrites_condition_cols(self, adata):
        adata.uns["nvwa_meta"] = {"condition_cols": ["old_col"]}
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert "old_col" not in _nvwa(adata)["condition_cols"]
        assert "treatment" in _nvwa(adata)["condition_cols"]

    def test_overwrites_ambiguous_cols(self, adata):
        adata.uns["nvwa_meta"] = {"ambiguous_cols": ["old_ambiguous"]}
        adata.obs["sample"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        assert "old_ambiguous" not in _nvwa(adata)["ambiguous_cols"]
        assert "sample" in _nvwa(adata)["ambiguous_cols"]

    def test_preserves_arbitrary_extra_key(self, adata):
        adata.uns["nvwa_meta"] = {"custom_tag": "my_project"}
        adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        assert _nvwa(adata)["custom_tag"] == "my_project"

    def test_initializes_user_resolutions_if_absent(self, adata):
        adata.uns["nvwa_meta"] = {}
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert "user_resolutions" in _nvwa(adata)
        assert isinstance(_nvwa(adata)["user_resolutions"], dict)

    def test_does_not_overwrite_existing_user_resolutions(self, adata):
        adata.uns["nvwa_meta"] = {"user_resolutions": {"sample": "batch"}}
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert _nvwa(adata)["user_resolutions"]["sample"] == "batch"

    def test_idempotent_second_call(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        result1 = dict(_nvwa(adata))
        classify(adata)
        result2 = dict(_nvwa(adata))
        assert result1["schema"] == result2["schema"]
        assert result1["condition_cols"] == result2["condition_cols"]


# ---------------------------------------------------------------------------
# 10. SchemaEntry serialization — plain-str role, not Enum instance
# ---------------------------------------------------------------------------

class TestSchemaEntrySerialization:
    """Roles stored in uns must be plain str values (not Role enum instances)."""

    def test_role_is_plain_string_not_enum(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        role_val = _schema(adata, "treatment")["role"]
        # Must be a plain str, not a Role enum instance.
        # Role is a (str, Enum) subclass, so the distinguishing check is type() identity.
        assert type(role_val) is str, (
            f"Expected plain str, got {type(role_val)!r} — role was not unwrapped to .value"
        )

    def test_all_roles_in_schema_are_strings(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        adata.obs["sample"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        classify(adata)
        for col, entry in _nvwa(adata)["schema"].items():
            assert isinstance(entry["role"], str), f"Col {col!r} has non-str role: {type(entry['role'])}"

    def test_role_value_matches_enum_value(self, adata):
        adata.obs["batch"] = pd.Categorical(["b1", "b2"] * (adata.n_obs // 2))
        classify(adata)
        role_val = _schema(adata, "batch")["role"]
        assert role_val == Role.BATCH.value
        assert role_val == "batch"

    def test_schema_entry_has_required_fields(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        entry = _schema(adata, "treatment")
        for field in ("role", "confidence", "n_unique"):
            assert field in entry, f"Missing required field: {field}"

    def test_schema_entry_confidence_is_string(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        entry = _schema(adata, "treatment")
        assert entry["confidence"] in ("high", "low")

    def test_n_unique_is_int(self, adata):
        adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
        classify(adata)
        entry = _schema(adata, "leiden")
        assert isinstance(entry["n_unique"], int)


# ---------------------------------------------------------------------------
# 11. Derived output lists (condition_cols, ambiguous_cols)
# ---------------------------------------------------------------------------

class TestDerivedLists:
    """condition_cols contains high-confidence CONDITION cols; ambiguous_cols mirrors AMBIGUOUS."""

    def test_condition_cols_populated(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert "treatment" in _nvwa(adata)["condition_cols"]

    def test_condition_cols_high_confidence_only(self, adata):
        # Ambiguous-hint cols must NOT appear in condition_cols even if candidate is CONDITION
        adata.obs["group"] = pd.Categorical(["ctrl", "treat"] * (adata.n_obs // 2))
        classify(adata)
        assert "group" not in _nvwa(adata)["condition_cols"]

    def test_ambiguous_cols_populated(self, adata):
        adata.obs["ident"] = pd.Categorical(["A", "B"] * (adata.n_obs // 2))
        classify(adata)
        assert "ident" in _nvwa(adata)["ambiguous_cols"]

    def test_multiple_condition_cols(self, adata):
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        adata.obs["sex"] = pd.Categorical(["M", "F"] * (adata.n_obs // 2))
        classify(adata)
        cond = _nvwa(adata)["condition_cols"]
        assert "treatment" in cond and "sex" in cond

    def test_multiple_ambiguous_cols(self, adata):
        adata.obs["sample"] = pd.Categorical(["s1", "s2"] * (adata.n_obs // 2))
        adata.obs["group"] = pd.Categorical(["ctrl", "treat"] * (adata.n_obs // 2))
        classify(adata)
        ambig = _nvwa(adata)["ambiguous_cols"]
        assert "sample" in ambig and "group" in ambig


# ---------------------------------------------------------------------------
# 12. Edge cases — empty obs, all-NaN column
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge cases: no columns in obs, all-NaN series, no nvwa_meta pre-existing."""

    def test_empty_obs_columns(self, adata):
        # Strip all obs columns
        for col in list(adata.obs.columns):
            del adata.obs[col]
        classify(adata)
        meta = _nvwa(adata)
        assert meta["schema"] == {}
        assert meta["condition_cols"] == []
        assert meta["ambiguous_cols"] == []

    def test_all_nan_column_produces_sensible_role(self, adata):
        # A column of all NaN: dtype=float64, nunique(dropna=True)=0.
        # Falls into numeric branch (step 5), 0 <= 20 → AMBIGUOUS with n_unique=0.
        # Contract does not specify NaN-only behavior; this test pins current behavior
        # as a sentinel — if the output ever changes this will catch it.
        adata.obs["mystery_col"] = np.nan
        classify(adata)
        entry = _schema(adata, "mystery_col")
        assert entry["n_unique"] == 0
        assert entry["role"] == Role.AMBIGUOUS.value

    def test_all_nan_column_does_not_crash(self, adata):
        # classify() must not raise even for all-NaN columns
        adata.obs["mystery_col"] = np.nan
        classify(adata)  # should not raise
        assert "mystery_col" in _nvwa(adata)["schema"]

    def test_no_pre_existing_nvwa_meta(self, adata):
        # adata.uns starts without nvwa_meta — classify() must create it
        adata.uns.pop("nvwa_meta", None)
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        assert "nvwa_meta" in adata.uns
        assert "schema" in _nvwa(adata)

    def test_nvwa_meta_not_dict_is_reset(self, adata):
        # If nvwa_meta is not a dict (corrupted), classify() must handle gracefully
        adata.uns["nvwa_meta"] = "corrupt_string"
        adata.obs["treatment"] = pd.Categorical(["ctrl", "treated"] * (adata.n_obs // 2))
        classify(adata)
        meta = _nvwa(adata)
        assert isinstance(meta, dict)
        assert "schema" in meta


# ---------------------------------------------------------------------------
# 13. Numeric dtype path — explicit branching
# ---------------------------------------------------------------------------

class TestNumericDtype:
    """Numeric dtype falls through known-name / QC patterns into cardinality check."""

    def test_high_card_numeric_unknown_name_is_qc(self, adata):
        # > 20 unique floats, unknown name → QC_METRIC
        adata.obs["expression_score"] = np.random.default_rng(7).uniform(0, 1, adata.n_obs)
        classify(adata)
        assert _schema(adata, "expression_score")["role"] == Role.QC_METRIC.value

    def test_low_card_numeric_unknown_name_is_ambiguous(self, adata):
        # Exactly 3 unique values, unknown name, numeric → AMBIGUOUS
        adata.obs["cluster_num"] = np.array([0.0, 1.0, 2.0] * (adata.n_obs // 3 + 1))[:adata.n_obs]
        classify(adata)
        entry = _schema(adata, "cluster_num")
        assert entry["role"] == Role.AMBIGUOUS.value
        assert Role.CONDITION.value in entry["candidates"]
        assert Role.CLUSTERING.value in entry["candidates"]

    def test_low_card_numeric_ambiguous_has_values(self, adata):
        adata.obs["cluster_num"] = np.array([0.0, 1.0, 2.0] * (adata.n_obs // 3 + 1))[:adata.n_obs]
        classify(adata)
        entry = _schema(adata, "cluster_num")
        assert entry["values"] is not None

    def test_high_card_numeric_no_values_field(self, adata):
        adata.obs["expression_score"] = np.random.default_rng(8).uniform(0, 1, adata.n_obs)
        classify(adata)
        entry = _schema(adata, "expression_score")
        # High-cardinality QC_METRIC from numeric path: values not set
        assert entry.get("values") is None
