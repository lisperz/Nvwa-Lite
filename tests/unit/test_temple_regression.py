"""Temple cardiac dataset regression tests — runs only when DATASET_NAME=temple.

These tests assert biology-specific facts and reproduce known CEO bug cases.
They are intentionally hardcoded to the Temple PI dataset (GSE223414_slim.h5ad)
because they document real production incidents that must never regress.

To run:
  DATASET_PATH=/path/to/GSE223414_slim.h5ad DATASET_NAME=temple \\
      PYTHONPATH=. uv run pytest tests/unit/test_temple_regression.py -v

Dataset facts used in assertions:
  n_obs          = 45,460 cells
  n_vars         = 25,639 genes
  cell_type      = 12 types  ("Early cardiomyocyte" = 10,234 cells)
  orig.ident     = 8 conditions (Control-D5/10/14/30, PA-IVS-1v-D5/10/14/30)
  seurat_clusters= 20 clusters (strings "0"–"19")
  Key genes      = TNNT2, MKI67, NKX2-5, CD3D
"""

from __future__ import annotations

import pytest

from tests.conftest import _SAVE_PLOTS, valid_png

# Skip the entire module if not running against the Temple dataset
pytestmark = pytest.mark.skipif(
    __import__("os").environ.get("DATASET_NAME", "").lower() != "temple",
    reason="Temple regression tests only run when DATASET_NAME=temple",
)

# ---------------------------------------------------------------------------
# Known Temple dataset constants
# ---------------------------------------------------------------------------
_TEMPLE_N_CELLS = 45_460
_TEMPLE_N_CELL_TYPES = 12
_TEMPLE_N_CONDITIONS = 8
_TEMPLE_EARLY_CARDIO_COUNT = 10_234
_TEMPLE_CELL_TYPE_COL = "cell_type"
_TEMPLE_CONDITION_COL = "orig.ident"
_TEMPLE_CLUSTER_COL = "seurat_clusters"
_TEMPLE_KNOWN_CELL_TYPES = {
    "Atrial-like cardiomyocyte",
    "Cardiac fibroblast",
    "Cardiac mesoderm",
    "Differentiating iPSC",
    "Early cardiomyocyte",
    "Epicardial progenitor",
    "FHF progenitor",
    "Left ventricular-like cardiomyocyte",
    "Right ventricular-like cardiomyocyte",
    "SHF progenitor",
    "Valvular endothelial cell",
    "Vascular endothelial cell",
}


# ===========================================================================
# Section 1 — Exact Dataset Counts
# ===========================================================================

class TestTempleExactCounts:

    def test_n_cells(self, adata):
        """Dataset has exactly 45,460 cells."""
        assert adata.n_obs == _TEMPLE_N_CELLS

    def test_n_cell_types(self, adata):
        """Dataset has exactly 12 cell types."""
        assert adata.obs[_TEMPLE_CELL_TYPE_COL].nunique() == _TEMPLE_N_CELL_TYPES

    def test_n_conditions(self, adata):
        """Dataset has exactly 8 conditions."""
        assert adata.obs[_TEMPLE_CONDITION_COL].nunique() == _TEMPLE_N_CONDITIONS

    def test_cross_tabulate_shape(self, adata):
        """cell_type × orig.ident produces (12, 8) table."""
        from src.analysis.composition import cross_tabulate_metadata
        ct = cross_tabulate_metadata(adata, row_key=_TEMPLE_CELL_TYPE_COL, col_key=_TEMPLE_CONDITION_COL)
        assert ct.shape == (_TEMPLE_N_CELL_TYPES, _TEMPLE_N_CONDITIONS)

    def test_cross_tabulate_total_cells(self, adata):
        """Sum of all cells in crosstab equals 45,460."""
        from src.analysis.composition import cross_tabulate_metadata
        ct = cross_tabulate_metadata(adata, row_key=_TEMPLE_CELL_TYPE_COL, col_key=_TEMPLE_CONDITION_COL)
        assert ct.sum().sum() == _TEMPLE_N_CELLS

    def test_early_cardiomyocyte_count(self, adata):
        """Early cardiomyocyte row sums to exactly 10,234."""
        from src.analysis.composition import cross_tabulate_metadata
        ct = cross_tabulate_metadata(adata, row_key=_TEMPLE_CELL_TYPE_COL, col_key=_TEMPLE_CONDITION_COL)
        assert ct.loc["Early cardiomyocyte"].sum() == _TEMPLE_EARLY_CARDIO_COUNT


# ===========================================================================
# Section 2 — Temple Biology Assertions
# ===========================================================================

class TestTempleBiology:

    def test_tnnt2_highest_in_cardiomyocyte(self, adata):
        """TNNT2 highest expression is in a cardiomyocyte cell type."""
        from src.analysis.calculations import find_top_expressing_cluster
        cluster, _ = find_top_expressing_cluster(adata, gene="TNNT2", groupby=_TEMPLE_CELL_TYPE_COL)
        assert "cardiomyocyte" in cluster.lower(), (
            f"Expected a cardiomyocyte cell type for TNNT2, got: '{cluster}'"
        )

    def test_known_genes_in_var_names(self, adata):
        """Key Temple genes exist in var_names."""
        for gene in ("TNNT2", "MKI67", "NKX2-5", "CD3D"):
            assert gene in adata.var_names, f"Expected gene '{gene}' not found in Temple dataset"

    def test_qc_mito_stats(self, adata):
        """percent.mito mean matches known Temple value (~7.31)."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        stats = get_obs_column_statistics(adata, "percent.mito")
        assert stats["n_cells"] == _TEMPLE_N_CELLS
        assert abs(stats["mean"] - 7.31) < 0.5
        assert stats["min"] == pytest.approx(0.0, abs=0.01)
        assert stats["max"] == pytest.approx(24.99, abs=0.5)

    def test_known_cell_types_present(self, adata):
        """All 12 known Temple cell types are present."""
        actual = set(adata.obs[_TEMPLE_CELL_TYPE_COL].unique())
        assert actual == _TEMPLE_KNOWN_CELL_TYPES


# ===========================================================================
# Section 3 — CEO Bug Regressions (preserved verbatim)
# ===========================================================================

class TestCEOBugRegressions:

    def test_ceo_state_contamination(self, adata):
        """Running the same query twice returns identical totals.

        CEO test A4: same query returned 45,460 vs 43,915 on two runs.
        Root cause was state contamination between test runs.
        """
        from src.analysis.composition import cross_tabulate_metadata

        result1 = cross_tabulate_metadata(adata, row_key=_TEMPLE_CELL_TYPE_COL, col_key=_TEMPLE_CONDITION_COL)
        result2 = cross_tabulate_metadata(adata, row_key=_TEMPLE_CELL_TYPE_COL, col_key=_TEMPLE_CONDITION_COL)

        assert result1.sum().sum() == result2.sum().sum() == _TEMPLE_N_CELLS
        assert result1.equals(result2)

    def test_ceo_pairwise_de_cell_type_filter(self, adata):
        """DE within Early cardiomyocyte between two conditions completes.

        CEO test D2: de_within_celltype for 'Early cardiomyocyte'
        PA-IVS vs Control looped 4 times with no result.
        Root cause: cell type filter was broken.

        Note: Early cardiomyocyte subset only contains D10/D14 timepoints,
        not D5. We use the conditions that actually exist in this subset.
        """
        from src.analysis.differential import run_pairwise_de

        adata_sub = adata[adata.obs[_TEMPLE_CELL_TYPE_COL] == "Early cardiomyocyte"].copy()
        available = adata_sub.obs[_TEMPLE_CONDITION_COL].unique().tolist()
        control = next(c for c in available if c.startswith("Control"))
        disease = next(c for c in available if c.startswith("PA-IVS"))

        result = run_pairwise_de(
            adata_sub,
            group1=control,
            group2=disease,
            groupby=_TEMPLE_CONDITION_COL,
        )

        assert result.group == control
        assert result.reference == disease
        assert len(result.results_df) > 0

    def test_ceo_top_n_per_cell_type_not_global(self, adata_with_de):
        """Each cell type has its own top-3 genes, not the same 3 global genes.

        CEO test C6: asked for top 3 per cell type; returned 3 global genes
        (TGFBI, S100A11, IGFBP7) instead of per-cell-type lists.
        """
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster_exact

        result = get_top_marker_genes_per_cluster_exact(
            adata_with_de, n_genes=3, groupby=_TEMPLE_CELL_TYPE_COL
        )

        all_gene_lists = list(result.values())
        assert len(all_gene_lists) == _TEMPLE_N_CELL_TYPES

        first_genes = [genes[0] for genes in all_gene_lists if genes]
        assert len(set(first_genes)) > 1, (
            "All cell types returned the same top gene — this is the global-gene bug"
        )

    def test_ceo_mki67_returns_cell_type_not_cluster_id(self, adata):
        """find_top_expressing_cluster with groupby=cell_type returns a cell type name.

        CEO test C1: asked which cell type expresses MKI67 highest;
        returned cluster IDs ("0", "Cluster 3") instead of cell type names.
        """
        from src.analysis.calculations import find_top_expressing_cluster

        cluster, _ = find_top_expressing_cluster(adata, gene="MKI67", groupby=_TEMPLE_CELL_TYPE_COL)

        assert not cluster.isdigit(), (
            f"Got cluster ID '{cluster}' instead of a cell type name — MKI67 dimension bug"
        )
        assert cluster in _TEMPLE_KNOWN_CELL_TYPES, (
            f"'{cluster}' is not a known cell type in the Temple dataset"
        )


# ===========================================================================
# Section 4 — Temple-Specific Plots
# ===========================================================================

class TestTemplePlots:

    def test_umap_by_cell_type(self, adata, plot_dir):
        """plot_umap colored by cell_type returns valid PNG."""
        from src.plotting.executor import plot_umap
        result = plot_umap(adata, color=_TEMPLE_CELL_TYPE_COL)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / "temple_umap_by_cell_type.png").write_bytes(result.image)

    def test_violin_tnnt2_by_cell_type(self, adata, plot_dir):
        """plot_violin for TNNT2 grouped by cell_type returns valid PNG."""
        from src.plotting.executor import plot_violin
        result = plot_violin(adata, genes=["TNNT2"], groupby=_TEMPLE_CELL_TYPE_COL)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / "temple_violin_TNNT2_by_cell_type.png").write_bytes(result.image)

    def test_feature_tnnt2(self, adata, plot_dir):
        """plot_feature for TNNT2 returns valid PNG."""
        from src.plotting.executor import plot_feature
        result = plot_feature(adata, gene="TNNT2")
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / "temple_feature_TNNT2.png").write_bytes(result.image)

    def test_feature_tnnt2_split_by_condition(self, adata, plot_dir):
        """plot_feature for TNNT2 split by condition produces multi-panel PNG."""
        from src.plotting.executor import plot_feature
        result = plot_feature(adata, gene="TNNT2", split_by=_TEMPLE_CONDITION_COL)
        assert valid_png(result.image)
        assert "split" in result.message.lower()
        if _SAVE_PLOTS:
            (plot_dir / "temple_feature_TNNT2_split_by_condition.png").write_bytes(result.image)

    def test_ceo_umap_legend_present(self, adata, plot_dir):
        """plot_umap with show_legend=True produces larger image than without.

        CEO test B1: agent removed labels but claimed success.
        """
        from src.plotting.executor import plot_umap
        with_legend = plot_umap(adata, color=_TEMPLE_CELL_TYPE_COL, show_legend=True)
        without_legend = plot_umap(adata, color=_TEMPLE_CELL_TYPE_COL, show_legend=False)
        assert valid_png(with_legend.image)
        assert valid_png(without_legend.image)
        assert len(with_legend.image) >= len(without_legend.image)
        if _SAVE_PLOTS:
            (plot_dir / "temple_ceo_b1_umap_WITH_legend.png").write_bytes(with_legend.image)
            (plot_dir / "temple_ceo_b1_umap_WITHOUT_legend.png").write_bytes(without_legend.image)

    def test_ceo_umap_split_by_and_color_by_independent(self, adata, plot_dir):
        """plot_umap with split_by and color set independently does not crash.

        CEO test B2: split_by and color_by treated as the same parameter.
        """
        from src.plotting.executor import plot_umap
        result = plot_umap(adata, color=_TEMPLE_CELL_TYPE_COL, split_by=_TEMPLE_CONDITION_COL)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / "temple_ceo_b2_umap_split_by_condition.png").write_bytes(result.image)
