"""Layer 1 unit tests — Temple PI dataset (GSE223414_slim.h5ad).

Tests every Layer 1 tool against the real Temple dataset.
No mocking. No synthetic data. All assertions grounded in known dataset facts.

Dataset facts used in assertions:
  n_obs          = 45,460 cells
  n_vars         = 25,639 genes
  cell_type      = 12 types  (e.g. "Early cardiomyocyte" = 10,234 cells)
  orig.ident     = 8 conditions (Control-D5/10/14/30, PA-IVS-1v-D5/10/14/30)
  seurat_clusters= 20 clusters (strings "0"–"19")
  QC column      = "percent.mito"  (NOT pct_counts_mt)
  Embeddings     = X_pca, X_umap both present
  Key genes      = TNNT2, NKX2-5, MKI67, CD3D  (all present)

Run with:
  PYTHONPATH=. uv run pytest tests/unit/test_layer1_tools.py -v
"""

from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Ensure project root is importable regardless of cwd
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

DATASET_PATH = Path(
    os.environ.get("TEMPLE_DATASET_PATH", "/Users/zhuchen/Downloads/GSE223414_slim.h5ad")
)

PNG_MAGIC = b"\x89PNG"
MIN_PLOT_BYTES = 5_000  # same threshold as integration runner

# Plot images are saved here for CEO/CTO visual review
_PLOT_DIR = (
    _REPO_ROOT
    / "test-records"
    / f"{datetime.now().strftime('%Y-%m-%d')}_layer1-plots"
)


# ---------------------------------------------------------------------------
# Shared fixture — loaded once for the entire module
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def adata():
    """Load the Temple PI dataset once and share across all tests.

    The Temple dataset stores embeddings as 'PCA' and 'UMAP' (not 'X_pca'/'X_umap'),
    and both are DataFrames rather than numpy arrays. We convert and copy them to
    the scanpy-standard keys so all plotting functions work correctly.
    """
    import anndata as ad
    import numpy as np

    adata = ad.read_h5ad(DATASET_PATH)

    # Convert DataFrame embeddings to numpy arrays under scanpy-standard keys
    if "UMAP" in adata.obsm and "X_umap" not in adata.obsm:
        adata.obsm["X_umap"] = np.array(adata.obsm["UMAP"])
    if "PCA" in adata.obsm and "X_pca" not in adata.obsm:
        adata.obsm["X_pca"] = np.array(adata.obsm["PCA"])

    return adata


@pytest.fixture(scope="module")
def adata_with_de(adata):
    """AnnData with DE results pre-computed (one-vs-rest on cell_type).

    Stored at module scope so DE is only run once and reused by all
    tests that need rank_genes_groups in adata.uns.
    """
    import anndata as ad
    from src.analysis.differential import run_differential_expression

    # Work on a copy so the base fixture stays clean
    adata_copy = adata.copy()
    run_differential_expression(adata_copy, groupby="cell_type", n_genes=20)
    return adata_copy


def _valid_png(image: bytes) -> bool:
    """Return True if image is a real PNG above the minimum size threshold."""
    return (
        isinstance(image, (bytes, bytearray))
        and image[:4] == PNG_MAGIC
        and len(image) >= MIN_PLOT_BYTES
    )


@pytest.fixture(scope="module")
def plot_dir() -> Path:
    """Create and return the dated plot output directory in test-records/.

    All plot images generated during tests are saved here for CEO/CTO review.
    """
    _PLOT_DIR.mkdir(parents=True, exist_ok=True)
    return _PLOT_DIR


# ===========================================================================
# Section 1 — QC Metrics  (src/analysis/qc_metrics.py)
# ===========================================================================

class TestQCMetrics:

    def test_resolve_percent_mito(self, adata):
        """percent.mito is the mito column name in this dataset."""
        from src.analysis.qc_metrics import resolve_qc_metric_column
        result = resolve_qc_metric_column(adata, "percent.mito")
        assert result == "percent.mito"

    def test_resolve_ncount_rna(self, adata):
        """nCount_RNA resolves correctly."""
        from src.analysis.qc_metrics import resolve_qc_metric_column
        result = resolve_qc_metric_column(adata, "nCount_RNA")
        assert result == "nCount_RNA"

    def test_resolve_nfeature_rna(self, adata):
        """nFeature_RNA resolves correctly."""
        from src.analysis.qc_metrics import resolve_qc_metric_column
        result = resolve_qc_metric_column(adata, "nFeature_RNA")
        assert result == "nFeature_RNA"

    def test_resolve_unknown_metric_returns_none(self, adata):
        """Unknown metric name returns None, not an error."""
        from src.analysis.qc_metrics import resolve_qc_metric_column
        result = resolve_qc_metric_column(adata, "nonexistent_metric_xyz")
        assert result is None

    def test_get_obs_column_statistics_mito(self, adata):
        """percent.mito stats match known dataset values."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        stats = get_obs_column_statistics(adata, "percent.mito")

        assert stats["n_cells"] == 45_460
        assert abs(stats["mean"] - 7.31) < 0.5
        assert stats["min"] == pytest.approx(0.0, abs=0.01)
        assert stats["max"] == pytest.approx(24.99, abs=0.5)
        assert "median" in stats
        assert "q25" in stats
        assert "q75" in stats

    def test_get_obs_column_statistics_ncount(self, adata):
        """nCount_RNA stats match known dataset values."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        stats = get_obs_column_statistics(adata, "nCount_RNA")

        assert stats["min"] == pytest.approx(499, abs=1)
        assert stats["max"] == pytest.approx(19_991, abs=1)
        assert stats["n_cells"] == 45_460

    def test_get_obs_column_statistics_invalid_column(self, adata):
        """Non-existent column raises ValueError."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        with pytest.raises(ValueError):
            get_obs_column_statistics(adata, "nonexistent_col")

    def test_summarize_qc_metrics_returns_dataframe(self, adata):
        """summarize_qc_metrics returns a non-empty DataFrame."""
        import pandas as pd
        from src.analysis.qc_metrics import summarize_qc_metrics
        df = summarize_qc_metrics(adata)

        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "metric" in df.columns
        assert "mean" in df.columns
        assert "median" in df.columns

    def test_summarize_qc_metrics_covers_all_three(self, adata):
        """All three QC metrics are present when passed explicitly.

        Note: percent.mito uses a dot-separated name not in the auto-detect
        candidates list, so we pass all three explicitly.
        """
        from src.analysis.qc_metrics import summarize_qc_metrics
        df = summarize_qc_metrics(adata, metric_names=["nCount_RNA", "nFeature_RNA", "percent.mito"])
        metrics = df["metric"].tolist()

        assert "nCount_RNA" in metrics
        assert "nFeature_RNA" in metrics
        assert "percent.mito" in metrics


# ===========================================================================
# Section 2 — Composition Analysis  (src/analysis/composition.py)
# ===========================================================================

class TestComposition:

    def test_cross_tabulate_shape(self, adata):
        """cell_type × orig.ident produces (12, 8) table."""
        from src.analysis.composition import cross_tabulate_metadata
        result = cross_tabulate_metadata(adata, row_key="cell_type", col_key="orig.ident")

        assert result.shape == (12, 8)

    def test_cross_tabulate_total_cells(self, adata):
        """Sum of all cells in crosstab equals 45,460."""
        from src.analysis.composition import cross_tabulate_metadata
        result = cross_tabulate_metadata(adata, row_key="cell_type", col_key="orig.ident")

        assert result.sum().sum() == 45_460

    def test_cross_tabulate_known_cell_type_count(self, adata):
        """Early cardiomyocyte row sums to 10,234."""
        from src.analysis.composition import cross_tabulate_metadata
        result = cross_tabulate_metadata(adata, row_key="cell_type", col_key="orig.ident")

        assert result.loc["Early cardiomyocyte"].sum() == 10_234

    def test_cross_tabulate_no_negative_values(self, adata):
        """All cell counts are non-negative."""
        from src.analysis.composition import cross_tabulate_metadata
        result = cross_tabulate_metadata(adata, row_key="cell_type", col_key="orig.ident")

        assert (result >= 0).all().all()

    def test_cross_tabulate_invalid_row_key_raises(self, adata):
        """Invalid row_key raises ValueError."""
        from src.analysis.composition import cross_tabulate_metadata
        with pytest.raises(ValueError):
            cross_tabulate_metadata(adata, row_key="nonexistent_col", col_key="orig.ident")

    def test_cross_tabulate_invalid_col_key_raises(self, adata):
        """Invalid col_key raises ValueError."""
        from src.analysis.composition import cross_tabulate_metadata
        with pytest.raises(ValueError):
            cross_tabulate_metadata(adata, row_key="cell_type", col_key="nonexistent_col")

    # --- CEO bug: state contamination (A4) ---
    def test_ceo_state_contamination(self, adata):
        """Running the same query twice returns identical totals.

        CEO test A4: same query returned 45,460 vs 43,915 on two runs.
        Root cause was state contamination between test runs.
        """
        from src.analysis.composition import cross_tabulate_metadata

        result1 = cross_tabulate_metadata(adata, row_key="cell_type", col_key="orig.ident")
        result2 = cross_tabulate_metadata(adata, row_key="cell_type", col_key="orig.ident")

        assert result1.sum().sum() == result2.sum().sum() == 45_460
        assert result1.equals(result2)


# ===========================================================================
# Section 3 — Differential Expression  (src/analysis/differential.py)
# ===========================================================================

class TestDifferentialExpression:

    def test_run_de_all_groups_returns_deresult(self, adata):
        """run_differential_expression returns a DEResult."""
        from src.analysis.differential import DEResult, run_differential_expression
        result = run_differential_expression(adata.copy(), groupby="cell_type", n_genes=5)

        assert isinstance(result, DEResult)
        assert result.group == "all"
        assert result.reference == "rest"

    def test_run_de_result_has_expected_columns(self, adata):
        """Results DataFrame has gene, log2fc, pval, pval_adj columns."""
        import pandas as pd
        from src.analysis.differential import run_differential_expression
        result = run_differential_expression(adata.copy(), groupby="cell_type", n_genes=5)

        assert isinstance(result.results_df, pd.DataFrame)
        for col in ("gene", "log2fc", "pval", "pval_adj"):
            assert col in result.results_df.columns

    def test_run_de_target_group(self, adata):
        """target_group='Early cardiomyocyte' returns result for that group only."""
        from src.analysis.differential import run_differential_expression
        result = run_differential_expression(
            adata.copy(),
            groupby="cell_type",
            n_genes=10,
            target_group="Early cardiomyocyte",
        )

        assert result.group == "Early cardiomyocyte"
        assert result.reference == "rest"
        assert len(result.results_df) == 10

    def test_run_de_invalid_groupby_raises(self, adata):
        """Invalid groupby raises ValueError."""
        from src.analysis.differential import run_differential_expression
        with pytest.raises(ValueError):
            run_differential_expression(adata.copy(), groupby="nonexistent_col")

    def test_run_de_invalid_target_group_raises(self, adata):
        """Invalid target_group raises ValueError."""
        from src.analysis.differential import run_differential_expression
        with pytest.raises(ValueError):
            run_differential_expression(
                adata.copy(), groupby="cell_type", target_group="Fake Cell Type"
            )

    def test_run_pairwise_de_two_conditions(self, adata):
        """Pairwise DE between Control-D5 and PA-IVS-1v-D5 returns valid result."""
        from src.analysis.differential import DEResult, run_pairwise_de
        result = run_pairwise_de(
            adata.copy(),
            group1="Control-D5",
            group2="PA-IVS-1v-D5",
            groupby="orig.ident",
        )

        assert isinstance(result, DEResult)
        assert result.group == "Control-D5"
        assert result.reference == "PA-IVS-1v-D5"
        assert len(result.results_df) > 0

    def test_run_pairwise_de_invalid_group_raises(self, adata):
        """Pairwise DE with non-existent group raises ValueError."""
        from src.analysis.differential import run_pairwise_de
        with pytest.raises(ValueError):
            run_pairwise_de(
                adata.copy(),
                group1="Fake Condition",
                group2="Control-D5",
                groupby="orig.ident",
            )

    def test_get_de_dataframe_after_run(self, adata_with_de):
        """get_de_dataframe returns non-empty DataFrame after DE is run."""
        import pandas as pd
        from src.analysis.differential import get_de_dataframe
        df = get_de_dataframe(adata_with_de, group="Early cardiomyocyte")

        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "gene" in df.columns

    def test_get_de_dataframe_no_results_raises(self, adata):
        """get_de_dataframe raises ValueError when no DE has been run."""
        from src.analysis.differential import get_de_dataframe
        with pytest.raises(ValueError):
            get_de_dataframe(adata.copy(), group="Early cardiomyocyte")

    def test_run_de_n_genes_zero_returns_all(self, adata):
        """n_genes=0 sentinel returns all genes (not just top 20)."""
        from src.analysis.differential import run_differential_expression
        result = run_differential_expression(
            adata.copy(),
            groupby="cell_type",
            n_genes=0,
            target_group="Early cardiomyocyte",
        )

        assert len(result.results_df) == adata.n_vars

    # --- CEO bug: pairwise DE cell type filter (D2) ---
    def test_ceo_pairwise_de_cell_type_filter(self, adata):
        """DE within Early cardiomyocyte between two conditions completes.

        CEO test D2: de_within_celltype for 'Early cardiomyocyte'
        PA-IVS vs Control looped 4 times with no result.
        Root cause: cell type filter was broken.

        Note: Early cardiomyocyte subset only contains D10/D14 timepoints,
        not D5. We use the conditions that actually exist in this subset.
        """
        from src.analysis.differential import run_pairwise_de

        adata_sub = adata[adata.obs["cell_type"] == "Early cardiomyocyte"].copy()
        # Verify which conditions exist in this subset
        available = adata_sub.obs["orig.ident"].unique().tolist()
        control = next(c for c in available if c.startswith("Control"))
        disease = next(c for c in available if c.startswith("PA-IVS"))

        result = run_pairwise_de(
            adata_sub,
            group1=control,
            group2=disease,
            groupby="orig.ident",
        )

        assert result.group == control
        assert result.reference == disease
        assert len(result.results_df) > 0


# ===========================================================================
# Section 4 — Marker Genes  (src/analysis/marker_genes.py)
# ===========================================================================

class TestMarkerGenes:

    def test_get_top_markers_per_cluster_returns_list(self, adata_with_de):
        """get_top_marker_genes_per_cluster returns a list of strings."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        genes = get_top_marker_genes_per_cluster(adata_with_de, n_genes=3, groupby="cell_type")

        assert isinstance(genes, list)
        assert len(genes) > 0
        assert all(isinstance(g, str) for g in genes)

    def test_get_top_markers_per_cluster_count(self, adata_with_de):
        """Top 3 per cluster across 12 cell types yields at most 36 unique genes."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        genes = get_top_marker_genes_per_cluster(adata_with_de, n_genes=3, groupby="cell_type")

        # At most 12 cell types × 3 genes = 36 (may be fewer due to deduplication)
        assert 1 <= len(genes) <= 36

    def test_get_top_markers_per_cluster_exact_structure(self, adata_with_de):
        """get_top_marker_genes_per_cluster_exact returns dict with one entry per cell type."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster_exact
        result = get_top_marker_genes_per_cluster_exact(adata_with_de, n_genes=3, groupby="cell_type")

        assert isinstance(result, dict)
        assert len(result) == 12  # 12 cell types

    def test_get_top_markers_overall_returns_list(self, adata_with_de):
        """get_top_marker_genes_overall returns a list of the requested length."""
        from src.analysis.marker_genes import get_top_marker_genes_overall
        genes = get_top_marker_genes_overall(adata_with_de, n_genes=5)

        assert isinstance(genes, list)
        assert len(genes) == 5

    def test_no_de_results_raises(self, adata):
        """Marker gene functions raise ValueError when no DE has been run."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        with pytest.raises(ValueError):
            get_top_marker_genes_per_cluster(adata.copy(), n_genes=3)

    # --- CEO bug: dot plot top-N per cell type, not global (C6) ---
    def test_ceo_top_n_per_cell_type_not_global(self, adata_with_de):
        """Each cell type has its own top-3 genes, not the same 3 global genes.

        CEO test C6: asked for top 3 per cell type; returned 3 global genes
        (TGFBI, S100A11, IGFBP7) instead of per-cell-type lists.
        """
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster_exact
        result = get_top_marker_genes_per_cluster_exact(adata_with_de, n_genes=3, groupby="cell_type")

        # Every cell type must have its own gene list
        all_gene_lists = list(result.values())
        assert len(all_gene_lists) == 12

        # Not all cell types should have the exact same top gene
        first_genes = [genes[0] for genes in all_gene_lists if genes]
        assert len(set(first_genes)) > 1, (
            "All cell types returned the same top gene — this is the global-gene bug"
        )


# ===========================================================================
# Section 5 — Calculations  (src/analysis/calculations.py)
# ===========================================================================

class TestCalculations:

    def test_calculate_cluster_averages_shape(self, adata):
        """calculate_cluster_averages returns one row per cell type."""
        import pandas as pd
        from src.analysis.calculations import calculate_cluster_averages
        df = calculate_cluster_averages(adata, gene="TNNT2", groupby="cell_type")

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 12
        assert "cluster" in df.columns
        assert "mean_expression" in df.columns
        assert "cell_count" in df.columns

    def test_calculate_cluster_averages_invalid_gene(self, adata):
        """Non-existent gene raises ValueError."""
        from src.analysis.calculations import calculate_cluster_averages
        with pytest.raises(ValueError):
            calculate_cluster_averages(adata, gene="FAKEGENE999", groupby="cell_type")

    def test_calculate_cluster_averages_invalid_groupby(self, adata):
        """Non-existent groupby raises ValueError."""
        from src.analysis.calculations import calculate_cluster_averages
        with pytest.raises(ValueError):
            calculate_cluster_averages(adata, gene="TNNT2", groupby="nonexistent_col")

    def test_find_top_expressing_cluster_returns_tuple(self, adata):
        """find_top_expressing_cluster returns (str, float)."""
        from src.analysis.calculations import find_top_expressing_cluster
        cluster, mean_expr = find_top_expressing_cluster(adata, gene="TNNT2", groupby="cell_type")

        assert isinstance(cluster, str)
        assert isinstance(mean_expr, float)
        assert mean_expr >= 0.0

    def test_find_top_expressing_cluster_tnnt2_is_cardiomyocyte(self, adata):
        """TNNT2 highest expression is in a cardiomyocyte cell type."""
        from src.analysis.calculations import find_top_expressing_cluster
        cluster, _ = find_top_expressing_cluster(adata, gene="TNNT2", groupby="cell_type")

        assert "cardiomyocyte" in cluster.lower(), (
            f"Expected a cardiomyocyte cell type for TNNT2, got: '{cluster}'"
        )

    def test_get_metadata_summary_returns_dataframe(self, adata):
        """get_metadata_summary returns a DataFrame with categorical columns."""
        import pandas as pd
        from src.analysis.calculations import get_metadata_summary
        df = get_metadata_summary(adata)

        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "column_name" in df.columns

    def test_get_metadata_summary_includes_cell_type(self, adata):
        """cell_type and orig.ident appear in the metadata summary."""
        from src.analysis.calculations import get_metadata_summary
        df = get_metadata_summary(adata)
        cols = df["column_name"].tolist()

        assert "cell_type" in cols
        assert "orig.ident" in cols

    # --- CEO bug: MKI67 dimension — cell type not cluster ID (C1) ---
    def test_ceo_mki67_returns_cell_type_not_cluster_id(self, adata):
        """find_top_expressing_cluster with groupby=cell_type returns a cell type name.

        CEO test C1: asked which cell type expresses MKI67 highest;
        returned cluster IDs ("0", "Cluster 3") instead of cell type names.
        """
        from src.analysis.calculations import find_top_expressing_cluster
        cluster, _ = find_top_expressing_cluster(adata, gene="MKI67", groupby="cell_type")

        # Must be a real cell type name, not a numeric cluster ID
        assert not cluster.isdigit(), (
            f"Got cluster ID '{cluster}' instead of a cell type name — this is the MKI67 dimension bug"
        )
        # Must be one of the known cell types
        known_cell_types = {
            "Atrial-like cardiomyocyte", "Cardiac fibroblast", "Cardiac mesoderm",
            "Differentiating iPSC", "Early cardiomyocyte", "Epicardial progenitor",
            "FHF progenitor", "Left ventricular-like cardiomyocyte",
            "Right ventricular-like cardiomyocyte", "SHF progenitor",
            "Valvular endothelial cell", "Vascular endothelial cell",
        }
        assert cluster in known_cell_types, (
            f"'{cluster}' is not a known cell type in the Temple dataset"
        )


# ===========================================================================
# Section 6 — Plotting  (src/plotting/executor.py, comparison.py, volcano.py)
# ===========================================================================

class TestPlotting:

    def test_plot_umap_by_cell_type(self, adata, plot_dir):
        """plot_umap colored by cell_type returns valid PNG."""
        from src.plotting.executor import PlotResult, plot_umap
        result = plot_umap(adata, color="cell_type")

        assert isinstance(result, PlotResult)
        assert _valid_png(result.image)
        (plot_dir / "umap_by_cell_type.png").write_bytes(result.image)

    def test_plot_umap_by_seurat_clusters(self, adata, plot_dir):
        """plot_umap colored by seurat_clusters returns valid PNG."""
        from src.plotting.executor import plot_umap
        result = plot_umap(adata, color="seurat_clusters")

        assert _valid_png(result.image)
        (plot_dir / "umap_by_seurat_clusters.png").write_bytes(result.image)

    def test_plot_umap_invalid_color_raises(self, adata):
        """plot_umap with non-existent color key raises ValueError."""
        from src.plotting.executor import plot_umap
        with pytest.raises(ValueError):
            plot_umap(adata, color="nonexistent_col_xyz")

    def test_plot_violin_single_gene(self, adata, plot_dir):
        """plot_violin for TNNT2 grouped by cell_type returns valid PNG."""
        from src.plotting.executor import plot_violin
        result = plot_violin(adata, genes=["TNNT2"], groupby="cell_type")

        assert _valid_png(result.image)
        (plot_dir / "violin_TNNT2_by_cell_type.png").write_bytes(result.image)

    def test_plot_violin_qc_metric(self, adata, plot_dir):
        """plot_violin accepts QC obs columns (percent.mito) as input."""
        from src.plotting.executor import plot_violin
        result = plot_violin(adata, genes=["percent.mito"], groupby="cell_type")

        assert _valid_png(result.image)
        (plot_dir / "violin_pct_mito_by_cell_type.png").write_bytes(result.image)

    def test_plot_violin_invalid_gene_raises(self, adata):
        """plot_violin with non-existent gene raises ValueError."""
        from src.plotting.executor import plot_violin
        with pytest.raises(ValueError):
            plot_violin(adata, genes=["FAKEGENE999"], groupby="cell_type")

    def test_plot_feature_returns_valid_png(self, adata, plot_dir):
        """plot_feature for TNNT2 returns valid PNG."""
        from src.plotting.executor import plot_feature
        result = plot_feature(adata, gene="TNNT2")

        assert _valid_png(result.image)
        (plot_dir / "feature_TNNT2.png").write_bytes(result.image)

    def test_plot_feature_invalid_gene_raises(self, adata):
        """plot_feature with non-existent gene raises ValueError."""
        from src.plotting.executor import plot_feature
        with pytest.raises(ValueError):
            plot_feature(adata, gene="FAKEGENE999")

    def test_plot_dotplot_returns_valid_png(self, adata, plot_dir):
        """plot_dotplot for TNNT2 + NKX2-5 grouped by cell_type returns valid PNG."""
        from src.plotting.comparison import plot_dotplot
        result = plot_dotplot(adata, genes=["TNNT2", "NKX2-5"], groupby="cell_type")

        assert _valid_png(result.image)
        (plot_dir / "dotplot_TNNT2_NKX25_by_cell_type.png").write_bytes(result.image)

    def test_plot_heatmap_returns_valid_png(self, adata, plot_dir):
        """plot_heatmap for TNNT2 + NKX2-5 grouped by cell_type returns valid PNG."""
        from src.plotting.comparison import plot_heatmap
        result = plot_heatmap(adata, genes=["TNNT2", "NKX2-5"], groupby="cell_type")

        assert _valid_png(result.image)
        (plot_dir / "heatmap_TNNT2_NKX25_by_cell_type.png").write_bytes(result.image)

    def test_plot_scatter_qc_metrics(self, adata, plot_dir):
        """plot_scatter of nCount_RNA vs nFeature_RNA returns valid PNG."""
        from src.plotting.comparison import plot_scatter
        result = plot_scatter(adata, gene_x="nCount_RNA", gene_y="nFeature_RNA")

        assert _valid_png(result.image)
        (plot_dir / "scatter_nCount_vs_nFeature.png").write_bytes(result.image)

    def test_plot_volcano_returns_valid_png(self, adata_with_de, plot_dir):
        """plot_volcano from DE results returns valid PNG."""
        from src.analysis.differential import get_de_dataframe
        from src.plotting.volcano import plot_volcano

        de_df = get_de_dataframe(adata_with_de, group="Early cardiomyocyte")
        result = plot_volcano(de_df, group="Early cardiomyocyte")

        assert _valid_png(result.image)
        (plot_dir / "volcano_early_cardiomyocyte.png").write_bytes(result.image)

    # --- CEO bug: UMAP legend present (B1) ---
    def test_ceo_umap_legend_present(self, adata, plot_dir):
        """plot_umap with show_legend=True produces a larger image than without.

        CEO test B1: agent removed labels but claimed success.
        A plot with a legend is always larger in bytes than one without.
        """
        from src.plotting.executor import plot_umap

        with_legend = plot_umap(adata, color="cell_type", show_legend=True)
        without_legend = plot_umap(adata, color="cell_type", show_legend=False)

        assert _valid_png(with_legend.image)
        assert _valid_png(without_legend.image)
        assert len(with_legend.image) >= len(without_legend.image)
        (plot_dir / "ceo_b1_umap_WITH_legend.png").write_bytes(with_legend.image)
        (plot_dir / "ceo_b1_umap_WITHOUT_legend.png").write_bytes(without_legend.image)

    # --- CEO bug: split_by and color_by are independent (B2) ---
    def test_ceo_umap_split_by_and_color_by_independent(self, adata, plot_dir):
        """plot_umap with both split_by and color set independently does not crash.

        CEO test B2: split_by and color_by treated as the same parameter,
        required 3 correction rounds.
        """
        from src.plotting.executor import plot_umap

        result = plot_umap(adata, color="cell_type", split_by="orig.ident")

        assert _valid_png(result.image)
        assert "cell_type" in result.message or "split" in result.message
        (plot_dir / "ceo_b2_umap_split_by_condition_color_by_cell_type.png").write_bytes(result.image)

    def test_plot_feature_split_by(self, adata, plot_dir):
        """plot_feature with split_by produces multi-panel PNG."""
        from src.plotting.executor import plot_feature
        result = plot_feature(adata, gene="TNNT2", split_by="orig.ident")

        assert _valid_png(result.image)
        assert "split" in result.message.lower()
        (plot_dir / "feature_TNNT2_split_by_condition.png").write_bytes(result.image)
