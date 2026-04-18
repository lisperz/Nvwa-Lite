"""Layer 1 portable unit tests — works on any well-formed h5ad dataset.

Tests assert structural contracts and mathematical invariants of the analysis
layer. No biology-specific assertions (no hardcoded gene names, cell types,
or exact counts). Tests that require a feature the dataset lacks (e.g. UMAP
embeddings) are skipped rather than failed.

Run with:
  DATASET_PATH=/path/to/any.h5ad PYTHONPATH=. uv run pytest tests/unit/test_layer1_portable.py -v

Fixtures (adata, profile, adata_with_de, plot_dir) come from tests/conftest.py.
"""

from __future__ import annotations

import pytest

# conftest.py exports valid_png and _SAVE_PLOTS at module level
from tests.conftest import _SAVE_PLOTS, valid_png


# ===========================================================================
# Section 1 — QC Metrics  (src/analysis/qc_metrics.py)
# ===========================================================================

class TestQCMetrics:

    def test_resolve_unknown_metric_returns_none(self, adata):
        """Unknown metric name returns None, not an error."""
        from src.analysis.qc_metrics import resolve_qc_metric_column
        assert resolve_qc_metric_column(adata, "nonexistent_metric_xyz_abc") is None

    def test_resolve_known_metric_returns_string(self, adata, profile):
        """Any discovered QC column resolves to a non-empty string."""
        from src.analysis.qc_metrics import resolve_qc_metric_column
        col = profile.qc_counts_column or profile.qc_mito_column or profile.qc_genes_column
        if col is None:
            pytest.skip("No QC columns found in this dataset")
        result = resolve_qc_metric_column(adata, col)
        assert isinstance(result, str) and len(result) > 0

    def test_obs_statistics_n_cells_matches_dataset(self, adata, profile):
        """Statistics n_cells equals adata.n_obs."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        col = profile.qc_counts_column or profile.qc_mito_column or profile.qc_genes_column
        if col is None:
            pytest.skip("No numeric QC column found")
        stats = get_obs_column_statistics(adata, col)
        assert stats["n_cells"] == adata.n_obs

    def test_obs_statistics_mean_in_range(self, adata, profile):
        """Mean is between min and max."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        col = profile.qc_counts_column or profile.qc_mito_column or profile.qc_genes_column
        if col is None:
            pytest.skip("No numeric QC column found")
        stats = get_obs_column_statistics(adata, col)
        assert stats["min"] <= stats["mean"] <= stats["max"]

    def test_obs_statistics_quartiles_ordered(self, adata, profile):
        """q25 <= median <= q75."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        col = profile.qc_counts_column or profile.qc_mito_column or profile.qc_genes_column
        if col is None:
            pytest.skip("No numeric QC column found")
        stats = get_obs_column_statistics(adata, col)
        assert stats["q25"] <= stats["median"] <= stats["q75"]

    def test_obs_statistics_invalid_column_raises(self, adata):
        """Non-existent column raises ValueError."""
        from src.analysis.qc_metrics import get_obs_column_statistics
        with pytest.raises(ValueError):
            get_obs_column_statistics(adata, "nonexistent_col_xyz")

    def test_summarize_qc_returns_dataframe(self, adata, profile):
        """summarize_qc_metrics returns a non-empty DataFrame with expected columns."""
        import pandas as pd
        from src.analysis.qc_metrics import summarize_qc_metrics
        cols = [c for c in [profile.qc_counts_column, profile.qc_genes_column, profile.qc_mito_column] if c]
        if not cols:
            pytest.skip("No QC columns found in this dataset")
        df = summarize_qc_metrics(adata, metric_names=cols)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        for expected_col in ("metric", "mean", "median"):
            assert expected_col in df.columns


# ===========================================================================
# Section 2 — Composition Analysis  (src/analysis/composition.py)
# ===========================================================================

class TestComposition:

    def test_cross_tabulate_sum_equals_n_cells(self, adata, profile):
        """Cross-tabulation total equals dataset cell count."""
        from src.analysis.composition import cross_tabulate_metadata
        if not profile.cell_type_column or not profile.condition_column:
            pytest.skip("Need both cell type and condition columns")
        ct = cross_tabulate_metadata(adata, row_key=profile.cell_type_column, col_key=profile.condition_column)
        assert ct.values.sum() == adata.n_obs

    def test_cross_tabulate_shape_matches_uniques(self, adata, profile):
        """Shape matches unique value counts of each column."""
        from src.analysis.composition import cross_tabulate_metadata
        if not profile.cell_type_column or not profile.condition_column:
            pytest.skip("Need both cell type and condition columns")
        ct = cross_tabulate_metadata(adata, row_key=profile.cell_type_column, col_key=profile.condition_column)
        n_types = adata.obs[profile.cell_type_column].nunique()
        n_conds = adata.obs[profile.condition_column].nunique()
        assert ct.shape == (n_types, n_conds)

    def test_cross_tabulate_no_negative_values(self, adata, profile):
        """All cell counts are non-negative."""
        from src.analysis.composition import cross_tabulate_metadata
        if not profile.cell_type_column or not profile.condition_column:
            pytest.skip("Need both cell type and condition columns")
        ct = cross_tabulate_metadata(adata, row_key=profile.cell_type_column, col_key=profile.condition_column)
        assert (ct >= 0).all().all()

    def test_cross_tabulate_invalid_key_raises(self, adata, profile):
        """Invalid row_key raises ValueError."""
        from src.analysis.composition import cross_tabulate_metadata
        col = profile.cell_type_column or profile.cluster_column or adata.obs.columns[0]
        with pytest.raises(ValueError):
            cross_tabulate_metadata(adata, row_key="nonexistent_col_xyz", col_key=col)

    def test_cross_tabulate_same_query_idempotent(self, adata, profile):
        """Running the same query twice returns identical results (state isolation)."""
        from src.analysis.composition import cross_tabulate_metadata
        if not profile.cell_type_column or not profile.condition_column:
            pytest.skip("Need both cell type and condition columns")
        r1 = cross_tabulate_metadata(adata, row_key=profile.cell_type_column, col_key=profile.condition_column)
        r2 = cross_tabulate_metadata(adata, row_key=profile.cell_type_column, col_key=profile.condition_column)
        assert r1.equals(r2)
        assert r1.values.sum() == r2.values.sum()


# ===========================================================================
# Section 3 — Differential Expression  (src/analysis/differential.py)
# ===========================================================================

class TestDifferentialExpression:

    def test_run_de_returns_deresult(self, adata, profile):
        """run_differential_expression returns a DEResult."""
        from src.analysis.differential import DEResult, run_differential_expression
        result = run_differential_expression(adata.copy(), groupby=profile.primary_groupby, n_genes=5)
        assert isinstance(result, DEResult)

    def test_run_de_result_has_expected_columns(self, adata, profile):
        """Results DataFrame has gene, log2fc, pval, pval_adj columns."""
        import pandas as pd
        from src.analysis.differential import run_differential_expression
        result = run_differential_expression(adata.copy(), groupby=profile.primary_groupby, n_genes=5)
        assert isinstance(result.results_df, pd.DataFrame)
        for col in ("gene", "log2fc", "pval", "pval_adj"):
            assert col in result.results_df.columns

    def test_run_de_genes_are_in_var_names(self, adata, profile):
        """All returned gene names exist in adata.var_names or adata.raw.var_names."""
        from src.analysis.differential import run_differential_expression
        result = run_differential_expression(adata.copy(), groupby=profile.primary_groupby, n_genes=5)
        raw_var_names = set(adata.raw.var_names) if adata.raw is not None else set()
        all_var_names = set(adata.var_names) | raw_var_names
        for gene in result.results_df["gene"]:
            assert gene in all_var_names, f"Gene '{gene}' not in var_names or raw.var_names"

    def test_run_de_invalid_groupby_raises(self, adata):
        """Invalid groupby raises ValueError."""
        from src.analysis.differential import run_differential_expression
        with pytest.raises(ValueError):
            run_differential_expression(adata.copy(), groupby="nonexistent_col_xyz")

    def test_run_de_n_genes_zero_returns_all(self, adata, profile):
        """n_genes=0 sentinel returns all genes."""
        from src.analysis.differential import run_differential_expression
        groups = adata.obs[profile.primary_groupby].unique()
        first_group = str(groups[0])
        result = run_differential_expression(
            adata.copy(), groupby=profile.primary_groupby, n_genes=0, target_group=first_group
        )
        assert len(result.results_df) == adata.n_vars

    def test_get_de_dataframe_no_results_raises(self, adata, profile):
        """get_de_dataframe raises ValueError when no DE has been run."""
        from src.analysis.differential import get_de_dataframe
        groups = adata.obs[profile.primary_groupby].unique()
        with pytest.raises(ValueError):
            get_de_dataframe(adata.copy(), group=str(groups[0]))

    def test_get_de_dataframe_after_run(self, adata_with_de, profile):
        """get_de_dataframe returns non-empty DataFrame after DE is run."""
        import pandas as pd
        from src.analysis.differential import get_de_dataframe
        groups = adata_with_de.obs[profile.primary_groupby].unique()
        df = get_de_dataframe(adata_with_de, group=str(groups[0]))
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "gene" in df.columns

    def test_pairwise_de_returns_deresult(self, adata, profile):
        """Pairwise DE between two groups returns a valid DEResult."""
        from src.analysis.differential import DEResult, run_pairwise_de
        groups = adata.obs[profile.primary_groupby].unique()
        if len(groups) < 2:
            pytest.skip("Need at least 2 groups for pairwise DE")
        result = run_pairwise_de(
            adata.copy(),
            group1=str(groups[0]),
            group2=str(groups[1]),
            groupby=profile.primary_groupby,
        )
        assert isinstance(result, DEResult)
        assert len(result.results_df) > 0


# ===========================================================================
# Section 4 — Marker Genes  (src/analysis/marker_genes.py)
# ===========================================================================

class TestMarkerGenes:

    def test_get_top_markers_per_cluster_returns_list(self, adata_with_de, profile):
        """get_top_marker_genes_per_cluster returns a non-empty list of strings."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        genes = get_top_marker_genes_per_cluster(adata_with_de, n_genes=3, groupby=profile.primary_groupby)
        assert isinstance(genes, list)
        assert len(genes) > 0
        assert all(isinstance(g, str) for g in genes)

    def test_get_top_markers_are_valid_gene_names(self, adata_with_de, profile):
        """All returned marker genes exist in adata.var_names or adata.raw.var_names."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        genes = get_top_marker_genes_per_cluster(adata_with_de, n_genes=3, groupby=profile.primary_groupby)
        raw_var_names = set(adata_with_de.raw.var_names) if adata_with_de.raw is not None else set()
        all_var_names = set(adata_with_de.var_names) | raw_var_names
        for gene in genes:
            assert gene in all_var_names, f"Marker gene '{gene}' not in var_names or raw.var_names"

    def test_get_top_markers_count_bounded(self, adata_with_de, profile):
        """Top 3 per group yields at most n_groups * 3 unique genes."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        n_groups = adata_with_de.obs[profile.primary_groupby].nunique()
        genes = get_top_marker_genes_per_cluster(adata_with_de, n_genes=3, groupby=profile.primary_groupby)
        assert 1 <= len(genes) <= n_groups * 3

    def test_get_top_markers_exact_structure(self, adata_with_de, profile):
        """get_top_marker_genes_per_cluster_exact returns dict with one entry per group."""
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster_exact
        n_groups = adata_with_de.obs[profile.primary_groupby].nunique()
        result = get_top_marker_genes_per_cluster_exact(adata_with_de, n_genes=3, groupby=profile.primary_groupby)
        assert isinstance(result, dict)
        assert len(result) == n_groups

    def test_top_markers_not_all_same_first_gene(self, adata_with_de, profile):
        """Each group has its own top gene — not all groups return the same gene.

        Portable version of CEO bug C6 (top-N per cell type, not global).
        """
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster_exact
        n_groups = adata_with_de.obs[profile.primary_groupby].nunique()
        if n_groups < 2:
            pytest.skip("Need at least 2 groups to check per-group diversity")
        result = get_top_marker_genes_per_cluster_exact(adata_with_de, n_genes=3, groupby=profile.primary_groupby)
        first_genes = [genes[0] for genes in result.values() if genes]
        assert len(set(first_genes)) > 1, "All groups returned the same top gene — global-gene bug"

    def test_no_de_results_raises(self, adata, profile):
        """Marker gene functions raise ValueError when no DE has been run.

        Skip if the dataset already contains pre-baked rank_genes_groups (e.g. PBMC
        tutorial h5ad files ship with DE results in uns).
        """
        import anndata as ad
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster
        fresh = ad.AnnData(X=adata.X, obs=adata.obs.copy(), var=adata.var.copy())
        with pytest.raises(ValueError):
            get_top_marker_genes_per_cluster(fresh, n_genes=3)


# ===========================================================================
# Section 5 — Calculations  (src/analysis/calculations.py)
# ===========================================================================

class TestCalculations:

    def test_calculate_cluster_averages_columns(self, adata, profile):
        """calculate_cluster_averages returns DataFrame with required columns."""
        import pandas as pd
        from src.analysis.calculations import calculate_cluster_averages
        gene = profile.sample_genes[0] if profile.sample_genes else None
        if gene is None:
            pytest.skip("No sample genes available")
        df = calculate_cluster_averages(adata, gene=gene, groupby=profile.primary_groupby)
        assert isinstance(df, pd.DataFrame)
        for col in ("cluster", "mean_expression", "cell_count"):
            assert col in df.columns

    def test_calculate_cluster_averages_cell_count_sum(self, adata, profile):
        """Cell counts in cluster averages sum to total n_cells."""
        from src.analysis.calculations import calculate_cluster_averages
        gene = profile.sample_genes[0] if profile.sample_genes else None
        if gene is None:
            pytest.skip("No sample genes available")
        df = calculate_cluster_averages(adata, gene=gene, groupby=profile.primary_groupby)
        assert df["cell_count"].sum() == adata.n_obs

    def test_calculate_cluster_averages_invalid_gene_raises(self, adata, profile):
        """Non-existent gene raises ValueError."""
        from src.analysis.calculations import calculate_cluster_averages
        with pytest.raises(ValueError):
            calculate_cluster_averages(adata, gene="FAKEGENE_XYZ_999", groupby=profile.primary_groupby)

    def test_find_top_expressing_cluster_returns_valid_group(self, adata, profile):
        """find_top_expressing_cluster returns a group that exists in the dataset."""
        from src.analysis.calculations import find_top_expressing_cluster
        gene = profile.sample_genes[0] if profile.sample_genes else None
        if gene is None:
            pytest.skip("No sample genes available")
        cluster, mean_expr = find_top_expressing_cluster(adata, gene=gene, groupby=profile.primary_groupby)
        assert isinstance(cluster, str)
        assert isinstance(mean_expr, float)
        assert mean_expr >= 0.0
        valid_groups = adata.obs[profile.primary_groupby].astype(str).unique().tolist()
        assert cluster in valid_groups, f"Returned cluster '{cluster}' not in {profile.primary_groupby}"

    def test_get_cluster_statistics_percentages_sum(self, adata, profile):
        """Cluster percentages sum to approximately 100%."""
        from src.analysis.calculations import get_cluster_statistics
        df = get_cluster_statistics(adata, groupby=profile.primary_groupby)
        assert abs(df["percentage"].sum() - 100.0) < 0.01

    def test_get_metadata_summary_returns_dataframe(self, adata):
        """get_metadata_summary returns a DataFrame with column_name column."""
        import pandas as pd
        from src.analysis.calculations import get_metadata_summary
        df = get_metadata_summary(adata)
        assert isinstance(df, pd.DataFrame)
        assert "column_name" in df.columns

    def test_find_top_cluster_not_numeric_when_groupby_is_annotation(self, adata, profile):
        """When groupby is a cell type annotation column, result is not a bare numeric ID.

        Portable version of CEO bug C1 (MKI67 returned cluster ID instead of cell type).
        """
        from src.analysis.calculations import find_top_expressing_cluster
        if not profile.cell_type_column:
            pytest.skip("No cell type annotation column found")
        gene = profile.sample_genes[0] if profile.sample_genes else None
        if gene is None:
            pytest.skip("No sample genes available")
        cluster, _ = find_top_expressing_cluster(adata, gene=gene, groupby=profile.cell_type_column)
        assert not cluster.isdigit(), (
            f"Got numeric cluster ID '{cluster}' when groupby='{profile.cell_type_column}' — dimension bug"
        )


# ===========================================================================
# Section 6 — Plotting  (src/plotting/executor.py, comparison.py, volcano.py)
# ===========================================================================

class TestPlotting:

    def test_plot_umap_by_groupby(self, adata, profile, plot_dir):
        """plot_umap colored by primary_groupby returns valid PNG."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        from src.plotting.executor import PlotResult, plot_umap
        result = plot_umap(adata, color=profile.primary_groupby)
        assert isinstance(result, PlotResult)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / f"umap_by_{profile.primary_groupby}.png").write_bytes(result.image)

    def test_plot_umap_by_gene(self, adata, profile, plot_dir):
        """plot_umap colored by a gene returns valid PNG."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        if not profile.sample_genes:
            pytest.skip("No sample genes available")
        from src.plotting.executor import plot_umap
        result = plot_umap(adata, color=profile.sample_genes[0])
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / f"umap_gene_{profile.sample_genes[0]}.png").write_bytes(result.image)

    def test_plot_umap_split_by_condition(self, adata, profile, plot_dir):
        """plot_umap with split_by returns valid PNG."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        if not profile.condition_column:
            pytest.skip("Dataset has no condition column")
        from src.plotting.executor import plot_umap
        result = plot_umap(adata, color=profile.primary_groupby, split_by=profile.condition_column)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / f"umap_split_by_{profile.condition_column}.png").write_bytes(result.image)

    def test_plot_umap_invalid_color_raises(self, adata, profile):
        """plot_umap with non-existent color key raises ValueError."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        from src.plotting.executor import plot_umap
        with pytest.raises(ValueError):
            plot_umap(adata, color="nonexistent_col_xyz_abc")

    def test_plot_violin_single_gene(self, adata, profile, plot_dir):
        """plot_violin for a sample gene returns valid PNG."""
        if not profile.sample_genes:
            pytest.skip("No sample genes available")
        from src.plotting.executor import plot_violin
        result = plot_violin(adata, genes=[profile.sample_genes[0]], groupby=profile.primary_groupby)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / f"violin_{profile.sample_genes[0]}.png").write_bytes(result.image)

    def test_plot_violin_invalid_gene_raises(self, adata, profile):
        """plot_violin with non-existent gene raises ValueError."""
        from src.plotting.executor import plot_violin
        with pytest.raises(ValueError):
            plot_violin(adata, genes=["FAKEGENE_XYZ_999"], groupby=profile.primary_groupby)

    def test_plot_feature_returns_valid_png(self, adata, profile, plot_dir):
        """plot_feature for a sample gene returns valid PNG."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        if not profile.sample_genes:
            pytest.skip("No sample genes available")
        from src.plotting.executor import plot_feature
        result = plot_feature(adata, gene=profile.sample_genes[0])
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / f"feature_{profile.sample_genes[0]}.png").write_bytes(result.image)

    def test_plot_feature_invalid_gene_raises(self, adata, profile):
        """plot_feature with non-existent gene raises ValueError."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        from src.plotting.executor import plot_feature
        with pytest.raises(ValueError):
            plot_feature(adata, gene="FAKEGENE_XYZ_999")

    def test_plot_feature_split_by(self, adata, profile, plot_dir):
        """plot_feature with split_by produces multi-panel PNG."""
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        if not profile.condition_column:
            pytest.skip("Dataset has no condition column")
        if not profile.sample_genes:
            pytest.skip("No sample genes available")
        from src.plotting.executor import plot_feature
        result = plot_feature(adata, gene=profile.sample_genes[0], split_by=profile.condition_column)
        assert valid_png(result.image)
        assert "split" in result.message.lower()
        if _SAVE_PLOTS:
            (plot_dir / f"feature_{profile.sample_genes[0]}_split.png").write_bytes(result.image)

    def test_plot_dotplot_returns_valid_png(self, adata, profile, plot_dir):
        """plot_dotplot for sample genes returns valid PNG."""
        if len(profile.sample_genes) < 2:
            pytest.skip("Need at least 2 sample genes")
        from src.plotting.comparison import plot_dotplot
        result = plot_dotplot(adata, genes=profile.sample_genes[:2], groupby=profile.primary_groupby)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / "dotplot_sample_genes.png").write_bytes(result.image)

    def test_plot_umap_with_legend_larger_than_without(self, adata, profile, plot_dir):
        """plot_umap with show_legend=True produces larger image than without.

        Portable version of CEO bug B1 (legend removed but agent claimed success).
        """
        if not profile.has_umap:
            pytest.skip("Dataset has no X_umap embedding")
        from src.plotting.executor import plot_umap
        with_legend = plot_umap(adata, color=profile.primary_groupby, show_legend=True)
        without_legend = plot_umap(adata, color=profile.primary_groupby, show_legend=False)
        assert valid_png(with_legend.image)
        assert valid_png(without_legend.image)
        assert len(with_legend.image) >= len(without_legend.image)
        if _SAVE_PLOTS:
            (plot_dir / "umap_with_legend.png").write_bytes(with_legend.image)
            (plot_dir / "umap_without_legend.png").write_bytes(without_legend.image)

    def test_plot_volcano_returns_valid_png(self, adata_with_de, profile, plot_dir):
        """plot_volcano from DE results returns valid PNG."""
        from src.analysis.differential import get_de_dataframe
        from src.plotting.volcano import plot_volcano
        groups = adata_with_de.obs[profile.primary_groupby].unique()
        group = str(groups[0])
        de_df = get_de_dataframe(adata_with_de, group=group)
        result = plot_volcano(de_df, group=group)
        assert valid_png(result.image)
        if _SAVE_PLOTS:
            (plot_dir / f"volcano_{group[:20]}.png").write_bytes(result.image)


# ===========================================================================
# Section 7 — Pairwise DE + Volcano Integration (Issue #35)
# ===========================================================================

class TestPairwiseVolcanoIntegration:
    """Tests for issue #35: volcano plot after pairwise DE.

    Verifies:
    1. Pairwise DE does not overwrite rank_genes_groups
    2. Volcano works for regular DE groups
    3. Volcano works for pairwise "X vs Y" keys
    4. Invalid volcano group raises ValueError (logs as error, not success)
    """

    def test_pairwise_de_preserves_rank_genes_groups(self, adata, profile):
        """Pairwise DE does not overwrite prior all-vs-rest DE results."""
        from src.analysis.differential import run_differential_expression, run_pairwise_de

        groups = adata.obs[profile.primary_groupby].unique()
        if len(groups) < 2:
            pytest.skip("Need at least 2 groups for pairwise DE")

        adata_copy = adata.copy()

        # Step 1: Run all-vs-rest DE
        run_differential_expression(
            adata_copy, groupby=profile.primary_groupby, method="wilcoxon", n_genes=10
        )
        original_groups = list(adata_copy.uns["rank_genes_groups"]["names"].dtype.names)
        assert len(original_groups) >= 2

        # Step 2: Run pairwise DE
        run_pairwise_de(
            adata_copy,
            group1=str(groups[0]),
            group2=str(groups[1]),
            groupby=profile.primary_groupby,
        )

        # Step 3: Verify rank_genes_groups still has all original groups
        preserved_groups = list(adata_copy.uns["rank_genes_groups"]["names"].dtype.names)
        assert preserved_groups == original_groups, (
            f"Pairwise DE corrupted rank_genes_groups. "
            f"Expected {len(original_groups)} groups, got {len(preserved_groups)}"
        )

    def test_volcano_works_for_regular_de_group(self, adata, profile):
        """Volcano plot works for regular DE group after all-vs-rest."""
        from src.analysis.differential import get_de_dataframe, run_differential_expression
        from src.plotting.volcano import plot_volcano

        groups = adata.obs[profile.primary_groupby].unique()
        if len(groups) < 1:
            pytest.skip("Need at least 1 group for DE")

        adata_copy = adata.copy()
        run_differential_expression(
            adata_copy, groupby=profile.primary_groupby, method="wilcoxon", n_genes=10
        )

        group = str(groups[0])
        de_df = get_de_dataframe(adata_copy, group=group)
        result = plot_volcano(de_df, group=group)
        assert valid_png(result.image)

    def test_volcano_works_for_pairwise_key(self, adata, profile):
        """Volcano plot works for pairwise 'X vs Y' key after compare_groups_de."""
        from src.analysis.differential import get_de_dataframe, run_pairwise_de

        groups = adata.obs[profile.primary_groupby].unique()
        if len(groups) < 2:
            pytest.skip("Need at least 2 groups for pairwise DE")

        adata_copy = adata.copy()
        group1 = str(groups[0])
        group2 = str(groups[1])

        # Run pairwise DE and store result in pairwise_de_result
        result = run_pairwise_de(
            adata_copy, group1=group1, group2=group2, groupby=profile.primary_groupby
        )
        adata_copy.uns["pairwise_de_result"] = {
            "group1": group1,
            "group2": group2,
            "results_df": result.results_df,
        }

        # Volcano should work with "X vs Y" key
        pairwise_key = f"{group1} vs {group2}"
        de_df = get_de_dataframe(adata_copy, group=pairwise_key)
        assert len(de_df) > 0
        assert "gene" in de_df.columns
        assert "log2fc" in de_df.columns

    def test_invalid_volcano_group_raises_valueerror(self, adata, profile):
        """Invalid volcano group raises ValueError (not returns error string)."""
        from src.analysis.differential import get_de_dataframe, run_differential_expression

        adata_copy = adata.copy()
        run_differential_expression(
            adata_copy, groupby=profile.primary_groupby, method="wilcoxon", n_genes=10
        )

        # Invalid group should raise, not return error string
        with pytest.raises(ValueError, match="not in DE results"):
            get_de_dataframe(adata_copy, group="nonexistent_group_xyz_abc")

