"""LangChain tool definitions for the visualization agent.

Each tool wraps a plotting/analysis function, providing validation,
error handling, and structured output for the agent.
"""

from __future__ import annotations

import io
import logging
from typing import TYPE_CHECKING, Callable

import pandas as pd
from langchain_core.tools import tool

from src.agent import analysis_tools
from src.agent import subset_tools
from src.agent.viz_state import update_viz_state
from src.analysis.calculations import calculate_mito_percentage, get_metadata_summary
from src.analysis.gene_lookup import lookup_gene_name
from src.analysis.composition import cross_tabulate_metadata
from src.analysis.differential import get_de_dataframe, run_differential_expression, run_pairwise_de, get_all_de_results
from src.analysis.marker_genes import get_top_marker_genes_per_cluster
from src.analysis.preprocessing import run_preprocessing
from src.analysis.qc_metrics import get_obs_column_statistics, resolve_qc_metric_column, summarize_qc_metrics
from src.plotting.comparison import plot_dotplot, plot_heatmap, plot_scatter
from src.plotting.composition import plot_composition
from src.plotting.executor import PlotResult, TableResult, plot_feature, plot_umap, plot_violin
from src.plotting.volcano import plot_volcano
from src.types import DatasetState, detect_dataset_state

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

_adata: AnnData | None = None
_dataset_state: DatasetState | None = None
_plot_results: list[PlotResult] = []
_table_results: list[TableResult] = []
_adata_replaced_callback: Callable[[AnnData], None] | None = None
_plot_generated_callback: Callable[[], None] | None = None


def bind_dataset(adata: AnnData) -> None:
    """Bind a dataset so tools can access it."""
    global _adata  # noqa: PLW0603
    _adata = adata
    # Also bind to analysis_tools and subset_tools modules
    analysis_tools._adata = adata
    subset_tools._adata = adata


def bind_dataset_state(state: DatasetState) -> None:
    """Bind the dataset state for tools to reference."""
    global _dataset_state  # noqa: PLW0603
    _dataset_state = state


def set_adata_replaced_callback(cb: Callable[[AnnData], None]) -> None:
    """Register a callback invoked when preprocessing replaces adata."""
    global _adata_replaced_callback  # noqa: PLW0603
    _adata_replaced_callback = cb


def set_plot_generated_callback(cb: Callable[[], None]) -> None:
    """Register a callback invoked when any plot is generated."""
    global _plot_generated_callback  # noqa: PLW0603
    _plot_generated_callback = cb
    subset_tools._plot_generated_callback = cb


def _get_adata() -> AnnData:
    if _adata is None:
        raise RuntimeError("No dataset loaded.")
    return _adata


def _get_cluster_key() -> str:
    """Auto-detect the clustering key from the dataset."""
    if _dataset_state and _dataset_state.cluster_key:
        return _dataset_state.cluster_key

    # Fallback detection
    adata = _get_adata()
    for key in ("leiden", "louvain"):
        if key in adata.obs.columns:
            return key

    # If no clustering found, return leiden as default
    return "leiden"


def _update_state() -> None:
    """Re-detect dataset state after a mutation."""
    global _dataset_state  # noqa: PLW0603
    if _adata is not None and _dataset_state is not None:
        _dataset_state = detect_dataset_state(
            _adata, _dataset_state.source, _dataset_state.filename,
        )


def _store_and_return(result: PlotResult) -> str:
    """Append a PlotResult to the buffer and return a text summary."""
    _plot_results.append(result)
    logger.info(f"[FEEDBACK] _store_and_return called, callback exists: {_plot_generated_callback is not None}")
    if _plot_generated_callback:
        logger.info("[FEEDBACK] Invoking plot generated callback")
        _plot_generated_callback()
    else:
        logger.warning("[FEEDBACK] No plot generated callback registered!")
    return f"Plot generated successfully.\nCode: {result.code}\n{result.message}"


def _store_table_and_return(result: TableResult) -> str:
    """Append a TableResult to the buffer and return a text summary."""
    _table_results.append(result)
    return f"Table generated successfully.\nCode: {result.code}\n{result.message}"


def get_plot_results() -> list[PlotResult]:
    """Retrieve all plot results from this turn (called by UI)."""
    # Combine results from all modules
    all_results = list(_plot_results)
    all_results.extend(analysis_tools._plot_results)
    all_results.extend(subset_tools._plot_results)
    return all_results


def get_table_results() -> list[TableResult]:
    """Retrieve all table results from this turn (called by UI)."""
    return list(_table_results)


def clear_plot_results() -> None:
    """Clear the plot result buffer."""
    _plot_results.clear()
    analysis_tools._plot_results.clear()
    subset_tools._plot_results.clear()


def clear_table_results() -> None:
    """Clear the table result buffer."""
    _table_results.clear()


@tool
def umap_plot(color_by: str, show_labels: bool = False, show_legend: bool = True, split_by: str = "") -> str:
    """Generate a UMAP plot colored by an observation key or gene name.

    Args:
        color_by: The observation column (e.g. 'louvain') or gene name (e.g. 'CD3E').
        show_labels: Whether to show cluster labels directly on the plot. Defaults to False.
        show_legend: Whether to show the legend. Defaults to True.
        split_by: Optional observation key to split the plot into separate panels (e.g. 'louvain' to show each cluster separately).
    """
    adata = _get_adata()

    # Auto-detect clustering key for split_by if requested
    if split_by == "cluster" or split_by == "clusters":
        split_by = _get_cluster_key()

    try:
        result = _store_and_return(plot_umap(
            adata, color=color_by, show_labels=show_labels, show_legend=show_legend,
            split_by=split_by if split_by else None
        ))
        update_viz_state("umap", color_by=color_by, split_by=split_by or None,
                         show_labels=show_labels, show_legend=show_legend)
        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("UMAP plot failed")
        return f"Unexpected error: {e}"


@tool
def violin_plot(genes: str, groupby: str = "") -> str:
    """Generate a violin plot showing gene expression distribution across groups.

    Args:
        genes: Gene name(s) to plot. Can be a single gene (e.g. 'MS4A1') or
               comma-separated genes (e.g. 'CD3E,CD8A').
        groupby: The observation key to group cells by. If empty, auto-detects clustering key.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        # Handle both single gene and multiple genes
        if "," in genes:
            gene_list = [g.strip() for g in genes.split(",")]
        else:
            gene_list = genes.strip()
        result = _store_and_return(plot_violin(adata, genes=gene_list, groupby=groupby))
        genes_list = gene_list if isinstance(gene_list, list) else [gene_list]
        update_viz_state("violin", groupby=groupby, genes=genes_list)
        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Violin plot failed")
        return f"Unexpected error: {e}"


@tool
def dotplot(genes: str, groupby: str = "") -> str:
    """Generate a dot plot for one or more genes across cell groups.

    Args:
        genes: Comma-separated gene names (e.g. 'CD3E,MS4A1,NKG7').
        groupby: The observation key to group cells by. If empty, auto-detects clustering key.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    gene_list = [g.strip() for g in genes.split(",")]
    try:
        result = _store_and_return(plot_dotplot(adata, genes=gene_list, groupby=groupby))
        update_viz_state("dotplot", groupby=groupby, genes=gene_list)
        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Dot plot failed")
        return f"Unexpected error: {e}"


@tool
def dotplot_combined(genes: str, row_key: str, col_key: str) -> str:
    """Generate a dot plot showing gene expression across TWO dimensions simultaneously.

    Creates a combined grouping (row_key + col_key) and plots as a single vertical list.
    The result is a long list of combined labels (e.g., "Cardiomyocyte + Control-D10").

    For cross-cell-type comparisons, prefer dotplot_matrix which creates a cleaner
    hierarchical matrix layout. Use this as a fallback when dotplot_matrix is insufficient.

    Args:
        genes: Comma-separated gene names (e.g., "TNNT2,NKX2-5")
        row_key: First grouping variable (e.g., "cell_type")
        col_key: Second grouping variable (e.g., "orig.ident")

    Returns:
        Success message with plot details, or error message.
    """
    adata = _get_adata()

    # Validate columns exist
    if row_key not in adata.obs.columns:
        return f"Error: '{row_key}' not found in observations. Available: {', '.join(adata.obs.columns)}"
    if col_key not in adata.obs.columns:
        return f"Error: '{col_key}' not found in observations. Available: {', '.join(adata.obs.columns)}"

    gene_list = [g.strip() for g in genes.split(",")]
    combined_key = f"{row_key}_x_{col_key}"

    try:
        # Create combined grouping column
        adata.obs[combined_key] = (
            adata.obs[row_key].astype(str) + " + " + adata.obs[col_key].astype(str)
        )

        # Generate dotplot with combined grouping
        result = plot_dotplot(adata, genes=gene_list, groupby=combined_key)

        # Clean up temporary column
        del adata.obs[combined_key]

        # Enhance result with context
        enhanced_result = PlotResult(
            image=result.image,
            code=f'# Combined grouping: {row_key} × {col_key}\n' + result.code,
            message=f"Dot plot of {', '.join(gene_list)} across {row_key} × {col_key} (combined grouping)."
        )

        result_str = _store_and_return(enhanced_result)
        update_viz_state("dotplot_combined", row_key=row_key, col_key=col_key, genes=gene_list)
        return result_str

    except ValueError as e:
        # Clean up on error
        if combined_key in adata.obs.columns:
            del adata.obs[combined_key]
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Combined dot plot failed")
        if combined_key in adata.obs.columns:
            del adata.obs[combined_key]
        return f"Unexpected error: {e}"


@tool
def dotplot_matrix(genes: str, cell_type_key: str, condition_key: str) -> str:
    """Generate a matrix-style dot plot for cell_type × condition comparison.

    This creates a hierarchical dot plot where:
    - X-axis shows cell types grouped by condition (hierarchical layout)
    - Y-axis shows genes
    - Dot size = fraction of cells expressing the gene
    - Dot color = mean expression level

    Use this for questions like:
    - "How does [gene] vary across cell types in disease vs normal?"
    - "Compare [gene] expression across cell types and conditions"
    - "Show [gene] in different cell types split by condition"

    Args:
        genes: Comma-separated gene names (e.g., "NKX2-5" or "NKX2-5,TNNT2")
        cell_type_key: Column name for cell type grouping (e.g., "cell_type")
        condition_key: Column name for condition grouping (e.g., "orig.ident")

    Returns:
        Success message with plot details, or error message.

    Example:
        dotplot_matrix(genes="NKX2-5", cell_type_key="cell_type", condition_key="orig.ident")
    """
    import io as _io
    import matplotlib.pyplot as plt

    adata = _get_adata()
    gene_list = [g.strip() for g in genes.split(",")]

    # Validate keys
    if cell_type_key not in adata.obs.columns:
        return f"Error: '{cell_type_key}' not found in adata.obs"
    if condition_key not in adata.obs.columns:
        return f"Error: '{condition_key}' not found in adata.obs"

    # Validate genes (check var_names and raw)
    available_genes = set(adata.var_names)
    if adata.raw is not None:
        available_genes.update(adata.raw.var_names)
    missing = [g for g in gene_list if g not in available_genes]
    if missing:
        return f"Error: Genes not found: {missing}"

    try:
        import scanpy as sc

        fig = sc.pl.dotplot(
            adata,
            var_names=gene_list,
            groupby=[cell_type_key, condition_key],
            return_fig=True,
        )
        fig.savefig(buf := _io.BytesIO(), format="png", bbox_inches="tight")
        plt.close("all")
        buf.seek(0)
        image_bytes = buf.read()

        n_cell_types = adata.obs[cell_type_key].nunique()
        n_conditions = adata.obs[condition_key].nunique()

        code = (
            f'sc.pl.dotplot(adata, var_names={gene_list},\n'
            f'    groupby=["{cell_type_key}", "{condition_key}"])'
        )
        message = (
            f"Matrix dot plot of {', '.join(gene_list)}: "
            f"{n_cell_types} cell types × {n_conditions} conditions. "
            f"X-axis: cell types grouped by condition (hierarchical). "
            f"Y-axis: genes. Dot size = fraction expressing, dot color = mean expression."
        )

        result = PlotResult(image=image_bytes, code=code, message=message)
        result_str = _store_and_return(result)
        update_viz_state("dotplot_matrix", cell_type_key=cell_type_key,
                         condition_key=condition_key, genes=gene_list)
        return result_str

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Matrix dot plot failed")
        return f"Unexpected error: {e}"


@tool
def feature_plot(gene: str, split_by: str = "") -> str:
    """Generate a feature plot showing gene expression on UMAP (viridis colormap).

    Args:
        gene: The gene name to visualize.
        split_by: Optional observation key to split the plot into separate panels
                  (e.g. 'orig.ident' to show one panel per condition/sample).
                  Each panel uses the same color scale for direct comparison.
    """
    adata = _get_adata()
    try:
        result = _store_and_return(plot_feature(
            adata, gene=gene, split_by=split_by if split_by else None,
        ))
        update_viz_state("feature", genes=[gene], split_by=split_by or None)
        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Feature plot failed")
        return f"Unexpected error: {e}"


@tool
def heatmap_plot(genes: str, groupby: str = "", n_genes_per_cluster: int = 0) -> str:
    """Generate a heatmap for genes across cell groups.

    Args:
        genes: Comma-separated gene names (e.g. 'CD3E,MS4A1,NKG7,LYZ').
        groupby: The observation key to group cells by. If empty, auto-detects clustering key.
        n_genes_per_cluster: Optional context about how many genes per cluster were requested (for display purposes).
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    gene_list = [g.strip() for g in genes.split(",")]
    try:
        result = _store_and_return(plot_heatmap(
            adata, genes=gene_list, groupby=groupby, n_genes_per_cluster=n_genes_per_cluster
        ))
        update_viz_state("heatmap", groupby=groupby, genes=gene_list)
        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Heatmap plot failed")
        return f"Unexpected error: {e}"


@tool
def scatter_plot(gene_x: str, gene_y: str, color_by: str = "") -> str:
    """Generate a scatter plot showing correlation between two genes.

    Args:
        gene_x: Gene name for x-axis (e.g. 'CD3E').
        gene_y: Gene name for y-axis (e.g. 'CD8A').
        color_by: Optional observation key or gene to color points by. Leave empty for no coloring.
    """
    adata = _get_adata()
    try:
        color = color_by if color_by else None
        result = _store_and_return(plot_scatter(adata, gene_x=gene_x, gene_y=gene_y, color_by=color))
        update_viz_state("scatter", color_by=color_by or None, genes=[gene_x, gene_y])
        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Scatter plot failed")
        return f"Unexpected error: {e}"


@tool
def volcano_plot_tool(group: str) -> str:
    """Generate a volcano plot for a specific cluster from DE results.

    Args:
        group: The cluster/group name to plot (e.g. '0', '1').
    """
    adata = _get_adata()
    try:
        de_df = get_de_dataframe(adata, group=group)
        return _store_and_return(plot_volcano(de_df, group=group))
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Volcano plot failed")
        return f"Unexpected error: {e}"


@tool
def dataset_info() -> str:
    """Get dataset metadata: genes, observation keys, cell count, and processing state."""
    adata = _get_adata()
    obs_keys = ", ".join(sorted(adata.obs.columns))

    gene_set: set[str] = set(adata.var_names)
    if adata.raw is not None:
        gene_set.update(adata.raw.var_names)
    all_genes = sorted(gene_set)

    state_info = ""
    if _dataset_state is not None:
        state_info = f"\n- Processing state: {_dataset_state.summary()}"

    # Cap gene list to avoid blowing up LLM context on large datasets
    max_genes_shown = 50
    if len(all_genes) > max_genes_shown:
        gene_preview = ", ".join(all_genes[:max_genes_shown])
        gene_line = f"- Sample genes (first {max_genes_shown} of {len(all_genes)}): {gene_preview}"
    else:
        gene_line = f"- All genes: {', '.join(all_genes)}"

    return (
        f"Dataset summary:\n"
        f"- Cells: {adata.n_obs}\n"
        f"- Total genes: {len(all_genes)} "
        f"(processed: {adata.n_vars}, raw: {adata.raw.n_vars if adata.raw else 0})\n"
        f"- Observation keys: {obs_keys}\n"
        f"{gene_line}"
        f"{state_info}"
    )


@tool
def check_data_status() -> str:
    """Check what preprocessing steps have been applied to the current dataset."""
    if _dataset_state is None:
        return "No dataset state available."
    return _dataset_state.summary()


@tool
def inspect_metadata(max_unique_values: int = 50) -> str:
    """Inspect categorical metadata columns in adata.obs and return their statistics.

    This tool automatically scans adata.obs for categorical columns (like louvain,
    cell_type, batch, sample, etc.) and returns their unique values and counts.

    Use this to:
    - Discover what metadata is available in the dataset
    - See cluster/cell type distributions
    - Identify batch effects or sample groupings
    - Understand categorical annotations

    Args:
        max_unique_values: Maximum number of unique values to consider a column
                          as categorical (default: 50). Columns with more unique
                          values are excluded to avoid high-cardinality fields.

    Returns:
        Formatted summary of categorical metadata columns with their unique values
        and counts.
    """
    adata = _get_adata()

    try:
        df = get_metadata_summary(adata, max_unique_values=max_unique_values)

        if df.empty:
            return "No categorical metadata columns found in adata.obs."

        # Format output
        result = f"Found {len(df)} categorical metadata columns in adata.obs:\n\n"

        for _, row in df.iterrows():
            col_name = row["column_name"]
            dtype = row["dtype"]
            n_unique = row["n_unique"]
            value_counts = row["value_counts"]

            result += f"Column: {col_name}\n"
            result += f"  Type: {dtype}\n"
            result += f"  Unique values: {n_unique}\n"
            result += f"  Distribution:\n"

            # Sort by count (descending)
            sorted_counts = sorted(value_counts.items(), key=lambda x: -x[1])
            for value, count in sorted_counts:
                percentage = (count / adata.n_obs) * 100
                result += f"    {value}: {count:,} cells ({percentage:.1f}%)\n"

            result += "\n"

        return result.strip()

    except Exception as e:
        logger.exception("Metadata inspection failed")
        return f"Error inspecting metadata: {e}"


@tool
def composition_analysis(
    row_key: str,
    col_key: str,
    show_percentages: bool = False,
    plot_only: bool = False,
) -> str:
    """Analyze cell composition across two categorical variables.

    This tool computes the cross-tabulation ONCE and returns BOTH:
    1. Exact count table (with CSV download)
    2. Stacked bar chart visualization

    This ensures the numbers in the table match the plot exactly, preventing
    any discrepancies or "cannot access data" errors.

    Use this for queries like:
    - "How many cells per cell type in each condition?"
    - "Show me the composition"
    - "Cell type distribution across samples"

    IMPORTANT: Do NOT call this if the composition table was already displayed
    in a previous turn. Instead, refer the user to the existing table.

    This is equivalent to Seurat's:
    obj@meta.data %>% group_by(row_key, col_key) %>% summarise(n=n())

    Args:
        row_key: Metadata column for rows (e.g., 'orig.ident', 'condition').
        col_key: Metadata column for columns (e.g., 'cell_type', 'leiden').
        show_percentages: If True, also show percentage table.
        plot_only: If True, skip table output (for plot-only requests).

    Returns:
        Both table and plot results with exact numbers visible to agent.
    """
    adata = _get_adata()

    try:
        # SINGLE SOURCE OF TRUTH - compute once
        crosstab_counts = cross_tabulate_metadata(adata, row_key, col_key, normalize=False)

        results = []

        # 1. Store and return count table (unless plot_only)
        if not plot_only:
            csv_buffer = io.StringIO()
            crosstab_counts.to_csv(csv_buffer)

            table_result = TableResult(
                csv_data=csv_buffer.getvalue(),
                code=f'composition_analysis(adata, "{row_key}", "{col_key}")',
                message=f"Cell counts: {col_key} across {row_key}",
                display_df=crosstab_counts.to_markdown(),
            )
            _store_table_and_return(table_result)

        # 2. Optionally show percentages
        if show_percentages and not plot_only:
            crosstab_pct = cross_tabulate_metadata(adata, row_key, col_key, normalize=True)
            csv_buffer_pct = io.StringIO()
            crosstab_pct.to_csv(csv_buffer_pct)

            table_result_pct = TableResult(
                csv_data=csv_buffer_pct.getvalue(),
                code=f'composition_analysis(adata, "{row_key}", "{col_key}", show_percentages=True)',
                message=f"Cell percentages: {col_key} across {row_key}",
                display_df=crosstab_pct.round(2).to_markdown(),
            )
            _store_table_and_return(table_result_pct)

        # 3. Generate plot from the SAME data
        plot_result = plot_composition(
            crosstab=crosstab_counts, row_key=row_key, col_key=col_key, kind="count"
        )
        _store_and_return(plot_result)

        # 4. Return a SHORT confirmation to the agent — no raw data dumps
        # The UI will render the table and plot separately
        output = (
            f"Composition analysis complete for {col_key} across {row_key}.\n"
            f"- Count table: displayed with CSV download button\n"
            f"- Stacked bar chart: displayed\n"
        )
        if show_percentages:
            output += "- Percentage table: displayed with CSV download button\n"
        output += (
            "\nThe exact cell counts are shown in the table above. "
            "If the user asks for exact numbers, refer them to the displayed table and CSV download. "
            "Do NOT regenerate the analysis — just reference the existing results."
        )

        update_viz_state("composition", row_key=row_key, col_key=col_key)
        return output

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Composition analysis failed")
        return f"Unexpected error: {e}"


@tool
def preprocess_data(
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mito: float = 20.0,
    n_top_genes: int = 2000,
    resolution: float = 0.5,
) -> str:
    """Run the full preprocessing pipeline: QC → filter → normalize → HVG → PCA → UMAP → clustering.

    IMPORTANT: Only use this on RAW count data. If data is already preprocessed, this will fail.

    Args:
        min_genes: Minimum genes per cell for filtering.
        min_cells: Minimum cells per gene for filtering.
        max_pct_mito: Maximum mitochondrial gene percentage.
        n_top_genes: Number of highly variable genes to select.
        resolution: Leiden clustering resolution.
    """
    global _adata  # noqa: PLW0603
    adata = _get_adata()

    # Check if data is already preprocessed
    if _dataset_state is not None:
        if _dataset_state.has_umap and _dataset_state.has_clustering:
            return (
                f"Error: Data is already preprocessed!\n\n"
                f"Current state:\n{_dataset_state.summary()}\n\n"
                f"The dataset already has UMAP coordinates and clustering ('{_dataset_state.cluster_key}').\n"
                f"You can directly use visualization tools like umap_plot, violin_plot, etc.\n"
                f"No need to preprocess again."
            )

    try:
        new_adata, result = run_preprocessing(
            adata,
            min_genes=min_genes,
            min_cells=min_cells,
            max_pct_mito=max_pct_mito,
            n_top_genes=n_top_genes,
            resolution=resolution,
        )
        _adata = new_adata
        if _adata_replaced_callback:
            _adata_replaced_callback(new_adata)
        _update_state()
        return result.message
    except Exception as e:
        logger.exception("Preprocessing failed")
        return f"Preprocessing error: {e}"


@tool
def differential_expression(groupby: str = "", method: str = "wilcoxon", n_genes: int = 20) -> str:
    """Run marker gene analysis to find distinguishing genes for ALL clusters/cell types.

    This performs one-vs-rest analysis for ALL groups (like Seurat's FindAllMarkers).
    Use this when you want to identify what genes define each cell type/cluster.

    IMPORTANT - When to use this tool:
    - "Find marker genes for all clusters"
    - "Full analysis report for all clusters"
    - "Differential expression for the entire dataset"
    - "What genes distinguish each cluster?"

    IMPORTANT - When NOT to use this tool:
    - For ONE specific cluster: use get_cluster_degs instead
    - For comparing TWO clusters: use compare_groups_de instead

    IMPORTANT - n_genes parameter:
    - Default is 20 (top 20 marker genes per cluster).
    - When user asks for "all markers", "complete marker table", "all genes per cluster",
      or "show all markers for each cluster", pass n_genes=0 to compute ALL markers.

    Args:
        groupby: The observation key for grouping. If empty, auto-detects clustering key.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
        n_genes: Number of top genes per cluster. Default 20. Use 0 for all genes.

    Returns:
        Summary message. Results are stored and can be accessed via get_top_markers
        or get_de_results_table.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        result = run_differential_expression(adata, groupby=groupby, method=method, n_genes=n_genes, target_group=None)
        _update_state()
        return result.message
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Marker gene analysis failed")
        return f"Marker gene analysis error: {e}"


@tool
def get_cluster_degs(
    cluster: str,
    groupby: str = "",
    method: str = "wilcoxon",
    n_genes: int = 20,
) -> str:
    """Get differentially expressed genes for ONE specific cluster (one-vs-rest).

    This performs one-vs-rest analysis for a SINGLE cluster/group, identifying genes
    that distinguish this cluster from all other clusters combined.

    IMPORTANT - When to use this tool:
    - "DEGs for Cluster 3"
    - "Marker genes for Cluster 4"
    - "What genes define B cells?"
    - "Differentially expressed genes for CD4 T cells"
    - "Full analysis report for cluster 6"

    IMPORTANT - When NOT to use this tool:
    - For ALL clusters: use differential_expression instead
    - For comparing TWO clusters: use compare_groups_de instead

    IMPORTANT - n_genes parameter:
    - Default is 20 (top 20 marker genes).
    - When user asks for "all markers", "all DEGs", or "complete marker list",
      pass n_genes=0 to compute ALL markers for this cluster.

    Args:
        cluster: The cluster/group identifier (e.g., "3", "Cluster 4", "B cells", "CD4 T cells").
        groupby: Optional observation key for grouping. Leave empty for auto-detection.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
        n_genes: Number of top genes. Default 20. Use 0 for all genes.

    Returns:
        Summary message with DEG results. Results are stored and can be accessed
        via get_de_results_table.
    """
    from src.analysis.cluster_resolution import resolve_analysis_scope

    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        # Resolve the cluster identifier
        resolved_cluster, resolved_groupby = resolve_analysis_scope(adata, cluster, groupby)

        if resolved_cluster is None:
            return "Error: Cluster identifier resolved to 'all clusters'. Use differential_expression for all clusters."

        # Run one-vs-rest DE for this specific cluster
        result = run_differential_expression(
            adata,
            groupby=resolved_groupby,
            method=method,
            n_genes=n_genes,
            target_group=resolved_cluster,
        )
        _update_state()
        return result.message
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Cluster DEG analysis failed")
        return f"Cluster DEG analysis error: {e}"


@tool
def compare_groups_de(
    group1: str,
    group2: str,
    groupby: str = "",
    method: str = "wilcoxon"
) -> str:
    """Run pairwise differential expression analysis between two specific groups.

    This performs direct comparison between two cell types/clusters (like Seurat's
    FindMarkers with ident.1 and ident.2). Use this when the user explicitly asks
    to compare two specific groups.

    IMPORTANT - Intelligent Cluster Resolution:
    This tool automatically handles BOTH numeric cluster IDs and cell type annotations:
    - Numeric IDs: "0", "1", "Cluster 2", "cluster_5" -> resolves to cluster column
    - Cell types: "CD4 T cells", "B cells" -> resolves to annotation column
    - You do NOT need to specify groupby - it auto-detects the correct column

    Examples of when to use this:
    - "Compare CD4 T cells vs CD8 T cells"
    - "Find DEGs between cluster 0 and cluster 1"
    - "Compare Cluster 1 vs Cluster 5"
    - "What genes differ between B cells and T cells"

    Args:
        group1: First group identifier (e.g., "CD4 T cells", "0", "Cluster 1").
        group2: Second group identifier (e.g., "CD8 T cells", "1", "Cluster 5").
        groupby: Optional observation key for grouping. Leave empty for auto-detection.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').

    Returns:
        Summary message with comparison results. Positive log2FC means higher in group1,
        negative log2FC means higher in group2.
    """
    adata = _get_adata()

    # Pass groupby as None if empty string (for auto-detection)
    groupby_param = groupby if groupby else None

    try:
        result = run_pairwise_de(
            adata,
            group1=group1,
            group2=group2,
            groupby=groupby_param,
            method=method
        )
        # Store results in uns for potential table export
        adata.uns["pairwise_de_result"] = {
            "group1": group1,
            "group2": group2,
            "results_df": result.results_df
        }
        return result.message
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Pairwise DE analysis failed")
        return f"Pairwise DE error: {e}"


@tool
def get_top_markers(n_genes_per_cluster: int = 10, groupby: str = "") -> str:
    """Get top marker genes for each cluster from DE results.

    Use this when you need to get marker genes for visualization (e.g., heatmap).
    Returns a comma-separated list of unique genes.

    Args:
        n_genes_per_cluster: Number of top genes to get from each cluster. Defaults to 10.
        groupby: The observation key for grouping. If empty, auto-detects clustering key.

    Returns:
        Comma-separated list of top marker genes (duplicates removed).
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        from src.analysis.marker_genes import get_top_marker_genes_per_cluster_exact

        cluster_genes = get_top_marker_genes_per_cluster_exact(
            adata, n_genes=n_genes_per_cluster, groupby=groupby
        )
        genes = get_top_marker_genes_per_cluster(adata, n_genes=n_genes_per_cluster, groupby=groupby)

        # Show per-cluster breakdown so the agent knows ALL genes must be used
        result = f"Top {n_genes_per_cluster} marker genes per cluster ({len(cluster_genes)} clusters):\n\n"
        for cluster, cluster_gene_list in cluster_genes.items():
            result += f"  {cluster}: {', '.join(cluster_gene_list)}\n"
        result += f"\nTotal unique genes: {len(genes)}\n"
        result += "IMPORTANT: Use ALL genes below for dotplot/heatmap (not just one cluster's genes):\n"
        result += ",".join(genes)

        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Get top markers failed")
        return f"Unexpected error: {e}"


@tool
def calculate_mito_pct() -> str:
    """Calculate mitochondrial percentage for all cells if not already present.

    This identifies genes starting with "MT-" and calculates the percentage of
    counts from these genes for each cell. Useful for quality control and filtering.

    Use this when users ask about "mito percentage", "mitochondrial percentage",
    or "pct_counts_mt" for violin plots or analysis.

    Returns:
        Summary message about the calculation.
    """
    adata = _get_adata()

    try:
        if "pct_counts_mt" in adata.obs.columns:
            n_cells = adata.n_obs
            mean_pct = adata.obs["pct_counts_mt"].mean()
            return (
                f"Mitochondrial percentage already calculated.\n"
                f"Cells: {n_cells}\n"
                f"Mean mito %: {mean_pct:.2f}%\n"
                f"Range: {adata.obs['pct_counts_mt'].min():.2f}% - {adata.obs['pct_counts_mt'].max():.2f}%\n\n"
                f"You can now use 'pct_counts_mt' in violin plots or other visualizations."
            )

        calculate_mito_percentage(adata)

        n_cells = adata.n_obs
        mean_pct = adata.obs["pct_counts_mt"].mean()
        return (
            f"Mitochondrial percentage calculated successfully.\n"
            f"Cells: {n_cells}\n"
            f"Mean mito %: {mean_pct:.2f}%\n"
            f"Range: {adata.obs['pct_counts_mt'].min():.2f}% - {adata.obs['pct_counts_mt'].max():.2f}%\n\n"
            f"You can now use 'pct_counts_mt' in violin plots or other visualizations."
        )
    except Exception as e:
        logger.exception("Calculate mito percentage failed")
        return f"Error calculating mitochondrial percentage: {e}"


@tool
def summarize_obs_column(column_name: str) -> str:
    """Compute descriptive statistics for a numeric column in adata.obs.

    Use this tool when users ask for summaries of QC metrics or cell-level metadata,
    such as:
    - "Summarize total counts"
    - "What's the median gene count per cell?"
    - "Show statistics for n_genes_by_counts"
    - "Describe the UMI distribution"

    This tool handles common naming variations automatically:
    - total_counts, n_counts, nCount_RNA (for UMI/read depth)
    - n_genes_by_counts, n_genes, nFeature_RNA (for gene complexity)
    - pct_counts_mt, percent.mt (for mitochondrial percentage)

    Args:
        column_name: The column name or metric name to summarize.
                    Can be the exact column name or a common alias.

    Returns:
        Formatted summary with descriptive statistics:
        - Number of cells (non-missing values)
        - Mean, median, standard deviation
        - Min, max, 25th and 75th percentiles
    """
    adata = _get_adata()

    try:
        # Try to resolve the column name
        resolved_col = resolve_qc_metric_column(adata, column_name)

        if resolved_col is None:
            # Column not found, provide helpful error
            numeric_cols = [
                col for col in adata.obs.columns
                if pd.api.types.is_numeric_dtype(adata.obs[col])
            ]
            return (
                f"Error: Could not find column '{column_name}' in adata.obs.\n\n"
                f"Available numeric columns:\n" + "\n".join(f"  - {col}" for col in numeric_cols)
            )

        # Compute statistics
        stats = get_obs_column_statistics(adata, resolved_col)

        # Format output
        result = f"Summary statistics for '{resolved_col}':\n\n"

        if resolved_col != column_name:
            result += f"(Resolved from '{column_name}')\n\n"

        result += f"Cells (non-missing): {stats['n_cells']:,}\n"
        if stats['n_missing'] > 0:
            result += f"Missing values: {stats['n_missing']:,}\n"

        result += f"\nDescriptive Statistics:\n"
        result += f"  Mean:   {stats['mean']:.2f}\n"
        result += f"  Median: {stats['median']:.2f}\n"
        result += f"  Std:    {stats['std']:.2f}\n"
        result += f"  Min:    {stats['min']:.2f}\n"
        result += f"  25%:    {stats['q25']:.2f}\n"
        result += f"  75%:    {stats['q75']:.2f}\n"
        result += f"  Max:    {stats['max']:.2f}\n"

        return result

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Summarize obs column failed")
        return f"Unexpected error: {e}"


@tool
def summarize_qc_metrics_tool() -> str:
    """Summarize all common QC metrics in the dataset.

    Use this tool when users ask for:
    - "Show QC summary"
    - "Summarize quality control metrics"
    - "What are the QC statistics?"
    - "Show me total counts and gene counts"

    Automatically detects and summarizes common QC metrics:
    - Total counts / UMI depth (total_counts, n_counts, nCount_RNA)
    - Number of genes expressed (n_genes_by_counts, n_genes, nFeature_RNA)
    - Mitochondrial percentage (pct_counts_mt, percent.mt)

    Returns:
        Formatted table with statistics for all detected QC metrics.
    """
    adata = _get_adata()

    try:
        df = summarize_qc_metrics(adata)

        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)

        table = TableResult(
            csv_data=csv_buffer.getvalue(),
            code="summarize_qc_metrics_tool()",
            message=f"QC metrics summary for {adata.n_obs:,} cells ({len(df)} metrics detected).",
            display_df=df.to_markdown(index=False),
        )
        return _store_table_and_return(table)

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Summarize QC metrics failed")
        return f"Unexpected error: {e}"


@tool
def get_de_results_table(
    groupby: str = "",
    top_n_per_cluster: int = 0,
    target_cluster: str = ""
) -> str:
    """Generate a comprehensive differential expression results TABLE with full statistics.

    USE THIS TOOL when users ask for:
    - "show the differential expression table"
    - "export DE results"
    - "download differential expression results"
    - "show DE statistics"
    - "differential expression results table"
    - "full analysis report for cluster X"
    - "marker table for all clusters"

    This generates a COMPLETE TABLE with ALL statistical information:
    - cluster: cluster ID
    - gene: gene name
    - scores: statistical scores
    - log2fc: log2 fold change (if available)
    - pval: p-value (if available)
    - pval_adj: adjusted p-value (if available)

    The table will be displayed in the UI with a CSV download button.

    IMPORTANT - Scope Handling:
    - If target_cluster is empty: exports ALL clusters
    - If target_cluster is specified: exports ONLY that cluster

    Args:
        groupby: The observation key for grouping. If empty, auto-detects clustering key.
        top_n_per_cluster: If > 0, only return top N genes per cluster. If 0, return all genes.
        target_cluster: If specified, only export this cluster. If empty, export all clusters.

    Returns:
        Summary message. The full table will be available for download in the UI.
    """
    from src.analysis.cluster_resolution import resolve_analysis_scope

    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        # Resolve target cluster if specified
        resolved_target = None
        if target_cluster:
            resolved_target, groupby = resolve_analysis_scope(adata, target_cluster, groupby)

        # Get all DE results
        df = get_all_de_results(adata)

        # Filter to specific cluster if requested
        if resolved_target is not None:
            df = df[df["cluster"] == resolved_target].copy()
            if len(df) == 0:
                return (
                    f"Error: No DE results found for cluster '{target_cluster}' "
                    f"(resolved to '{resolved_target}'). "
                    f"Please run get_cluster_degs first."
                )

        # Filter to top N per cluster if requested
        if top_n_per_cluster > 0:
            df = df.groupby("cluster").head(top_n_per_cluster).reset_index(drop=True)

        # Format the dataframe for display
        df_display = df.copy()
        if "log2fc" in df_display.columns:
            df_display["log2fc"] = df_display["log2fc"].round(3)
        if "pval" in df_display.columns:
            df_display["pval"] = df_display["pval"].apply(lambda x: f"{x:.2e}")
        if "pval_adj" in df_display.columns:
            df_display["pval_adj"] = df_display["pval_adj"].apply(lambda x: f"{x:.2e}")
        if "scores" in df_display.columns:
            df_display["scores"] = df_display["scores"].round(3)

        # Generate CSV data
        csv_data = df.to_csv(index=False)

        # Generate display representation (first 100 rows)
        display_rows = min(100, len(df_display))
        display_df_subset = df_display.head(display_rows)

        # Convert to markdown with proper formatting
        display_df = display_df_subset.to_markdown(index=False, tablefmt="pipe")

        # Add summary info at the top
        if len(df) > display_rows:
            display_df = f"**Showing first {display_rows} of {len(df)} total rows**\n\n{display_df}"

        # Create code representation
        scope_desc = f"cluster {target_cluster}" if target_cluster else "all clusters"
        code = f"""# Differential expression results for {scope_desc}
# Grouping by: {groupby}
# Total genes: {len(df)}
# Clusters: {df['cluster'].nunique()}

import pandas as pd
de_results = pd.read_csv('de_results.csv')
print(de_results.head())
"""

        # Create message
        n_clusters = df["cluster"].nunique()
        cluster_names = sorted(df["cluster"].unique())
        columns_list = ", ".join(df.columns)

        if resolved_target is not None:
            message = (
                f"Differential expression results table generated for {target_cluster} "
                f"(resolved to '{resolved_target}').\n"
                f"Grouping: {groupby}\n"
                f"Total genes: {len(df)}\n"
                f"{'Top ' + str(top_n_per_cluster) + ' genes' if top_n_per_cluster > 0 else 'All genes included'}\n\n"
                f"Columns: {columns_list}\n"
                f"The full table is available for download as CSV."
            )
        else:
            message = (
                f"Differential expression results table generated for ALL clusters.\n"
                f"Grouping: {groupby}\n"
                f"Clusters: {n_clusters} ({', '.join(str(c) for c in cluster_names[:5])}"
                f"{'...' if n_clusters > 5 else ''})\n"
                f"Total genes: {len(df)}\n"
                f"{'Top ' + str(top_n_per_cluster) + ' genes per cluster' if top_n_per_cluster > 0 else 'All genes included'}\n\n"
                f"Columns: {columns_list}\n"
                f"The full table is available for download as CSV."
            )

        # Store the table result
        table_result = TableResult(
            csv_data=csv_data,
            code=code,
            message=message,
            display_df=display_df
        )
        return _store_table_and_return(table_result)

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Get DE results table failed")
        return f"Unexpected error: {e}"


@tool
def get_pairwise_de_table(top_n: int = 0, padj_threshold: float = 0.05) -> str:
    """Generate a table for the most recent pairwise differential expression comparison.

    USE THIS TOOL after running compare_groups_de when users want to see or export
    the detailed results table.

    Args:
        top_n: If > 0, only return top N genes (sorted by adjusted p-value). If 0, return all.
        padj_threshold: Only include genes with adjusted p-value below this threshold.

    Returns:
        Summary message. The full table will be available for download in the UI.
    """
    adata = _get_adata()

    try:
        # Check if pairwise DE results exist
        if "pairwise_de_result" not in adata.uns:
            return (
                "Error: No pairwise DE results found. "
                "Please run compare_groups_de first to compare two groups."
            )

        pairwise_result = adata.uns["pairwise_de_result"]
        group1 = pairwise_result["group1"]
        group2 = pairwise_result["group2"]
        df = pairwise_result["results_df"].copy()

        # Filter by adjusted p-value threshold
        df = df[df["pval_adj"] < padj_threshold]

        # Filter to top N if requested
        if top_n > 0:
            df = df.head(top_n)

        if len(df) == 0:
            return (
                f"No genes found with adjusted p-value < {padj_threshold}. "
                f"Try increasing the threshold or check if the comparison was successful."
            )

        # Format the dataframe for display
        df_display = df.copy()
        df_display["log2fc"] = df_display["log2fc"].round(3)
        df_display["pval"] = df_display["pval"].apply(lambda x: f"{x:.2e}")
        df_display["pval_adj"] = df_display["pval_adj"].apply(lambda x: f"{x:.2e}")

        # Generate CSV data
        csv_data = df.to_csv(index=False)

        # Generate display representation (first 100 rows)
        display_rows = min(100, len(df_display))
        display_df_subset = df_display.head(display_rows)

        # Convert to markdown with proper formatting
        from tabulate import tabulate
        display_df = tabulate(display_df_subset, headers="keys", tablefmt="pipe", showindex=False)

        # Add summary info at the top
        if len(df) > display_rows:
            display_df = f"**Showing first {display_rows} of {len(df)} total rows**\n\n{display_df}"

        # Create code representation
        code = f"""# Pairwise differential expression: {group1} vs {group2}
# Total significant genes (padj < {padj_threshold}): {len(df)}
# Positive log2FC: higher in {group1}
# Negative log2FC: higher in {group2}

import pandas as pd
de_results = pd.read_csv('pairwise_de_{group1}_vs_{group2}.csv')
print(de_results.head())
"""

        # Create message
        n_upregulated = (df["log2fc"] > 0).sum()
        n_downregulated = (df["log2fc"] < 0).sum()
        message = (
            f"Pairwise DE results table: {group1} vs {group2}\n"
            f"Significant genes (padj < {padj_threshold}): {len(df)}\n"
            f"  - Higher in {group1}: {n_upregulated} genes\n"
            f"  - Higher in {group2}: {n_downregulated} genes\n\n"
            f"{'Top ' + str(top_n) + ' genes shown' if top_n > 0 else 'All significant genes included'}\n"
            f"The full table is available for download as CSV."
        )

        # Store the table result
        table_result = TableResult(
            csv_data=csv_data,
            code=code,
            message=message,
            display_df=display_df
        )
        return _store_table_and_return(table_result)

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Get pairwise DE table failed")
        return f"Unexpected error: {e}"


@tool
def get_cluster_mapping(groupby: str = "") -> str:
    """Get the numeric index to cluster/cell type name mapping.

    Use this tool when you need to understand which numeric cluster ID corresponds
    to which cell type name. This is especially useful when the dataset has
    annotated cell types instead of numeric cluster IDs.

    Examples of when to use this:
    - User asks "which cell type is cluster 1?"
    - User asks "what are the cluster names?"
    - Before comparing clusters by numeric ID (e.g., "compare cluster 1 vs cluster 5")

    Args:
        groupby: The observation key for grouping. If empty, auto-detects clustering key.

    Returns:
        A formatted mapping showing numeric index -> cluster/cell type name.
    """
    from src.analysis.cluster_resolution import create_cluster_index_mapping

    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        index_mapping = create_cluster_index_mapping(adata, groupby)

        # Format the mapping nicely
        mapping_lines = [f"Cluster index mapping for '{groupby}':"]
        for idx, name in index_mapping.items():
            mapping_lines.append(f"  {idx}: {name}")

        return "\n".join(mapping_lines)

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Get cluster mapping failed")
        return f"Unexpected error: {e}"


@tool
def lookup_gene(gene_name: str) -> str:
    """Look up a gene name in the loaded dataset, resolving typos, case mismatches, and species prefixes.

    Use this tool when:
    - A gene is not found by other tools
    - The user may have entered a protein name, alias, or wrong-species name
    - You want to confirm a gene exists before plotting

    If resolved, returns the exact gene name to use and proceeds.
    If unresolved, returns top fuzzy candidates for the user to confirm.

    Args:
        gene_name: The gene name to look up (e.g. 'TNNT2', 'Tnnt2', 'tnnt2').
    """
    adata = _get_adata()
    result = lookup_gene_name(adata, gene_name)
    return result.message


def get_all_tools() -> list:
    """Return all tools for the agent."""
    return [
        # Plotting tools
        umap_plot, violin_plot, dotplot, dotplot_combined, dotplot_matrix, feature_plot,
        heatmap_plot, scatter_plot, volcano_plot_tool,
        # Core tools
        dataset_info, check_data_status, inspect_metadata, composition_analysis,
        preprocess_data, differential_expression, get_cluster_degs, compare_groups_de, get_top_markers,
        # QC and statistics tools
        calculate_mito_pct, summarize_obs_column, summarize_qc_metrics_tool,
        # Cluster resolution tools
        get_cluster_mapping,
        # Gene lookup tool
        lookup_gene,
        # Table tools
        get_de_results_table, get_pairwise_de_table,
        # Analysis/reasoning tools
        *analysis_tools.get_analysis_tools(),
        # Subset analysis tools
        *subset_tools.get_subset_tools(),
    ]
