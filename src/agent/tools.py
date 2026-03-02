"""LangChain tool definitions for the visualization agent.

Each tool wraps a plotting/analysis function, providing validation,
error handling, and structured output for the agent.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Callable

from langchain_core.tools import tool

from src.agent import analysis_tools
from src.analysis.calculations import calculate_mito_percentage
from src.analysis.differential import get_de_dataframe, run_differential_expression
from src.analysis.marker_genes import get_top_marker_genes_per_cluster
from src.analysis.preprocessing import run_preprocessing
from src.plotting.comparison import plot_dotplot, plot_heatmap, plot_scatter
from src.plotting.executor import PlotResult, plot_feature, plot_umap, plot_violin
from src.plotting.volcano import plot_volcano
from src.types import DatasetState, detect_dataset_state

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

_adata: AnnData | None = None
_dataset_state: DatasetState | None = None
_plot_results: list[PlotResult] = []
_adata_replaced_callback: Callable[[AnnData], None] | None = None


def bind_dataset(adata: AnnData) -> None:
    """Bind a dataset so tools can access it."""
    global _adata  # noqa: PLW0603
    _adata = adata
    # Also bind to analysis_tools module
    analysis_tools._adata = adata


def bind_dataset_state(state: DatasetState) -> None:
    """Bind the dataset state for tools to reference."""
    global _dataset_state  # noqa: PLW0603
    _dataset_state = state


def set_adata_replaced_callback(cb: Callable[[AnnData], None]) -> None:
    """Register a callback invoked when preprocessing replaces adata."""
    global _adata_replaced_callback  # noqa: PLW0603
    _adata_replaced_callback = cb


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
    return f"Plot generated successfully.\nCode: {result.code}\n{result.message}"


def get_plot_results() -> list[PlotResult]:
    """Retrieve all plot results from this turn (called by UI)."""
    # Combine results from both modules
    all_results = list(_plot_results)
    all_results.extend(analysis_tools._plot_results)
    return all_results


def clear_plot_results() -> None:
    """Clear the plot result buffer."""
    _plot_results.clear()
    analysis_tools._plot_results.clear()


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
        return _store_and_return(plot_umap(
            adata, color=color_by, show_labels=show_labels, show_legend=show_legend,
            split_by=split_by if split_by else None
        ))
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
        return _store_and_return(plot_violin(adata, genes=gene_list, groupby=groupby))
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
        return _store_and_return(plot_dotplot(adata, genes=gene_list, groupby=groupby))
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Dot plot failed")
        return f"Unexpected error: {e}"


@tool
def feature_plot(gene: str) -> str:
    """Generate a feature plot showing gene expression on UMAP (viridis colormap).

    Args:
        gene: The gene name to visualize.
    """
    adata = _get_adata()
    try:
        return _store_and_return(plot_feature(adata, gene=gene))
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
        return _store_and_return(plot_heatmap(
            adata, genes=gene_list, groupby=groupby, n_genes_per_cluster=n_genes_per_cluster
        ))
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
        return _store_and_return(plot_scatter(adata, gene_x=gene_x, gene_y=gene_y, color_by=color))
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
def differential_expression(groupby: str = "", method: str = "wilcoxon") -> str:
    """Run differential expression analysis to find marker genes per cluster.

    Args:
        groupby: The observation key for grouping. If empty, auto-detects clustering key.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        result = run_differential_expression(adata, groupby=groupby, method=method)
        _update_state()
        return result.message
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("DE analysis failed")
        return f"DE error: {e}"


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
        genes = get_top_marker_genes_per_cluster(adata, n_genes=n_genes_per_cluster, groupby=groupby)

        result = f"Top {n_genes_per_cluster} marker genes per cluster:\n"
        result += f"Total unique genes: {len(genes)}\n\n"
        result += "Genes (comma-separated for heatmap):\n"
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


def get_all_tools() -> list:
    """Return all tools for the agent."""
    return [
        # Plotting tools
        umap_plot, violin_plot, dotplot, feature_plot,
        heatmap_plot, scatter_plot, volcano_plot_tool,
        # Core tools
        dataset_info, check_data_status,
        preprocess_data, differential_expression, get_top_markers, calculate_mito_pct,
        # Analysis/reasoning tools
        *analysis_tools.get_analysis_tools(),
    ]
