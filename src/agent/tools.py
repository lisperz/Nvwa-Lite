"""LangChain tool definitions for the visualization agent.

Each tool wraps a plotting/analysis function, providing validation,
error handling, and structured output for the agent.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Callable

from langchain_core.tools import tool

from src.analysis.differential import get_de_dataframe, run_differential_expression
from src.analysis.preprocessing import run_preprocessing
from src.plotting.executor import (
    PlotResult,
    plot_dotplot,
    plot_feature,
    plot_heatmap,
    plot_umap,
    plot_violin,
)
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
    return list(_plot_results)


def clear_plot_results() -> None:
    """Clear the plot result buffer."""
    _plot_results.clear()


@tool
def umap_plot(color_by: str) -> str:
    """Generate a UMAP plot colored by an observation key or gene name.

    Args:
        color_by: The observation column (e.g. 'louvain') or gene name (e.g. 'CD3E').
    """
    adata = _get_adata()
    try:
        return _store_and_return(plot_umap(adata, color=color_by))
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("UMAP plot failed")
        return f"Unexpected error: {e}"


@tool
def violin_plot(gene: str, groupby: str = "louvain") -> str:
    """Generate a violin plot showing a gene's expression distribution across groups.

    Args:
        gene: The gene name to plot (e.g. 'MS4A1').
        groupby: The observation key to group cells by. Defaults to 'louvain'.
    """
    adata = _get_adata()
    try:
        return _store_and_return(plot_violin(adata, gene=gene, groupby=groupby))
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Violin plot failed")
        return f"Unexpected error: {e}"


@tool
def dotplot(genes: str, groupby: str = "louvain") -> str:
    """Generate a dot plot for one or more genes across cell groups.

    Args:
        genes: Comma-separated gene names (e.g. 'CD3E,MS4A1,NKG7').
        groupby: The observation key to group cells by. Defaults to 'louvain'.
    """
    adata = _get_adata()
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
def heatmap_plot(genes: str, groupby: str = "leiden") -> str:
    """Generate a heatmap for genes across cell groups.

    Args:
        genes: Comma-separated gene names (e.g. 'CD3E,MS4A1,NKG7,LYZ').
        groupby: The observation key to group cells by. Defaults to 'leiden'.
    """
    adata = _get_adata()
    gene_list = [g.strip() for g in genes.split(",")]
    try:
        return _store_and_return(plot_heatmap(adata, genes=gene_list, groupby=groupby))
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Heatmap plot failed")
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

    Args:
        min_genes: Minimum genes per cell for filtering.
        min_cells: Minimum cells per gene for filtering.
        max_pct_mito: Maximum mitochondrial gene percentage.
        n_top_genes: Number of highly variable genes to select.
        resolution: Leiden clustering resolution.
    """
    global _adata  # noqa: PLW0603
    adata = _get_adata()
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
def differential_expression(groupby: str = "leiden", method: str = "wilcoxon") -> str:
    """Run differential expression analysis to find marker genes per cluster.

    Args:
        groupby: The observation key for grouping (e.g. 'leiden').
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
    """
    adata = _get_adata()
    try:
        result = run_differential_expression(adata, groupby=groupby, method=method)
        _update_state()
        return result.message
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("DE analysis failed")
        return f"DE error: {e}"


def get_all_tools() -> list:
    """Return all tools for the agent."""
    return [
        umap_plot, violin_plot, dotplot, feature_plot,
        heatmap_plot, volcano_plot_tool,
        dataset_info, check_data_status,
        preprocess_data, differential_expression,
    ]
