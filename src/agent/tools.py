"""LangChain tool definitions for the visualization agent.

Each tool wraps a plotting executor function, providing validation,
error handling, and structured output for the agent.
"""

from __future__ import annotations

import base64
import logging
from typing import TYPE_CHECKING

from langchain_core.tools import tool

from src.plotting.executor import (
    PlotResult,
    plot_dotplot,
    plot_umap,
    plot_violin,
)

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

# Module-level reference to the active dataset, set by bind_dataset()
_adata: AnnData | None = None


def bind_dataset(adata: AnnData) -> None:
    """Bind a dataset so tools can access it."""
    global _adata  # noqa: PLW0603
    _adata = adata


def _get_adata() -> AnnData:
    if _adata is None:
        raise RuntimeError("No dataset loaded. Please upload or select a dataset first.")
    return _adata


def _result_to_dict(result: PlotResult) -> str:
    """Convert PlotResult to a string response for the agent.

    The image bytes are stored in session state for the UI to display.
    """
    # Store image in a module-level buffer for the UI to pick up
    global _last_plot_result  # noqa: PLW0603
    _last_plot_result = result
    return f"Plot generated successfully.\nCode: {result.code}\n{result.message}"


# Module-level buffer for the last plot result
_last_plot_result: PlotResult | None = None


def get_last_plot_result() -> PlotResult | None:
    """Retrieve the last plot result (called by the UI layer)."""
    return _last_plot_result


def clear_last_plot_result() -> None:
    """Clear the last plot result buffer."""
    global _last_plot_result  # noqa: PLW0603
    _last_plot_result = None


# ---------------------------------------------------------------------------
# Tool definitions
# ---------------------------------------------------------------------------

@tool
def umap_plot(color_by: str) -> str:
    """Generate a UMAP plot colored by an observation key or gene name.

    Args:
        color_by: The observation column (e.g. 'louvain') or gene name (e.g. 'CD3E') to color the UMAP by.
    """
    adata = _get_adata()
    try:
        result = plot_umap(adata, color=color_by)
        return _result_to_dict(result)
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("UMAP plot failed")
        return f"Unexpected error generating UMAP plot: {e}"


@tool
def violin_plot(gene: str, groupby: str = "louvain") -> str:
    """Generate a violin plot showing a gene's expression distribution across cell groups.

    Args:
        gene: The gene name to plot (e.g. 'MS4A1', 'CD3E').
        groupby: The observation key to group cells by. Defaults to 'louvain'.
    """
    adata = _get_adata()
    try:
        result = plot_violin(adata, gene=gene, groupby=groupby)
        return _result_to_dict(result)
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Violin plot failed")
        return f"Unexpected error generating violin plot: {e}"


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
        result = plot_dotplot(adata, genes=gene_list, groupby=groupby)
        return _result_to_dict(result)
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Dot plot failed")
        return f"Unexpected error generating dot plot: {e}"


@tool
def dataset_info() -> str:
    """Get dataset metadata: available genes, observation keys, and cell count."""
    adata = _get_adata()
    obs_keys = ", ".join(sorted(adata.obs.columns))

    # Include genes from both processed and raw layers
    gene_set: set[str] = set(adata.var_names)
    if adata.raw is not None:
        gene_set.update(adata.raw.var_names)
    all_genes = sorted(gene_set)

    return (
        f"Dataset summary:\n"
        f"- Cells: {adata.n_obs}\n"
        f"- Total genes: {len(all_genes)} (processed: {adata.n_vars}, raw: {adata.raw.n_vars if adata.raw else 0})\n"
        f"- Observation keys: {obs_keys}\n"
        f"- All genes: {', '.join(all_genes)}"
    )


def get_all_tools() -> list:
    """Return all tools for the agent."""
    return [umap_plot, violin_plot, dotplot, dataset_info]
