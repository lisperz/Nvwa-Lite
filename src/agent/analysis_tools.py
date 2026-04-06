"""Reasoning and analysis tools for the agent.

Provides tools for calculating expression, finding patterns, highlighting clusters,
and annotating cell types.
"""

from __future__ import annotations

import io
import logging
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from langchain_core.tools import tool

from src.analysis.calculations import (
    calculate_cluster_averages,
    find_top_expressing_cluster,
    get_cluster_cell_indices,
    rename_cluster_labels,
)
from src.plotting.executor import PlotResult

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

# Global state (shared with main tools.py)
_adata: AnnData | None = None
_plot_results: list[PlotResult] = []


def _get_adata() -> AnnData:
    """Get the bound dataset."""
    if _adata is None:
        raise RuntimeError("No dataset loaded.")
    return _adata


def _get_cluster_key() -> str:
    """Auto-detect the clustering key from the dataset."""
    adata = _get_adata()
    for key in ("leiden", "louvain", "seurat_clusters", "cluster", "clusters"):
        if key in adata.obs.columns:
            return key
    # If no clustering found, return leiden as default
    return "leiden"


def _figure_to_bytes() -> bytes:
    """Capture the current matplotlib figure as PNG bytes and close it."""
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    plt.close("all")
    buf.seek(0)
    return buf.read()


def _store_and_return(result: PlotResult) -> str:
    """Append a PlotResult to the buffer and return a text summary."""
    _plot_results.append(result)
    return f"Plot generated successfully.\nCode: {result.code}\n{result.message}"


@tool
def calculate_average_expression(gene: str, groupby: str = "") -> str:
    """Calculate average expression of a gene across all clusters.

    Use this to understand which clusters express a gene most highly.

    Args:
        gene: Gene name to analyze (e.g. 'LYZ', 'MS4A1').
        groupby: Observation key for grouping. If empty, auto-detects clustering key.

    Returns:
        Table showing mean expression and cell count per cluster, sorted by expression.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        df = calculate_cluster_averages(adata, gene, groupby)

        # Format as a readable table
        result = f"Average expression of {gene} across {groupby} clusters:\n\n"
        result += df.to_string(index=False)
        result += f"\n\nHighest expression: Cluster {df.iloc[0]['cluster']} "
        result += f"(mean={df.iloc[0]['mean_expression']:.2f})"

        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Calculate average expression failed")
        return f"Unexpected error: {e}"


@tool
def find_highest_expression(gene: str, groupby: str = "") -> str:
    """Find which group has the highest expression of a specific gene.

    CRITICAL: The groupby parameter determines the DIMENSION of your answer:
    - groupby="cell_type" -> answers "which cell type"
    - groupby="orig.ident" -> answers "which condition/sample"
    - groupby="leiden" -> answers "which cluster"

    IMPORTANT: Match the groupby to the user's question dimension!
    - User asks "which cell type" -> use cell_type column
    - User asks "which condition" -> use condition column
    - User asks "which cluster" -> use cluster column

    Args:
        gene: Gene name to analyze (e.g. 'MKI67', 'MS4A1').
        groupby: Observation key for grouping. MUST match the dimension user asked about.
                 If empty, auto-detects clustering key (WARNING: may not match user intent).

    Returns:
        The group with highest expression, mean value, and top 3 groups for context.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided (with warning)
    if not groupby:
        groupby = _get_cluster_key()
        logger.warning(f"No groupby specified for find_highest_expression, defaulting to {groupby}. This may not match user intent!")

    try:
        cluster_id, mean_expr = find_top_expressing_cluster(adata, gene, groupby)

        # Get cell count for context
        df = calculate_cluster_averages(adata, gene, groupby)
        top_row = df[df["cluster"] == cluster_id].iloc[0]
        cell_count = int(top_row["cell_count"])

        # Determine dimension type for answer formatting
        groupby_lower = groupby.lower()
        if "cell_type" in groupby_lower or "celltype" in groupby_lower or "annotation" in groupby_lower or "label" in groupby_lower:
            dimension_name = "cell type"
        elif "condition" in groupby_lower or "sample" in groupby_lower or "orig.ident" in groupby_lower or "batch" in groupby_lower or "treatment" in groupby_lower:
            dimension_name = "condition/sample"
        elif "cluster" in groupby_lower or "leiden" in groupby_lower or "louvain" in groupby_lower or "seurat" in groupby_lower:
            dimension_name = "cluster"
        else:
            dimension_name = "group"

        result = f"Highest {gene} expression by {dimension_name} (grouping column: {groupby}):\n\n"
        result += f"**ANSWER: {cluster_id}** has the highest expression.\n"
        result += f"Mean expression: {mean_expr:.2f}\n"
        result += f"Cells in this {dimension_name}: {cell_count:,}\n\n"
        result += f"Top 3 {dimension_name}s by {gene} expression:\n"
        result += df.head(3).to_string(index=False)
        result += f"\n\n=== IMPORTANT FOR YOUR RESPONSE ===\n"
        result += f"Your answer dimension is '{dimension_name}' (column: {groupby}).\n"
        result += f"If user asked 'which {dimension_name}', answer with: '{cluster_id}'\n"
        result += f"Do NOT answer with a different dimension (e.g., cluster number if they asked for cell type).\n"
        result += f"If this dimension doesn't match what the user asked, re-run this tool with the correct groupby parameter."

        return result
    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Find highest expression failed")
        return f"Unexpected error: {e}"


@tool
def highlight_cluster(cluster_id: str, color_by: str = "") -> str:
    """Generate a UMAP plot highlighting a specific cluster.

    Use this to visually emphasize a cluster of interest on the UMAP.

    Args:
        cluster_id: The cluster to highlight (e.g. '0', '3').
        color_by: Observation key for coloring. If empty, auto-detects clustering key.

    Returns:
        UMAP plot with the specified cluster highlighted.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not color_by:
        color_by = _get_cluster_key()

    try:
        if "X_umap" not in adata.obsm:
            return "Error: UMAP not computed. Run preprocessing first."

        if color_by not in adata.obs.columns:
            return f"Error: '{color_by}' not found in observations."

        # Get cell indices for the cluster
        cell_indices = get_cluster_cell_indices(adata, cluster_id, color_by)

        if len(cell_indices) == 0:
            return f"Error: Cluster '{cluster_id}' not found in '{color_by}'."

        # Create highlight mask
        highlight = np.zeros(adata.n_obs, dtype=bool)
        highlight[cell_indices] = True

        # Generate plot
        code = f'sc.pl.umap(adata, color="{color_by}", highlight=[cells_in_cluster_{cluster_id}])'
        logger.info("Executing: %s", code)

        # Plot with highlighting
        _fig, ax = plt.subplots(figsize=(8, 6))

        # Plot all cells in gray
        umap_coords = adata.obsm["X_umap"]
        ax.scatter(umap_coords[~highlight, 0], umap_coords[~highlight, 1],
                  c='lightgray', s=1, alpha=0.5, label='Other cells')

        # Plot highlighted cluster in red
        ax.scatter(umap_coords[highlight, 0], umap_coords[highlight, 1],
                  c='red', s=3, alpha=0.8, label=f'Cluster {cluster_id}')

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title(f'UMAP highlighting Cluster {cluster_id}')
        ax.legend()

        image = _figure_to_bytes()

        message = f"UMAP plot highlighting cluster {cluster_id} ({len(cell_indices):,} cells)."
        result = PlotResult(image=image, code=code, message=message)

        return _store_and_return(result)

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Highlight cluster failed")
        return f"Unexpected error: {e}"


@tool
def rename_cluster(old_name: str, new_name: str, groupby: str = "") -> str:
    """Rename a cluster to annotate it with a cell type name.

    Use this after identifying a cluster's cell type based on marker genes.
    This modifies the dataset in-place.

    Args:
        old_name: Current cluster ID (e.g. '0', '3').
        new_name: New name for the cluster (e.g. 'B cells', 'T cells').
        groupby: Observation key for grouping. If empty, auto-detects clustering key.

    Returns:
        Confirmation message with the renamed cluster.
    """
    adata = _get_adata()

    # Auto-detect clustering key if not provided
    if not groupby:
        groupby = _get_cluster_key()

    try:
        # Verify cluster exists
        if groupby not in adata.obs.columns:
            return f"Error: '{groupby}' not found in observations."

        clusters = adata.obs[groupby].astype(str).values
        if old_name not in clusters:
            available = ", ".join(sorted(set(clusters)))
            return f"Error: Cluster '{old_name}' not found. Available: {available}"

        # Count cells before rename
        cell_count = np.sum(clusters == old_name)

        # Perform rename
        rename_cluster_labels(adata, old_name, new_name, groupby)

        result = f"Successfully renamed cluster '{old_name}' to '{new_name}'.\n"
        result += f"Affected {cell_count:,} cells in '{groupby}'.\n\n"
        result += "Tip: Use umap_plot to visualize the updated labels."

        return result

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Rename cluster failed")
        return f"Unexpected error: {e}"


def get_analysis_tools() -> list:
    """Return all analysis/reasoning tools."""
    return [
        calculate_average_expression,
        find_highest_expression,
        highlight_cluster,
        rename_cluster,
    ]
