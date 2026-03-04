"""Basic plotting functions: UMAP, violin, and feature plots.

Each function validates inputs, generates a plot, and returns a PlotResult
containing the PNG image bytes, the code that was executed, and a description.
"""

from __future__ import annotations

import io
import logging
from dataclasses import dataclass

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from src.plotting.validation import gene_exists, validate_gene, validate_obs_key

logger = logging.getLogger(__name__)


@dataclass
class PlotResult:
    """Bundle of plot output: image bytes, source code, and description."""

    image: bytes
    code: str
    message: str


def _figure_to_bytes() -> bytes:
    """Capture the current matplotlib figure as PNG bytes and close it."""
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    plt.close("all")
    buf.seek(0)
    return buf.read()


def plot_umap(
    adata: AnnData,
    color: str,
    title: str | None = None,
    show_labels: bool = False,
    show_legend: bool = True,
    split_by: str | None = None,
) -> PlotResult:
    """Generate a UMAP plot colored by an obs key or gene.

    Args:
        adata: The annotated data matrix.
        color: An observation column name or gene name to color by.
        title: Optional plot title.
        show_labels: Whether to show cluster labels on the plot.
        show_legend: Whether to show the legend.
        split_by: Optional observation key to split the plot into separate panels.

    Returns:
        PlotResult with image, code, and description.
    """
    # Validate: color can be an obs key OR a gene name (including raw)
    is_obs = color in adata.obs.columns
    is_gene = gene_exists(adata, color)
    if not is_obs and not is_gene:
        # Try to give a helpful error
        err = validate_gene(adata, color)
        if err is None:
            err = validate_obs_key(adata, color)
        raise ValueError(err or f"'{color}' not found as gene or observation key.")

    # Validate split_by if provided
    if split_by:
        split_err = validate_obs_key(adata, split_by)
        if split_err:
            raise ValueError(split_err)

    # Build code string
    code = f'sc.pl.umap(adata, color="{color}"'
    if title:
        code += f', title="{title}"'
    if show_labels:
        code += ", label=True"
    if not show_legend:
        code += ", legend_loc=None"
    if split_by:
        code += f', groups="{split_by}"'
    code += ")"

    logger.info("Executing: %s", code)

    # Generate plot with split if requested
    if split_by:
        # Create separate panels for each group in split_by
        import numpy as np
        groups = adata.obs[split_by].unique()
        n_groups = len(groups)

        # Calculate grid layout
        n_cols = min(4, n_groups)  # Max 4 columns
        n_rows = (n_groups + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
        if n_groups == 1:
            axes = [axes]
        else:
            axes = axes.flatten() if n_groups > 1 else [axes]

        # Plot each group in a separate panel
        for idx, group in enumerate(sorted(groups)):
            ax = axes[idx]

            # Create a mask for this group
            mask = adata.obs[split_by] == group
            subset = adata[mask, :]

            # Plot UMAP for this subset
            sc.pl.umap(
                subset,
                color=color,
                ax=ax,
                title=f"{split_by}: {group}",
                legend_loc=None,
                show=False,
            )

        # Hide unused subplots
        for idx in range(n_groups, len(axes)):
            axes[idx].axis('off')

        plt.tight_layout()
    else:
        # Generate regular plot
        legend_loc = "right margin" if show_legend else None
        sc.pl.umap(
            adata,
            color=color,
            title=title,
            legend_loc=legend_loc,
            add_outline=show_labels,
            show=False,
        )

        # Add labels if requested
        if show_labels and is_obs:
            # Get cluster centers for labeling
            import numpy as np
            umap_coords = adata.obsm["X_umap"]
            clusters = adata.obs[color].values
            for cluster in np.unique(clusters):
                mask = clusters == cluster
                center_x = np.median(umap_coords[mask, 0])
                center_y = np.median(umap_coords[mask, 1])
                plt.text(center_x, center_y, str(cluster),
                        fontsize=12, fontweight='bold',
                        ha='center', va='center')

    image = _figure_to_bytes()

    label = "cell type" if color == "louvain" else color
    if split_by:
        message = f"UMAP plot colored by {label}, split by {split_by} into {len(groups)} panels."
    else:
        message = f"UMAP plot colored by {label}."
    return PlotResult(image=image, code=code, message=message)


def plot_violin(
    adata: AnnData,
    genes: str | list[str],
    groupby: str = "louvain",
    title: str | None = None,
) -> PlotResult:
    """Generate a violin plot for one or more genes across cell groups.

    Args:
        adata: The annotated data matrix.
        genes: Gene name(s) to plot. Can be a single gene or list of genes.
               Also supports QC metrics from adata.obs (e.g., 'pct_counts_mt', 'n_genes_by_counts').
        groupby: Observation key to group cells by.
        title: Optional plot title.

    Returns:
        PlotResult with image, code, and description.
    """
    # Handle both single gene and list of genes
    if isinstance(genes, str):
        gene_list = [genes]
    else:
        gene_list = genes

    # Validate all genes - check both var_names and obs.columns for QC metrics
    for gene in gene_list:
        # First check if it's a QC metric in obs.columns
        if gene in adata.obs.columns:
            continue  # Valid QC metric, skip gene validation
        # Otherwise validate as a gene
        gene_err = validate_gene(adata, gene)
        if gene_err:
            raise ValueError(gene_err)

    obs_err = validate_obs_key(adata, groupby)
    if obs_err:
        raise ValueError(obs_err)

    # Build code string
    if len(gene_list) == 1:
        code = f'sc.pl.violin(adata, keys="{gene_list[0]}", groupby="{groupby}"'
    else:
        code = f'sc.pl.violin(adata, keys={gene_list}, groupby="{groupby}"'
    if title:
        code += f', title="{title}"'
    code += ")"

    logger.info("Executing: %s", code)
    sc.pl.violin(adata, keys=gene_list, groupby=groupby, rotation=45, show=False)
    if title:
        plt.title(title)
    image = _figure_to_bytes()

    if len(gene_list) == 1:
        message = f"Violin plot of {gene_list[0]} expression across {groupby} groups."
    else:
        message = f"Violin plot of {', '.join(gene_list)} expression across {groupby} groups."
    return PlotResult(image=image, code=code, message=message)


def plot_feature(adata: AnnData, gene: str) -> PlotResult:
    """Generate a feature plot (UMAP colored by gene expression with viridis).

    Args:
        adata: The annotated data matrix (must have UMAP coordinates).
        gene: Gene name to visualize.
    """
    gene_err = validate_gene(adata, gene)
    if gene_err:
        raise ValueError(gene_err)

    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP not computed. Run preprocessing first.")

    code = f'sc.pl.umap(adata, color="{gene}", cmap="viridis")'

    logger.info("Executing: %s", code)
    sc.pl.umap(adata, color=gene, cmap="viridis", show=False)
    image = _figure_to_bytes()

    message = f"Feature plot of {gene} expression on UMAP."
    return PlotResult(image=image, code=code, message=message)
