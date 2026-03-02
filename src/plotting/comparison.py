"""Comparison plots: dotplot, heatmap, and scatter plots.

These plots compare gene expression across multiple genes or cell groups.
"""

from __future__ import annotations

import io
import logging
from dataclasses import dataclass

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from src.plotting.validation import validate_gene, validate_obs_key

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


def plot_dotplot(
    adata: AnnData,
    genes: list[str],
    groupby: str = "louvain",
) -> PlotResult:
    """Generate a dot plot for one or more genes across cell groups.

    Args:
        adata: The annotated data matrix.
        genes: List of gene names.
        groupby: Observation key to group cells by.

    Returns:
        PlotResult with image, code, and description.
    """
    errors: list[str] = []
    for g in genes:
        err = validate_gene(adata, g)
        if err:
            errors.append(err)
    if errors:
        raise ValueError(" ".join(errors))

    obs_err = validate_obs_key(adata, groupby)
    if obs_err:
        raise ValueError(obs_err)

    genes_str = str(genes)
    code = f'sc.pl.dotplot(adata, var_names={genes_str}, groupby="{groupby}")'

    logger.info("Executing: %s", code)
    sc.pl.dotplot(adata, var_names=genes, groupby=groupby, show=False)
    image = _figure_to_bytes()

    message = f"Dot plot of {', '.join(genes)} across {groupby} groups."
    return PlotResult(image=image, code=code, message=message)


def plot_heatmap(
    adata: AnnData,
    genes: list[str],
    groupby: str = "leiden",
) -> PlotResult:
    """Generate a heatmap for genes across cell groups.

    Args:
        adata: The annotated data matrix.
        genes: List of gene names.
        groupby: Observation key to group cells by.
    """
    errors: list[str] = []
    for g in genes:
        err = validate_gene(adata, g)
        if err:
            errors.append(err)
    if errors:
        raise ValueError(" ".join(errors))

    obs_err = validate_obs_key(adata, groupby)
    if obs_err:
        raise ValueError(obs_err)

    genes_str = str(genes)
    code = f'sc.pl.heatmap(adata, var_names={genes_str}, groupby="{groupby}", swap_axes=True)'

    logger.info("Executing: %s", code)
    # swap_axes=True puts groups on y-axis for better label readability
    # figsize adjusted for long group names
    n_groups = adata.obs[groupby].nunique()
    fig_height = max(6, n_groups * 0.5)
    fig_width = max(10, len(genes) * 1.2)
    sc.pl.heatmap(
        adata, var_names=genes, groupby=groupby,
        swap_axes=True, figsize=(fig_width, fig_height),
        dendrogram=True, show=False,
    )
    image = _figure_to_bytes()

    message = f"Heatmap of {', '.join(genes)} across {groupby} groups."
    return PlotResult(image=image, code=code, message=message)


def plot_scatter(
    adata: AnnData,
    gene_x: str,
    gene_y: str,
    color_by: str | None = None,
) -> PlotResult:
    """Generate a scatter plot showing correlation between two genes.

    Args:
        adata: The annotated data matrix.
        gene_x: Gene name for x-axis.
        gene_y: Gene name for y-axis.
        color_by: Optional observation key or gene to color points by.

    Returns:
        PlotResult with image, code, and description.
    """
    # Validate genes
    err_x = validate_gene(adata, gene_x)
    if err_x:
        raise ValueError(err_x)
    err_y = validate_gene(adata, gene_y)
    if err_y:
        raise ValueError(err_y)

    # Validate color_by if provided
    if color_by:
        from src.plotting.validation import gene_exists
        is_obs = color_by in adata.obs.columns
        is_gene = gene_exists(adata, color_by)
        if not is_obs and not is_gene:
            err = validate_gene(adata, color_by)
            if err is None:
                err = validate_obs_key(adata, color_by)
            raise ValueError(err or f"'{color_by}' not found as gene or observation key.")

    # Build code string
    code = f'sc.pl.scatter(adata, x="{gene_x}", y="{gene_y}"'
    if color_by:
        code += f', color="{color_by}"'
    code += ")"

    logger.info("Executing: %s", code)
    sc.pl.scatter(adata, x=gene_x, y=gene_y, color=color_by, show=False)
    image = _figure_to_bytes()

    if color_by:
        message = f"Scatter plot of {gene_x} vs {gene_y}, colored by {color_by}."
    else:
        message = f"Scatter plot showing correlation between {gene_x} and {gene_y}."
    return PlotResult(image=image, code=code, message=message)
