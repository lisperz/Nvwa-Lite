"""Scanpy plotting executor with input validation.

Each function validates inputs, generates a plot, and returns a PlotResult
containing the PNG image bytes, the code that was executed, and a description.
"""

from __future__ import annotations

import difflib
import io
import logging
from dataclasses import dataclass

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


@dataclass
class PlotResult:
    """Bundle of plot output: image bytes, source code, and description."""

    image: bytes
    code: str
    message: str


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def _gene_exists(adata: AnnData, gene: str) -> bool:
    """Check if a gene exists in adata.var_names or adata.raw.var_names."""
    if gene in adata.var_names:
        return True
    if adata.raw is not None and gene in adata.raw.var_names:
        return True
    return False


def _all_gene_names(adata: AnnData) -> list[str]:
    """Get all gene names from both processed and raw layers."""
    names = set(adata.var_names)
    if adata.raw is not None:
        names.update(adata.raw.var_names)
    return sorted(names)


def validate_gene(adata: AnnData, gene: str) -> str | None:
    """Return an error message if gene is not in the dataset, else None."""
    if _gene_exists(adata, gene):
        return None
    all_genes = _all_gene_names(adata)
    matches = difflib.get_close_matches(gene, all_genes, n=5, cutoff=0.5)
    if matches:
        suggestions = ", ".join(matches)
        return f"Gene '{gene}' not found. Did you mean: {suggestions}?"
    return f"Gene '{gene}' not found in this dataset ({len(all_genes)} genes available)."


def validate_obs_key(adata: AnnData, key: str) -> str | None:
    """Return an error message if key is not in adata.obs, else None."""
    if key in adata.obs.columns:
        return None
    available = ", ".join(sorted(adata.obs.columns))
    return f"Observation key '{key}' not found. Available keys: {available}."


def _figure_to_bytes() -> bytes:
    """Capture the current matplotlib figure as PNG bytes and close it."""
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    plt.close("all")
    buf.seek(0)
    return buf.read()


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_umap(adata: AnnData, color: str, title: str | None = None) -> PlotResult:
    """Generate a UMAP plot colored by an obs key or gene.

    Args:
        adata: The annotated data matrix.
        color: An observation column name or gene name to color by.
        title: Optional plot title.

    Returns:
        PlotResult with image, code, and description.
    """
    # Validate: color can be an obs key OR a gene name (including raw)
    is_obs = color in adata.obs.columns
    is_gene = _gene_exists(adata, color)
    if not is_obs and not is_gene:
        # Try to give a helpful error
        err = validate_gene(adata, color)
        if err is None:
            err = validate_obs_key(adata, color)
        raise ValueError(err or f"'{color}' not found as gene or observation key.")

    code = f'sc.pl.umap(adata, color="{color}"'
    if title:
        code += f', title="{title}"'
    code += ")"

    logger.info("Executing: %s", code)
    sc.pl.umap(adata, color=color, title=title, show=False)
    image = _figure_to_bytes()

    label = "cell type" if color == "louvain" else color
    message = f"UMAP plot colored by {label}."
    return PlotResult(image=image, code=code, message=message)


def plot_violin(
    adata: AnnData,
    gene: str,
    groupby: str = "louvain",
    title: str | None = None,
) -> PlotResult:
    """Generate a violin plot for a gene across cell groups.

    Args:
        adata: The annotated data matrix.
        gene: Gene name to plot.
        groupby: Observation key to group cells by.
        title: Optional plot title.

    Returns:
        PlotResult with image, code, and description.
    """
    gene_err = validate_gene(adata, gene)
    if gene_err:
        raise ValueError(gene_err)
    obs_err = validate_obs_key(adata, groupby)
    if obs_err:
        raise ValueError(obs_err)

    code = f'sc.pl.violin(adata, keys="{gene}", groupby="{groupby}"'
    if title:
        code += f', title="{title}"'
    code += ")"

    logger.info("Executing: %s", code)
    sc.pl.violin(adata, keys=gene, groupby=groupby, rotation=45, show=False)
    if title:
        plt.title(title)
    image = _figure_to_bytes()

    message = f"Violin plot of {gene} expression across {groupby} groups."
    return PlotResult(image=image, code=code, message=message)


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
