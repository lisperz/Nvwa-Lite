"""Input validation helpers for plotting functions.

Provides gene and observation key validation with fuzzy matching suggestions.
"""

from __future__ import annotations

import difflib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def gene_exists(adata: AnnData, gene: str) -> bool:
    """Check if a gene exists in adata.var_names or adata.raw.var_names."""
    if gene in adata.var_names:
        return True
    if adata.raw is not None and gene in adata.raw.var_names:
        return True
    return False


def all_gene_names(adata: AnnData) -> list[str]:
    """Get all gene names from both processed and raw layers."""
    names = set(adata.var_names)
    if adata.raw is not None:
        names.update(adata.raw.var_names)
    return sorted(names)


def validate_gene(adata: AnnData, gene: str) -> str | None:
    """Return an error message if gene is not in the dataset, else None."""
    if gene_exists(adata, gene):
        return None
    all_genes = all_gene_names(adata)
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


def validate_obs_or_gene(adata: AnnData, field: str) -> str | None:
    """Return an error message if field is not in adata.obs or adata.var_names, else None.

    This is useful for plotting functions that accept both QC metrics (from obs)
    and gene names (from var_names).

    Args:
        adata: The annotated data matrix.
        field: Field name to validate (can be obs column or gene name).

    Returns:
        None if field exists, error message otherwise.
    """
    # Check if it's an observation column (QC metric)
    if field in adata.obs.columns:
        return None

    # Check if it's a gene name
    if gene_exists(adata, field):
        return None

    # Not found - provide helpful error
    gene_err = validate_gene(adata, field)
    obs_err = validate_obs_key(adata, field)

    return f"'{field}' not found as gene or observation key. {gene_err} {obs_err}"
