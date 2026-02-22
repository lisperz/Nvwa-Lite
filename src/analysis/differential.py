"""Differential expression analysis wrappers.

Provides a clean interface around scanpy's rank_genes_groups.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


@dataclass
class DEResult:
    """Result of a differential expression analysis."""

    group: str
    reference: str
    results_df: pd.DataFrame
    message: str


def run_differential_expression(
    adata: AnnData,
    groupby: str = "leiden",
    *,
    method: str = "wilcoxon",
    n_genes: int = 100,
) -> DEResult:
    """Run DE analysis using rank_genes_groups.

    Args:
        adata: Processed AnnData with clustering.
        groupby: Observation key for grouping.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
        n_genes: Number of top genes per group.
    """
    if groupby not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(f"Key '{groupby}' not found. Available: {available}")

    logger.info("Running DE: groupby=%s, method=%s", groupby, method)
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    message = (
        f"Differential expression complete.\n"
        f"Method: {method} | Groups: {len(groups)} ({', '.join(groups[:5])}"
        f"{'...' if len(groups) > 5 else ''})\n"
        f"Top {n_genes} genes per group stored."
    )

    df = get_de_dataframe(adata, group=groups[0])
    return DEResult(group="all", reference="rest", results_df=df, message=message)


def get_de_dataframe(adata: AnnData, group: str) -> pd.DataFrame:
    """Extract DE results for a specific group as a tidy DataFrame.

    Returns DataFrame with columns: gene, log2fc, pval, pval_adj.
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("No DE results found. Run differential_expression first.")

    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)
    if group not in groups:
        raise ValueError(
            f"Group '{group}' not in DE results. Available: {', '.join(groups)}"
        )

    return pd.DataFrame({
        "gene": result["names"][group],
        "log2fc": result["logfoldchanges"][group],
        "pval": result["pvals"][group],
        "pval_adj": result["pvals_adj"][group],
    })
