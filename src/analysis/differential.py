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
    n_genes: int = 20,
) -> DEResult:
    """Run marker gene analysis across all groups (one-vs-rest for each group).

    This is equivalent to Seurat's FindAllMarkers(). Use this when you want to
    identify marker genes that distinguish each cell type/cluster from all others.

    Args:
        adata: Processed AnnData with clustering.
        groupby: Observation key for grouping.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
        n_genes: Number of top genes per group. Defaults to 20.
    """
    if groupby not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(f"Key '{groupby}' not found. Available: {available}")

    logger.info("Running marker gene analysis: groupby=%s, method=%s, n_genes=%d", groupby, method, n_genes)
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    message = (
        f"Marker gene analysis complete (one-vs-rest for all groups).\n"
        f"Method: {method} | Groups: {len(groups)} ({', '.join(groups[:5])}"
        f"{'...' if len(groups) > 5 else ''})\n"
        f"Top {n_genes} marker genes per group stored."
    )

    df = get_de_dataframe(adata, group=groups[0])
    return DEResult(group="all", reference="rest", results_df=df, message=message)


def run_pairwise_de(
    adata: AnnData,
    group1: str,
    group2: str,
    groupby: str = "leiden",
    *,
    method: str = "wilcoxon",
) -> DEResult:
    """Run pairwise differential expression between two specific groups.

    This is equivalent to Seurat's FindMarkers(ident.1=group1, ident.2=group2).
    Use this when you want to compare two specific cell types/clusters directly.

    Args:
        adata: Processed AnnData with clustering.
        group1: First group identifier (e.g., "CD4 T cells" or "0").
        group2: Second group identifier (e.g., "CD8 T cells" or "1").
        groupby: Observation key for grouping.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').

    Returns:
        DEResult with comparison results (positive log2FC = higher in group1).
    """
    if groupby not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(f"Key '{groupby}' not found. Available: {available}")

    # Get unique values in the groupby column
    available_groups = adata.obs[groupby].unique().tolist()
    available_groups_str = [str(g) for g in available_groups]

    # Validate group1 and group2
    if str(group1) not in available_groups_str:
        raise ValueError(
            f"Group '{group1}' not found in '{groupby}'. "
            f"Available groups: {', '.join(available_groups_str)}"
        )
    if str(group2) not in available_groups_str:
        raise ValueError(
            f"Group '{group2}' not found in '{groupby}'. "
            f"Available groups: {', '.join(available_groups_str)}"
        )

    logger.info(
        "Running pairwise DE: %s vs %s (groupby=%s, method=%s)",
        group1, group2, groupby, method
    )

    # Run rank_genes_groups with specific reference
    # This compares group1 vs group2 directly
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[str(group1)],
        reference=str(group2),
        method=method,
    )

    # Extract results
    result = adata.uns["rank_genes_groups"]
    genes = result["names"][str(group1)]
    logfc = result["logfoldchanges"][str(group1)]
    pvals = result["pvals"][str(group1)]
    pvals_adj = result["pvals_adj"][str(group1)]

    # Create DataFrame
    df = pd.DataFrame({
        "gene": genes,
        "log2fc": logfc,
        "pval": pvals,
        "pval_adj": pvals_adj,
    })

    # Sort by adjusted p-value
    df = df.sort_values("pval_adj")

    message = (
        f"Pairwise differential expression complete.\n"
        f"Comparison: {group1} vs {group2}\n"
        f"Method: {method}\n"
        f"Total genes analyzed: {len(df)}\n"
        f"Significant genes (padj < 0.05): {(df['pval_adj'] < 0.05).sum()}\n\n"
        f"Interpretation:\n"
        f"- Positive log2FC: Higher expression in {group1}\n"
        f"- Negative log2FC: Higher expression in {group2}"
    )

    return DEResult(group=str(group1), reference=str(group2), results_df=df, message=message)


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


def get_all_de_results(adata: AnnData) -> pd.DataFrame:
    """Extract all DE results for all groups as a comprehensive DataFrame.

    Returns DataFrame with columns: cluster, gene, and available statistics
    (log2fc, pval, pval_adj, scores - depending on what's available).
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("No DE results found. Run differential_expression first.")

    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)

    # Check which fields are available
    has_logfc = "logfoldchanges" in result
    has_pval = "pvals" in result
    has_pval_adj = "pvals_adj" in result
    has_scores = "scores" in result

    # Collect all results
    all_results = []
    for group in groups:
        data = {
            "cluster": group,
            "gene": result["names"][group],
        }

        # Add optional fields if available
        if has_scores:
            data["scores"] = result["scores"][group]
        if has_logfc:
            data["log2fc"] = result["logfoldchanges"][group]
        if has_pval:
            data["pval"] = result["pvals"][group]
        if has_pval_adj:
            data["pval_adj"] = result["pvals_adj"][group]

        df = pd.DataFrame(data)
        all_results.append(df)

    # Combine all groups
    combined_df = pd.concat(all_results, ignore_index=True)

    # Sort by cluster and then by best available metric
    if has_pval_adj:
        combined_df = combined_df.sort_values(["cluster", "pval_adj"])
    elif has_scores:
        combined_df = combined_df.sort_values(["cluster", "scores"], ascending=[True, False])
    else:
        combined_df = combined_df.sort_values("cluster")

    return combined_df
