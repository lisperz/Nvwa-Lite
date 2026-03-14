"""Differential expression analysis wrappers.

Provides a clean interface around scanpy's rank_genes_groups.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd
import scanpy as sc
from anndata import AnnData

from src.analysis.cluster_resolution import resolve_pairwise_groups

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
    target_group: str | None = None,
) -> DEResult:
    """Run marker gene analysis (one-vs-rest).

    This is equivalent to Seurat's FindAllMarkers() or FindMarkers(ident.1=group).

    - If target_group is None: runs one-vs-rest for ALL groups
    - If target_group is specified: runs one-vs-rest for ONLY that group

    Args:
        adata: Processed AnnData with clustering.
        groupby: Observation key for grouping.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').
        n_genes: Number of top genes per group. Defaults to 20.
        target_group: If specified, only analyze this group vs rest. If None, analyze all groups.
    """
    if groupby not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(f"Key '{groupby}' not found. Available: {available}")

    # If target_group specified, run one-vs-rest for that group only
    if target_group is not None:
        available_groups = adata.obs[groupby].astype(str).unique().tolist()
        if target_group not in available_groups:
            raise ValueError(
                f"Group '{target_group}' not found in '{groupby}'. "
                f"Available groups: {', '.join(available_groups)}"
            )

        logger.info(
            "Running one-vs-rest DE for single group: %s (groupby=%s, method=%s, n_genes=%d)",
            target_group, groupby, method, n_genes
        )
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            groups=[target_group],
            reference="rest",
            method=method,
            n_genes=n_genes
        )

        df = get_de_dataframe(adata, group=target_group)
        message = (
            f"Marker gene analysis complete (one-vs-rest for {target_group}).\n"
            f"Method: {method} | Grouping: {groupby}\n"
            f"Top {n_genes} marker genes identified for {target_group} vs all other groups.\n"
            f"Significant genes (padj < 0.05): {(df['pval_adj'] < 0.05).sum()}"
        )
        return DEResult(group=target_group, reference="rest", results_df=df, message=message)

    # Otherwise, run one-vs-rest for ALL groups
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
    groupby: str | None = None,
    *,
    method: str = "wilcoxon",
) -> DEResult:
    """Run pairwise differential expression between two specific groups.

    This is equivalent to Seurat's FindMarkers(ident.1=group1, ident.2=group2).
    Use this when you want to compare two specific cell types/clusters directly.

    Intelligently handles both numeric cluster IDs (0, 1, "Cluster 2") and
    cell type annotations ("CD4 T cells", "B cells").

    Args:
        adata: Processed AnnData with clustering.
        group1: First group identifier (e.g., "CD4 T cells", "0", "Cluster 1").
        group2: Second group identifier (e.g., "CD8 T cells", "1", "Cluster 5").
        groupby: Observation key for grouping. If None, auto-detects.
        method: Statistical method ('wilcoxon', 't-test', 'logreg').

    Returns:
        DEResult with comparison results (positive log2FC = higher in group1).
    """
    # Resolve group identifiers intelligently
    resolved_group1, resolved_group2, resolved_groupby = resolve_pairwise_groups(
        adata, group1, group2, groupby
    )

    logger.info(
        "Running pairwise DE: %s vs %s (resolved to: %s vs %s in '%s', method=%s)",
        group1, group2, resolved_group1, resolved_group2, resolved_groupby, method
    )

    # Run rank_genes_groups with specific reference
    # This compares group1 vs group2 directly
    sc.tl.rank_genes_groups(
        adata,
        groupby=resolved_groupby,
        groups=[resolved_group1],
        reference=resolved_group2,
        method=method,
    )

    # Extract results
    result = adata.uns["rank_genes_groups"]
    genes = result["names"][resolved_group1]
    logfc = result["logfoldchanges"][resolved_group1]
    pvals = result["pvals"][resolved_group1]
    pvals_adj = result["pvals_adj"][resolved_group1]

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
        f"Grouping column: {resolved_groupby}\n"
        f"Resolved identifiers: {resolved_group1} vs {resolved_group2}\n"
        f"Method: {method}\n"
        f"Total genes analyzed: {len(df)}\n"
        f"Significant genes (padj < 0.05): {(df['pval_adj'] < 0.05).sum()}\n\n"
        f"Interpretation:\n"
        f"- Positive log2FC: Higher expression in {group1}\n"
        f"- Negative log2FC: Higher expression in {group2}"
    )

    return DEResult(group=resolved_group1, reference=resolved_group2, results_df=df, message=message)


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
