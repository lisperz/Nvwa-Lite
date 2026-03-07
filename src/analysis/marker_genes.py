"""Helper functions for extracting marker genes from DE results.

Provides utilities to get top marker genes per cluster or overall.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


def get_top_marker_genes_per_cluster(
    adata: AnnData,
    n_genes: int = 20,
    groupby: str = "leiden",
) -> list[str]:
    """Get top N marker genes for each cluster, removing duplicates.

    Args:
        adata: AnnData with DE results.
        n_genes: Number of top genes per cluster.
        groupby: Grouping key used in DE analysis.

    Returns:
        List of unique gene names, preserving order.
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("No DE results found. Run differential_expression first.")

    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)

    # Collect top N genes from each cluster
    all_genes = []
    for group in groups:
        genes = result["names"][group][:n_genes]
        all_genes.extend(genes)

    # Remove duplicates while preserving order
    unique_genes = list(dict.fromkeys(all_genes))

    logger.info(
        "Collected %d genes from %d clusters, %d unique",
        len(all_genes), len(groups), len(unique_genes)
    )

    return unique_genes


def get_top_marker_genes_overall(
    adata: AnnData,
    n_genes: int = 20,
) -> list[str]:
    """Get top N marker genes overall, ranked by significance across all clusters.

    Args:
        adata: AnnData with DE results.
        n_genes: Number of top genes to return.

    Returns:
        List of top N gene names.
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("No DE results found. Run differential_expression first.")

    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)

    # Collect all genes with their scores
    gene_scores = {}
    for group in groups:
        genes = result["names"][group]
        scores = result["scores"][group] if "scores" in result else result["logfoldchanges"][group]

        for gene, score in zip(genes, scores):
            if gene not in gene_scores:
                gene_scores[gene] = []
            gene_scores[gene].append(abs(float(score)))

    # Rank genes by maximum score across clusters
    gene_max_scores = {gene: max(scores) for gene, scores in gene_scores.items()}
    sorted_genes = sorted(gene_max_scores.items(), key=lambda x: -x[1])

    top_genes = [gene for gene, _score in sorted_genes[:n_genes]]

    logger.info("Selected top %d genes from %d total unique genes", n_genes, len(gene_scores))

    return top_genes


def get_top_marker_genes_per_cluster_exact(
    adata: AnnData,
    n_genes: int = 20,
    groupby: str = "leiden",
) -> dict[str, list[str]]:
    """Get exactly N marker genes for each cluster (may include duplicates across clusters).

    Args:
        adata: AnnData with DE results.
        n_genes: Number of genes per cluster.
        groupby: Grouping key used in DE analysis.

    Returns:
        Dictionary mapping cluster ID to list of top genes.
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("No DE results found. Run differential_expression first.")

    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)

    cluster_genes = {}
    for group in groups:
        genes = result["names"][group][:n_genes]
        cluster_genes[str(group)] = list(genes)

    return cluster_genes
