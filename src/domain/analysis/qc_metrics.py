"""QC metrics calculation and summarization utilities.

Provides functions to compute descriptive statistics for cell-level QC metrics
stored in adata.obs columns.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


def resolve_qc_metric_column(
    adata: AnnData,
    metric_name: str,
) -> str | None:
    """Resolve a QC metric name to the actual column name in adata.obs.

    Handles common naming variations across different preprocessing pipelines.

    Args:
        adata: The annotated data matrix.
        metric_name: The metric name to resolve (e.g., "total_counts", "n_genes").

    Returns:
        The actual column name if found, None otherwise.
    """
    metric_lower = metric_name.lower().replace(" ", "_")

    # Direct match
    if metric_name in adata.obs.columns:
        return metric_name

    # Common aliases for total counts / UMI depth
    if any(x in metric_lower for x in ["total_count", "n_count", "ncount", "umi", "depth"]):
        for candidate in ["total_counts", "n_counts", "nCount_RNA", "nUMI", "n_umi"]:
            if candidate in adata.obs.columns:
                return candidate

    # Common aliases for number of genes expressed
    if any(x in metric_lower for x in ["n_gene", "ngene", "complexity", "feature"]):
        for candidate in ["n_genes_by_counts", "n_genes", "nFeature_RNA", "nGene"]:
            if candidate in adata.obs.columns:
                return candidate

    # Common aliases for mitochondrial percentage
    if any(x in metric_lower for x in ["mito", "mt", "percent_mt"]):
        for candidate in ["pct_counts_mt", "percent.mt", "percent_mito", "pct_mt"]:
            if candidate in adata.obs.columns:
                return candidate

    return None


def get_obs_column_statistics(
    adata: AnnData,
    column_name: str,
) -> dict[str, float | int]:
    """Compute descriptive statistics for a numeric column in adata.obs.

    Args:
        adata: The annotated data matrix.
        column_name: The column name in adata.obs.

    Returns:
        Dictionary with statistics: n_cells, mean, median, std, min, max, q25, q75.

    Raises:
        ValueError: If column not found or not numeric.
    """
    if column_name not in adata.obs.columns:
        raise ValueError(f"Column '{column_name}' not found in adata.obs")

    data = adata.obs[column_name]

    # Check if numeric
    if not pd.api.types.is_numeric_dtype(data):
        raise ValueError(f"Column '{column_name}' is not numeric (dtype: {data.dtype})")

    # Remove NaN values for statistics
    data_clean = data.dropna()

    if len(data_clean) == 0:
        raise ValueError(f"Column '{column_name}' contains only NaN values")

    stats = {
        "n_cells": int(len(data_clean)),
        "n_missing": int(len(data) - len(data_clean)),
        "mean": float(data_clean.mean()),
        "median": float(data_clean.median()),
        "std": float(data_clean.std()),
        "min": float(data_clean.min()),
        "max": float(data_clean.max()),
        "q25": float(data_clean.quantile(0.25)),
        "q75": float(data_clean.quantile(0.75)),
    }

    return stats


def summarize_qc_metrics(
    adata: AnnData,
    metric_names: list[str] | None = None,
    groupby: str | None = None,
) -> pd.DataFrame:
    """Summarize multiple QC metrics with descriptive statistics.

    Args:
        adata: The annotated data matrix.
        metric_names: List of metric names to summarize. If None, auto-detects
                     common QC metrics (total_counts, n_genes, pct_counts_mt).
        groupby: Optional obs column to group by. When set, returns one row per
                 (group × metric) with a leading 'group' column.

    Returns:
        DataFrame with one row per metric (global) or per group × metric (grouped).

    Raises:
        ValueError: If groupby column is not found in adata.obs.
    """
    if groupby is not None and groupby not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(
            f"Column '{groupby}' not found in adata.obs. Available columns: {available}"
        )

    if metric_names is None:
        # Auto-detect common QC metrics
        metric_names = []
        for candidate in ["total_counts", "n_counts", "nCount_RNA"]:
            if candidate in adata.obs.columns:
                metric_names.append(candidate)
                break

        for candidate in ["n_genes_by_counts", "n_genes", "nFeature_RNA"]:
            if candidate in adata.obs.columns:
                metric_names.append(candidate)
                break

        for candidate in ["pct_counts_mt", "percent.mt", "percent_mito"]:
            if candidate in adata.obs.columns:
                metric_names.append(candidate)
                break

    if not metric_names:
        raise ValueError("No QC metrics found in adata.obs")

    if groupby is None:
        return _summarize_subset(adata, metric_names)

    groups = sorted(adata.obs[groupby].dropna().unique().tolist(), key=str)
    group_frames: list[pd.DataFrame] = []
    for group_val in groups:
        mask = adata.obs[groupby] == group_val
        subset = adata[mask]
        group_df = _summarize_subset(subset, metric_names)
        group_df.insert(0, "group", str(group_val))
        group_frames.append(group_df)

    if not group_frames:
        raise ValueError(f"No groups found in column '{groupby}'")

    return pd.concat(group_frames, ignore_index=True)


def _summarize_subset(adata: AnnData, metric_names: list[str]) -> pd.DataFrame:
    """Compute per-metric statistics for a single (possibly subsetted) AnnData."""
    results = []
    for metric in metric_names:
        resolved = resolve_qc_metric_column(adata, metric)
        if resolved is None:
            logger.warning(f"Could not resolve metric '{metric}', skipping")
            continue

        try:
            stats = get_obs_column_statistics(adata, resolved)
            stats["metric"] = resolved
            results.append(stats)
        except ValueError as e:
            logger.warning(f"Could not compute statistics for '{resolved}': {e}")
            continue

    if not results:
        raise ValueError("No valid QC metrics could be summarized")

    df = pd.DataFrame(results)
    cols = ["metric"] + [c for c in df.columns if c != "metric"]
    return df[cols]
