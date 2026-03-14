"""Calculation utilities for cluster analysis and reasoning.

Provides functions to calculate average expression, find top expressing clusters,
and manipulate cluster annotations.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


def calculate_cluster_averages(
    adata: AnnData,
    gene: str,
    groupby: str = "leiden",
) -> pd.DataFrame:
    """Calculate average expression of a gene across clusters.

    Args:
        adata: The annotated data matrix.
        gene: Gene name to calculate averages for.
        groupby: Observation key for grouping (e.g. 'leiden').

    Returns:
        DataFrame with columns: cluster, mean_expression, cell_count
    """
    if gene not in adata.var_names and (adata.raw is None or gene not in adata.raw.var_names):
        raise ValueError(f"Gene '{gene}' not found in dataset.")

    if groupby not in adata.obs.columns:
        raise ValueError(f"Grouping key '{groupby}' not found in observations.")

    # Get expression data (prefer raw if available for genes)
    if adata.raw is not None and gene in adata.raw.var_names:
        gene_idx = list(adata.raw.var_names).index(gene)
        expression = adata.raw.X[:, gene_idx]
    else:
        gene_idx = list(adata.var_names).index(gene)
        expression = adata.X[:, gene_idx]

    # Convert to dense if sparse
    if hasattr(expression, "toarray"):
        expression = expression.toarray().flatten()
    else:
        expression = np.asarray(expression).flatten()

    # Calculate averages per cluster
    clusters = adata.obs[groupby].values
    unique_clusters = np.unique(clusters)

    results = []
    for cluster in unique_clusters:
        mask = clusters == cluster
        mean_expr = np.mean(expression[mask])
        cell_count = np.sum(mask)
        results.append({
            "cluster": str(cluster),
            "mean_expression": float(mean_expr),
            "cell_count": int(cell_count),
        })

    df = pd.DataFrame(results)
    df = df.sort_values("mean_expression", ascending=False)
    return df


def find_top_expressing_cluster(
    adata: AnnData,
    gene: str,
    groupby: str = "leiden",
) -> tuple[str, float]:
    """Find the cluster with highest average expression of a gene.

    Args:
        adata: The annotated data matrix.
        gene: Gene name to analyze.
        groupby: Observation key for grouping.

    Returns:
        Tuple of (cluster_id, mean_expression)
    """
    df = calculate_cluster_averages(adata, gene, groupby)
    top_row = df.iloc[0]
    return str(top_row["cluster"]), float(top_row["mean_expression"])


def get_cluster_cell_indices(
    adata: AnnData,
    cluster_id: str,
    groupby: str = "leiden",
) -> np.ndarray:
    """Get indices of cells belonging to a specific cluster.

    Args:
        adata: The annotated data matrix.
        cluster_id: The cluster identifier.
        groupby: Observation key for grouping.

    Returns:
        Array of cell indices.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"Grouping key '{groupby}' not found in observations.")

    mask = adata.obs[groupby].astype(str) == str(cluster_id)
    return np.where(mask)[0]


def rename_cluster_labels(
    adata: AnnData,
    old_name: str,
    new_name: str,
    groupby: str = "leiden",
) -> None:
    """Rename a cluster label in the AnnData object (in-place).

    Args:
        adata: The annotated data matrix.
        old_name: Current cluster name/ID.
        new_name: New name for the cluster.
        groupby: Observation key for grouping.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"Grouping key '{groupby}' not found in observations.")

    # Create a copy of the column to avoid SettingWithCopyWarning
    cluster_col = adata.obs[groupby].astype(str).copy()
    cluster_col[cluster_col == str(old_name)] = new_name
    adata.obs[groupby] = cluster_col

    logger.info("Renamed cluster '%s' to '%s' in '%s'", old_name, new_name, groupby)


def get_cluster_statistics(
    adata: AnnData,
    groupby: str = "leiden",
) -> pd.DataFrame:
    """Get summary statistics for all clusters.

    Args:
        adata: The annotated data matrix.
        groupby: Observation key for grouping.

    Returns:
        DataFrame with cluster statistics.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"Grouping key '{groupby}' not found in observations.")

    clusters = adata.obs[groupby].values
    unique_clusters = np.unique(clusters)

    results = []
    for cluster in unique_clusters:
        mask = clusters == cluster
        cell_count = np.sum(mask)
        percentage = (cell_count / len(clusters)) * 100

        results.append({
            "cluster": str(cluster),
            "cell_count": int(cell_count),
            "percentage": float(percentage),
        })

    df = pd.DataFrame(results)
    df = df.sort_values("cell_count", ascending=False)
    return df


def calculate_mito_percentage(adata: AnnData) -> None:
    """Calculate mitochondrial percentage if not already present.

    Identifies genes starting with "MT-" (case-insensitive) and calculates
    the percentage of counts from these genes for each cell.

    Args:
        adata: The annotated data matrix (modified in-place).

    Returns:
        None (modifies adata.obs in-place)
    """
    import scanpy as sc

    if "pct_counts_mt" in adata.obs.columns:
        logger.info("Mitochondrial percentage already calculated.")
        return

    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    n_mt_genes = adata.var["mt"].sum()

    if n_mt_genes == 0:
        logger.warning("No mitochondrial genes (MT-*) found in dataset.")
        adata.obs["pct_counts_mt"] = 0.0
        return

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True,
    )
    logger.info("Calculated mitochondrial percentage for %d cells (%d MT genes)", adata.n_obs, n_mt_genes)


def get_metadata_summary(
    adata: AnnData,
    max_unique_values: int = 50,
    exclude_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Scan adata.obs for categorical columns and return structured statistics.

    This function automatically inspects metadata columns in adata.obs and
    identifies categorical or category-like columns, returning their unique
    values and counts.

    Args:
        adata: The annotated data matrix.
        max_unique_values: Maximum number of unique values to consider a column
                          as categorical. Columns with more unique values are
                          excluded to avoid high-cardinality fields.
        exclude_columns: Optional list of column names to exclude from the summary.

    Returns:
        DataFrame with columns: column_name, dtype, n_unique, unique_values, value_counts
        Each row represents one categorical column from adata.obs.
    """
    if exclude_columns is None:
        exclude_columns = []

    results = []

    for col in adata.obs.columns:
        # Skip excluded columns
        if col in exclude_columns:
            continue

        # Get column data
        col_data = adata.obs[col]
        dtype = str(col_data.dtype)
        n_unique = col_data.nunique()

        # Only include categorical-like columns with reasonable cardinality
        # Skip high-cardinality columns (likely continuous or IDs)
        if n_unique > max_unique_values:
            logger.debug(f"Skipping column '{col}' with {n_unique} unique values (exceeds max_unique_values={max_unique_values})")
            continue

        # Skip numeric columns that look continuous (many unique values relative to total)
        if pd.api.types.is_numeric_dtype(col_data) and n_unique > min(20, len(col_data) * 0.1):
            logger.debug(f"Skipping numeric column '{col}' with {n_unique} unique values (likely continuous)")
            continue

        # Get unique values and counts
        value_counts = col_data.value_counts().to_dict()
        unique_values = list(value_counts.keys())

        # Convert to strings for consistent display
        unique_values_str = [str(v) for v in unique_values]
        value_counts_str = {str(k): v for k, v in value_counts.items()}

        results.append({
            "column_name": col,
            "dtype": dtype,
            "n_unique": n_unique,
            "unique_values": unique_values_str,
            "value_counts": value_counts_str,
        })

    # Create DataFrame
    if not results:
        logger.warning("No categorical columns found in adata.obs")
        return pd.DataFrame(columns=["column_name", "dtype", "n_unique", "unique_values", "value_counts"])

    df = pd.DataFrame(results)
    df = df.sort_values("n_unique")  # Sort by number of unique values

    logger.info(f"Found {len(df)} categorical columns in adata.obs")
    return df
