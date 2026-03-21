"""Cluster and group resolution utilities.

Provides intelligent resolution of cluster identifiers, handling both
numeric cluster IDs and annotated cell type names.
"""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


def detect_grouping_columns(adata: AnnData) -> dict[str, list[str]]:
    """Detect available grouping columns in adata.obs.

    Returns:
        Dictionary with keys:
        - "cluster_columns": Columns that look like cluster IDs (leiden, louvain, etc.)
        - "annotation_columns": Columns that look like cell type annotations
    """
    cluster_columns = []
    annotation_columns = []

    for col in adata.obs.columns:
        col_lower = col.lower()

        # Check if it's a clustering column
        if any(x in col_lower for x in ["leiden", "louvain", "cluster", "seurat_cluster"]):
            cluster_columns.append(col)
        # Check if it's an annotation column
        elif any(x in col_lower for x in ["cell_type", "celltype", "annotation", "label"]):
            annotation_columns.append(col)

    return {
        "cluster_columns": cluster_columns,
        "annotation_columns": annotation_columns,
    }


def create_cluster_index_mapping(adata: AnnData, groupby: str) -> dict[int, str]:
    """Create a mapping from numeric indices to cluster/cell type names.

    When a grouping column contains cell type names instead of numeric IDs,
    this creates an ordered mapping so users can reference them by index.

    Args:
        adata: The annotated data matrix.
        groupby: The grouping column name.

    Returns:
        Dictionary mapping numeric index to cluster/cell type name.
        Sorted alphabetically for consistency.
    """
    unique_groups = sorted(adata.obs[groupby].astype(str).unique())
    return {idx: group for idx, group in enumerate(unique_groups)}


def parse_cluster_identifier(identifier: str) -> tuple[str, bool]:
    """Parse a cluster identifier and determine if it's numeric.

    Args:
        identifier: The identifier string (e.g., "Cluster 1", "cluster_5", "1", "B cells")

    Returns:
        Tuple of (cleaned_identifier, is_numeric)
        - cleaned_identifier: The numeric part if numeric, original string otherwise
        - is_numeric: True if the identifier represents a numeric cluster ID
    """
    # Remove common prefixes
    cleaned = identifier.strip()
    cleaned = re.sub(r"^cluster[_\s]*", "", cleaned, flags=re.IGNORECASE)

    # Check if it's purely numeric
    if cleaned.isdigit():
        return cleaned, True

    # Check if original was numeric
    if identifier.strip().isdigit():
        return identifier.strip(), True

    return identifier, False


def resolve_group_identifier(
    adata: AnnData,
    identifier: str,
    groupby: str | None = None,
) -> tuple[str, str]:
    """Resolve a group identifier to the actual value and column name.

    This function intelligently handles:
    - Numeric cluster IDs (0, 1, 2, "Cluster 1", "cluster_5")
    - Cell type annotations ("CD4 T cells", "B cells")
    - Auto-detection of the appropriate grouping column
    - Mapping numeric indices to cell type names when needed

    Args:
        adata: The annotated data matrix.
        identifier: The group identifier to resolve.
        groupby: Optional grouping column. If None, auto-detects.

    Returns:
        Tuple of (resolved_value, column_name)

    Raises:
        ValueError: If the identifier cannot be resolved.
    """
    # Parse the identifier
    cleaned_id, is_numeric = parse_cluster_identifier(identifier)

    # If groupby is specified, use it directly
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise ValueError(f"Grouping column '{groupby}' not found in adata.obs")

        available_groups = adata.obs[groupby].astype(str).unique().tolist()

        # Try exact match first
        if identifier in available_groups:
            return identifier, groupby

        # Try cleaned version
        if cleaned_id in available_groups:
            return cleaned_id, groupby

        # If numeric, try to find numeric match
        if is_numeric:
            for group in available_groups:
                if group == cleaned_id or group == identifier:
                    return group, groupby

            # Try numeric index mapping (e.g., "1" -> "B cells" if B cells is at index 1)
            try:
                numeric_idx = int(cleaned_id)
                index_mapping = create_cluster_index_mapping(adata, groupby)
                if numeric_idx in index_mapping:
                    resolved_name = index_mapping[numeric_idx]
                    logger.info(
                        f"Resolved numeric index {numeric_idx} to '{resolved_name}' "
                        f"in column '{groupby}'"
                    )
                    return resolved_name, groupby
            except (ValueError, KeyError):
                pass

        raise ValueError(
            f"Group '{identifier}' not found in '{groupby}'. "
            f"Available groups: {', '.join(available_groups)}"
        )

    # Auto-detect grouping column
    grouping_cols = detect_grouping_columns(adata)

    # If numeric identifier, prefer cluster columns
    if is_numeric:
        search_order = grouping_cols["cluster_columns"] + grouping_cols["annotation_columns"]
    else:
        # If text identifier, prefer annotation columns
        search_order = grouping_cols["annotation_columns"] + grouping_cols["cluster_columns"]

    # Try to find the identifier in available columns
    for col in search_order:
        available_groups = adata.obs[col].astype(str).unique().tolist()

        # Try exact match
        if identifier in available_groups:
            return identifier, col

        # Try cleaned version
        if cleaned_id in available_groups:
            return cleaned_id, col

        # For numeric identifiers, try to match numeric values
        if is_numeric:
            for group in available_groups:
                if group == cleaned_id:
                    return group, col

            # Try numeric index mapping
            try:
                numeric_idx = int(cleaned_id)
                index_mapping = create_cluster_index_mapping(adata, col)
                if numeric_idx in index_mapping:
                    resolved_name = index_mapping[numeric_idx]
                    logger.info(
                        f"Resolved numeric index {numeric_idx} to '{resolved_name}' "
                        f"in column '{col}'"
                    )
                    return resolved_name, col
            except (ValueError, KeyError):
                pass

    # If we get here, we couldn't resolve the identifier
    all_cols = grouping_cols["cluster_columns"] + grouping_cols["annotation_columns"]
    if not all_cols:
        raise ValueError(
            f"No grouping columns found in adata.obs. "
            f"Available columns: {', '.join(adata.obs.columns)}"
        )

    # Provide helpful error message with index mapping info
    error_msg = f"Could not resolve group identifier '{identifier}'. "
    error_msg += f"Searched in columns: {', '.join(all_cols)}. "

    # Show available groups from the first column
    first_col = all_cols[0]
    available = adata.obs[first_col].astype(str).unique().tolist()
    error_msg += f"Available groups in '{first_col}': {', '.join(available)}"

    # If numeric identifier was provided, show the index mapping
    if is_numeric:
        try:
            index_mapping = create_cluster_index_mapping(adata, first_col)
            error_msg += f". Numeric index mapping: {index_mapping}"
        except Exception:
            pass

    raise ValueError(error_msg)


def resolve_analysis_scope(
    adata: AnnData,
    target: str | None,
    groupby: str | None = None,
) -> tuple[str | None, str]:
    """Resolve analysis scope for DE/report generation.

    Determines whether the request is for:
    - All clusters (target is None or "all")
    - One specific cluster by numeric ID or name

    Args:
        adata: The annotated data matrix.
        target: Target cluster/group identifier, or None for all clusters.
        groupby: Optional grouping column. If None, auto-detects.

    Returns:
        Tuple of (resolved_target_or_none, column_name)
        - resolved_target_or_none: None for "all clusters", or the resolved group name
        - column_name: The grouping column to use

    Raises:
        ValueError: If the target cannot be resolved.
    """
    # Auto-detect grouping column if not specified
    if groupby is None:
        grouping_cols = detect_grouping_columns(adata)
        all_cols = grouping_cols["cluster_columns"] + grouping_cols["annotation_columns"]
        if not all_cols:
            raise ValueError(
                f"No grouping columns found in adata.obs. "
                f"Available columns: {', '.join(adata.obs.columns)}"
            )
        groupby = all_cols[0]

    # Check if target is "all" or None
    if target is None or target.lower() in ["all", "all clusters", "all groups"]:
        return None, groupby

    # Resolve specific target
    resolved_target, resolved_col = resolve_group_identifier(adata, target, groupby)
    return resolved_target, resolved_col


def resolve_pairwise_groups(
    adata: AnnData,
    group1: str,
    group2: str,
    groupby: str | None = None,
) -> tuple[str, str, str]:
    """Resolve two group identifiers for pairwise comparison.

    Ensures both groups are in the same grouping column.

    Args:
        adata: The annotated data matrix.
        group1: First group identifier.
        group2: Second group identifier.
        groupby: Optional grouping column. If None, auto-detects.

    Returns:
        Tuple of (resolved_group1, resolved_group2, column_name)

    Raises:
        ValueError: If groups cannot be resolved or are in different columns.
    """
    # Resolve first group
    resolved1, col1 = resolve_group_identifier(adata, group1, groupby)

    # Resolve second group using the same column
    try:
        resolved2, col2 = resolve_group_identifier(adata, group2, col1)
    except ValueError:
        # If second group not found in same column, provide helpful error
        available = adata.obs[col1].astype(str).unique().tolist()
        raise ValueError(
            f"Group '{group2}' not found in '{col1}' (where '{group1}' was found). "
            f"Available groups: {', '.join(available)}"
        )

    if col1 != col2:
        raise ValueError(
            f"Groups are in different columns: '{group1}' in '{col1}', "
            f"'{group2}' in '{col2}'. Both groups must be in the same column."
        )

    return resolved1, resolved2, col1
