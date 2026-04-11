"""Cell composition analysis utilities.

Provides functions for cross-tabulating metadata columns to analyze
cell composition across conditions, samples, or other groupings.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


def cross_tabulate_metadata(
    adata: AnnData,
    row_key: str,
    col_key: str,
    normalize: bool = False,
) -> pd.DataFrame:
    """Cross-tabulate two metadata columns to show cell counts per group combination.

    This function performs a two-way groupby operation on adata.obs, equivalent to
    Seurat's: obj@meta.data %>% group_by(row_key, col_key) %>% summarise(n=n())

    Args:
        adata: The annotated data matrix.
        row_key: First metadata column (will be rows in output), e.g. 'orig.ident'.
        col_key: Second metadata column (will be columns in output), e.g. 'cell_type'.
        normalize: If True, return percentages instead of raw counts.
                  Percentages are calculated per row (e.g., percentage of each
                  cell type within each condition).

    Returns:
        DataFrame where:
        - Rows = unique values from row_key
        - Columns = unique values from col_key
        - Values = cell counts (or percentages if normalize=True)
        - Missing combinations are filled with 0

    Raises:
        ValueError: If either row_key or col_key is not found in adata.obs.

    Example:
        >>> ct = cross_tabulate_metadata(adata, "orig.ident", "cell_type")
        >>> print(ct)
                          B cells  T cells  NK cells
        Control-D5           1234     5678       901
        PA-IVS-1v-D5          567     2345       123
    """
    # Validate both keys exist
    if row_key not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(
            f"Row key '{row_key}' not found in metadata. "
            f"Available columns: {available}"
        )

    if col_key not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        raise ValueError(
            f"Column key '{col_key}' not found in metadata. "
            f"Available columns: {available}"
        )

    logger.info(
        "Cross-tabulating: %s (rows) × %s (columns), normalize=%s",
        row_key, col_key, normalize
    )

    # Perform cross-tabulation
    crosstab = adata.obs.groupby([row_key, col_key], observed=True).size().unstack(fill_value=0)

    if normalize:
        # Convert to percentages per row
        crosstab = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    return crosstab
