"""Cross-cutting subset primitive for plot tools.

Internal helper (NOT @register). Plot tools that absorb subset variants per
local/product/tool_migration_map.md call ``_apply_subset`` after the resolver
has canonicalized ``subset_value`` to filter adata before rendering.

The collision-dim guard is opt-in: tools whose primary visual axis collapses
to nothing useful when the subset reduces it to a single category (e.g.
violin's ``groupby``) pass that column as ``collision_dim``. Tools that
retain information at single-value (e.g. umap, where the point cloud itself
carries spatial meaning) omit it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from src.core.results import ToolExecutionError

if TYPE_CHECKING:
    from anndata import AnnData


class SubsetCollisionError(ToolExecutionError):
    """Subset collapses the caller-declared collision dimension to one value."""


def _apply_subset(
    adata: "AnnData",
    key: str,
    value: list[str],
    *,
    collision_dim: str | None = None,
) -> "AnnData":
    """Filter ``adata`` to rows where ``adata.obs[key]`` is in ``value``.

    Args:
        adata: Source dataset.
        key: obs column to filter on (already validated as present by caller
            or raised here).
        value: Canonical values to keep. Caller is responsible for canonicalizing
            via the resolver before calling — this helper does not do fuzzy
            matching.
        collision_dim: If provided, raise ``SubsetCollisionError`` when the
            resulting subset reduces ``adata.obs[collision_dim]`` to a single
            unique value.

    Returns:
        New AnnData (copy, not view) restricted to matching rows.
    """
    if key not in adata.obs.columns:
        raise ToolExecutionError(
            f"Subset column '{key}' not found in adata.obs.",
            tool_name="_apply_subset",
        )
    if not value:
        raise ToolExecutionError(
            f"Empty subset_value for column '{key}'.",
            tool_name="_apply_subset",
        )

    mask = adata.obs[key].isin(value)
    if int(mask.sum()) == 0:
        raise ToolExecutionError(
            f"No cells match {key} in {list(value)}.",
            tool_name="_apply_subset",
        )

    sub = adata[mask].copy()

    if collision_dim is not None and collision_dim in sub.obs.columns:
        if sub.obs[collision_dim].nunique() <= 1:
            raise SubsetCollisionError(
                f"Subset on '{key}' reduces '{collision_dim}' to a single "
                f"value; no comparison possible.",
                tool_name="_apply_subset",
            )

    return sub
