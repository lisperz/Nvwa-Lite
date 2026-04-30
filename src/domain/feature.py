"""Feature plot — Yalu Layer 1 §3A.

Single ``@register`` tool covering Yalu scenarios 3A.1–3A.4 via internal
branching on ``split_by`` + ``subset_key`` / ``subset_value``. Subset
variants (§3A.3, §3A.4) compose server-side per
``local/product/tool_migration_map.md``.

UMAP-colored-by-obs intent (Yalu §1) belongs to ``umap_plot``; this tool
colors the UMAP by gene expression only.
"""

from __future__ import annotations

import io
import logging
import math
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import scanpy as sc

from src.agent.viz_state import update_viz_state
from src.core.registry import register
from src.core.results import ArtifactResult, ToolExecutionError
from src.domain.subset import _apply_subset

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

_COLORMAP = "Reds"


@register(
    description=(
        "Feature plot — UMAP scatter colored by a single gene's expression "
        "(Reds colormap). Optional split_by renders per-category panels with "
        "shared expression range and shared UMAP axes. Optional subset_key + "
        "subset_value pre-filter to a single cell type before plotting "
        "(e.g. 'show MKI67 in cardiomyocytes only'). For coloring UMAP by an "
        "obs column (cell type / condition), use umap_plot."
    ),
    params={
        "gene": {
            "description": (
                "Single gene symbol to visualize (resolved against "
                "adata.var_names — case-insensitive, species-prefix aware)."
            ),
            "field_type": "gene",
        },
        "split_by": {
            "description": (
                "Optional obs column for per-category panels (e.g. condition "
                "column). Empty string = single panel."
            ),
            "field_type": "obs_column",
        },
        "subset_key": {
            "description": (
                "Optional obs column for cell-type subsetting before plotting. "
                "Pair with subset_value."
            ),
            "field_type": "obs_column",
        },
        "subset_value": {
            "description": (
                "Single cell-type value to keep (canonical name). §3A subset "
                "variants are single cell type only — broad terms that resolve "
                "to multiple subtypes raise a clarification turn upstream."
            ),
            "field_type": "cell_type",
        },
    },
)
def feature_plot(
    adata: "AnnData",
    gene: str,
    split_by: str = "",
    subset_key: str = "",
    subset_value: list[str] | None = None,
) -> ArtifactResult:
    """Render a feature plot (Yalu §3A.1–§3A.4).

    Param resolution: ``gene`` is canonicalized against ``adata.var_names``
    by ``gene_lookup`` upstream. ``split_by`` / ``subset_key`` are
    canonicalized by ``obs_column_lookup``. ``subset_value`` is canonicalized
    by ``cell_type_lookup`` — broad terms that match multiple subtypes raise
    a clarification turn before reaching this body.
    """
    sv = list(subset_value) if subset_value else []

    if bool(subset_key) != bool(sv):
        raise ToolExecutionError(
            "subset_key and subset_value must be provided together.",
            tool_name="feature_plot",
        )
    if sv and len(sv) != 1:
        raise ToolExecutionError(
            f"feature_plot subset is single cell type only; got {len(sv)}: {sv}.",
            tool_name="feature_plot",
        )

    if subset_key and sv:
        adata = _apply_subset(adata, subset_key, sv)

    if split_by:
        image_bytes, code = _feature_split(adata, gene, split_by)
    else:
        image_bytes, code = _feature_simple(adata, gene)

    update_viz_state(
        "feature",
        genes=[gene],
        split_by=split_by or None,
        subset_key=subset_key or None,
        subset_value=sv or None,
    )

    entities: list[str] = [gene]
    if split_by:
        entities.append(split_by)
    entities.extend(sv)

    parts = [f"Feature plot of {gene} on UMAP for {adata.n_obs:,} cells"]
    if split_by:
        parts.append(f", split by {split_by}")
    if sv:
        parts.append(f"; subsetted to {subset_key} ∈ {sv}")
    parts.append(".")

    return ArtifactResult(
        text="".join(parts),
        artifact_kind="image",
        tool_name="feature_plot",
        params_used={
            "gene": gene,
            "split_by": split_by or None,
            "subset_key": subset_key or None,
            "subset_value": sv or None,
        },
        entities_acted_on=entities,
        image_bytes=image_bytes,
        code=code,
    )


def _dot_size(n_obs: int) -> float:
    """Capped scanpy auto-formula. Floor=3 keeps dots visible at dpi=150
    on large datasets; cap=20 prevents auto-inflation on small subsets
    (matches scanpy default at ~6k cells)."""
    return max(3.0, min(20.0, 120000 / max(n_obs, 1)))


def _feature_simple(adata: "AnnData", gene: str) -> tuple[bytes, str]:
    size = _dot_size(adata.n_obs)
    sc.pl.umap(
        adata, color=[gene], color_map=_COLORMAP,
        size=size, frameon=False, show=False,
    )
    code = (
        f"sc.pl.umap(adata, color=[{gene!r}], color_map={_COLORMAP!r}, "
        f"frameon=False, show=False)"
    )
    return _capture_png(), code


def _feature_split(
    adata: "AnnData", gene: str, split_by: str
) -> tuple[bytes, str]:
    """Per-category panels with shared expression range + shared UMAP axes
    (Yalu Layer 3 §3A.2 spec). Per-panel colorbar via ``sc.pl.umap`` default
    — no shared bar.
    """
    conditions = sorted(adata.obs[split_by].dropna().unique().tolist())
    n = len(conditions)
    if n == 0:
        raise ToolExecutionError(
            f"Column '{split_by}' has no non-null values to split on.",
            tool_name="feature_plot",
        )

    if gene not in adata.var_names:
        raise ToolExecutionError(
            f"Gene '{gene}' not found in adata.var_names.",
            tool_name="feature_plot",
        )
    expr = adata[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    vmin, vmax = float(expr.min()), float(expr.max())

    plotted_mask = adata.obs[split_by].isin(conditions).to_numpy()
    coords = adata.obsm["X_umap"][plotted_mask]
    x_min, x_max = float(coords[:, 0].min()), float(coords[:, 0].max())
    y_min, y_max = float(coords[:, 1].min()), float(coords[:, 1].max())

    n_rows, n_cols = _grid_shape(n)
    size = _dot_size(adata.n_obs)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    axes_flat = _flatten_axes(axes, n_rows, n_cols)

    for idx, c in enumerate(conditions):
        ax = axes_flat[idx]
        sub = adata[adata.obs[split_by] == c]
        sc.pl.umap(
            sub, color=[gene], color_map=_COLORMAP,
            vmin=vmin, vmax=vmax, size=size,
            ax=ax, show=False, title=f"{gene} — {c}", frameon=False,
        )
        ax.set_xlim(x_min - 0.5, x_max + 0.5)
        ax.set_ylim(y_min - 0.5, y_max + 0.5)

    for idx in range(n, len(axes_flat)):
        axes_flat[idx].axis("off")

    plt.tight_layout()

    code = (
        f"# Per-category panels with shared expression range + UMAP axes\n"
        f"for c in adata.obs[{split_by!r}].dropna().unique():\n"
        f"    sc.pl.umap(adata[adata.obs[{split_by!r}] == c], "
        f"color=[{gene!r}], color_map={_COLORMAP!r}, vmin=vmin, vmax=vmax, "
        f"frameon=False, show=False)"
    )
    return _capture_png(), code


def _grid_shape(n: int) -> tuple[int, int]:
    """C4-locked layout: n≤3 → 1 row, 4-6 → 2 rows, ≥7 → 3 rows.
    Cols = ceil(n / n_rows). For n>9, rows stay at 3 and cols grow.
    """
    if n <= 3:
        return 1, n
    if n <= 6:
        return 2, math.ceil(n / 2)
    return 3, math.ceil(n / 3)


def _flatten_axes(axes, n_rows: int, n_cols: int) -> list:
    if n_rows == 1 and n_cols == 1:
        return [axes]
    if n_rows == 1 or n_cols == 1:
        return list(axes)
    return [ax for row in axes for ax in row]


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
