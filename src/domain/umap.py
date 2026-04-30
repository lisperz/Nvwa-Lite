"""UMAP visualization — Yalu Layer 1 §1.

Single ``@register`` tool covering Yalu scenarios 1.1, 1.2, 1.3, 1.4 via
internal branching on ``split_by`` + ``subset_key`` / ``subset_value``.
Subset variants compose server-side rather than as a multi-step LLM
workflow, per ``local/product/tool_migration_map.md``.

Gene-on-UMAP intent (Yalu §3A) belongs to ``feature_plot``; ``color_by``
here is restricted to obs columns.
"""

from __future__ import annotations

import io
import logging
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


@register(
    description=(
        "UMAP scatter plot colored by an obs column (cell type / cluster / "
        "condition). Optional split_by renders per-category panels with "
        "unified axes. Optional subset_key + subset_value pre-filter cells "
        "before plotting (e.g. 'show UMAP for cardiomyocytes only'). For "
        "gene expression on UMAP, use feature_plot."
    ),
    params={
        "color_by": {
            "description": (
                "obs column to color cells by (cell type / cluster / condition). "
                "Optional — defaults to the cell-type column, falling back to "
                "the cluster column."
            ),
            "field_type": "obs_column",
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
                "Cell-type values to keep (canonical names or broad terms; "
                "broad terms like 'cardiomyocytes' resolve via cell_type_lookup "
                "and trigger a clarification turn if ambiguous)."
            ),
            "field_type": "cell_type",
        },
    },
)
def umap_plot(
    adata: "AnnData",
    color_by: str = "",
    split_by: str = "",
    subset_key: str = "",
    subset_value: list[str] | None = None,
) -> ArtifactResult:
    """Render a UMAP scatter; covers Yalu §1.1–§1.4 via internal composition.

    Param resolution: ``color_by`` / ``split_by`` / ``subset_key`` are
    canonicalized by ``obs_column_lookup`` upstream — values reaching this
    body are guaranteed to be present in ``adata.obs.columns``. ``color_by``
    falls back to the cell-type column (cell_type → celltype) and then the
    cluster column (leiden → louvain → seurat_clusters) per Yalu §1.1.
    """
    sv = list(subset_value) if subset_value else []

    if bool(subset_key) != bool(sv):
        raise ToolExecutionError(
            "subset_key and subset_value must be provided together.",
            tool_name="umap_plot",
        )

    if not color_by:
        color_by = _default_color_by(adata)

    if subset_key and sv:
        adata = _apply_subset(adata, subset_key, sv)

    if split_by:
        image_bytes, code = _umap_split(adata, color_by, split_by)
    else:
        image_bytes, code = _umap_simple(adata, color_by)

    update_viz_state(
        "umap",
        color_by=color_by,
        split_by=split_by or None,
        subset_key=subset_key or None,
        subset_value=sv or None,
    )

    entities: list[str] = [color_by]
    if split_by:
        entities.append(split_by)
    entities.extend(sv)

    parts = [f"UMAP rendered for {adata.n_obs:,} cells colored by {color_by}"]
    if split_by:
        parts.append(f", split by {split_by}")
    if sv:
        parts.append(f"; subsetted to {subset_key} ∈ {sv}")
    parts.append(".")

    return ArtifactResult(
        text="".join(parts),
        artifact_kind="image",
        tool_name="umap_plot",
        params_used={
            "color_by": color_by,
            "split_by": split_by or None,
            "subset_key": subset_key or None,
            "subset_value": sv or None,
        },
        entities_acted_on=entities,
        image_bytes=image_bytes,
        code=code,
    )


def _umap_simple(adata: "AnnData", color_by: str) -> tuple[bytes, str]:
    code = f"sc.pl.umap(adata, color=[{color_by!r}], frameon=False, show=False)"
    sc.pl.umap(adata, color=[color_by], frameon=False, show=False)
    return _capture_png(), code


def _umap_split(
    adata: "AnnData", color_by: str, split_by: str
) -> tuple[bytes, str]:
    """Per-category panels; xlim/ylim shared across panels = range over the
    union of cells that actually appear in any panel (Q3 lock).

    NaN values in ``split_by`` are excluded from panel iteration AND from the
    axis-range computation.
    """
    conditions = sorted(adata.obs[split_by].dropna().unique().tolist())
    n = len(conditions)
    if n == 0:
        raise ToolExecutionError(
            f"Column '{split_by}' has no non-null values to split on.",
            tool_name="umap_plot",
        )

    plotted_mask = adata.obs[split_by].isin(conditions).to_numpy()
    coords = adata.obsm["X_umap"][plotted_mask]
    x_min, x_max = float(coords[:, 0].min()), float(coords[:, 0].max())
    y_min, y_max = float(coords[:, 1].min()), float(coords[:, 1].max())

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    axes_list = [axes] if n == 1 else list(axes)

    for ax, c in zip(axes_list, conditions):
        sub = adata[adata.obs[split_by] == c]
        sc.pl.umap(
            sub, color=[color_by], ax=ax, show=False,
            title=str(c), frameon=False,
        )
        ax.set_xlim(x_min - 0.5, x_max + 0.5)
        ax.set_ylim(y_min - 0.5, y_max + 0.5)

    plt.tight_layout()
    code = (
        f"# Per-category panels with unified axes\n"
        f"for c in adata.obs[{split_by!r}].dropna().unique():\n"
        f"    sc.pl.umap(adata[adata.obs[{split_by!r}] == c], "
        f"color=[{color_by!r}], frameon=False, show=False)"
    )
    return _capture_png(), code


def _default_color_by(adata: "AnnData") -> str:
    """Yalu §1.1 implicit default: cell-type column, then cluster column."""
    for cand in ("cell_type", "celltype", "annotation", "label"):
        if cand in adata.obs.columns:
            return cand
    for cand in ("leiden", "louvain", "seurat_clusters"):
        if cand in adata.obs.columns:
            return cand
    raise ToolExecutionError(
        "color_by is required and no cell-type or cluster column was "
        "detected in adata.obs. Please specify color_by explicitly.",
        tool_name="umap_plot",
    )


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
