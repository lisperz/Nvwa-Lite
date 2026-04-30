"""Dot plot visualization — Yalu Layer 1 §3C.

Single ``@register`` tool covering Yalu §3C.1–§3C.4 via internal branching
on ``split_by`` + ``subset_key`` / ``subset_value``. Subset variants compose
server-side via ``_apply_subset``.

split_by mode: per Yalu §3C.2, creates a combined ``<groupby>_<split_by>``
groupby column on a defensive adata copy and color-codes the resulting
y-axis tick labels by condition (via ``condition_color_map``). The
``reference`` param anchors blue. (Dotplot's groupby labels live on the
y-axis; cf. heatmap which uses ``swap_axes=True`` and labels x-axis.)

The "top marker genes" alternative in Yalu §3C.1/§3C.2 prompts (a workflow
that pulls top-N from prior find_all_markers state) is **deferred to
T-055** — implementing it cleanly requires either a state-aware extractor
or workflow infra. Yalu evidence supports the design (5 scenarios across
§3C/§3D reference "top marker genes" as gene-list alternative); see
[T-055](pending_tasks.md). This module covers the concrete-gene-list
branch only.
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
from src.domain.plot_style import condition_color_map
from src.domain.subset import _apply_subset

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


@register(
    description=(
        "Dot plot showing expression of multiple genes across cell groups "
        "(cell type / cluster / condition). Dot size = % cells expressing; "
        "dot color = mean expression. Yalu §3C default for multi-gene "
        "overview prompts. Optional split_by combines groupby × condition "
        "into a single axis with condition-colored tick labels (Yalu §3C.2). "
        "Optional subset_key + subset_value pre-filter cells (Yalu §3C.3). "
        "The reference param anchors the control condition to blue in "
        "split mode. Use heatmap_plot only when the user explicitly asks "
        "for a heatmap or implies clustering / dendrogram view."
    ),
    params={
        "genes": {
            "description": "Gene names to display in the dot plot.",
            "field_type": "gene",
        },
        "groupby": {
            "description": (
                "obs column to group cells by (cell type / cluster / "
                "condition). Optional — defaults to the cell-type column, "
                "falling back to the cluster column."
            ),
            "field_type": "obs_column",
        },
        "split_by": {
            "description": (
                "Optional condition column. When set, creates a combined "
                "<groupby>_<split_by> axis with condition-colored tick labels."
            ),
            "field_type": "obs_column",
        },
        "subset_key": {
            "description": (
                "Optional obs column for cell-type subsetting before "
                "plotting. Pair with subset_value."
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
        "reference": {
            "description": (
                "Optional reference / control condition value. In split mode, "
                "anchors blue in the condition color palette; the rest take "
                "red / green / purple in alphabetical order. Empty = "
                "alphabetically-first condition gets blue."
            ),
            "field_type": "condition",
        },
    },
)
def dot_plot(
    adata: "AnnData",
    genes: list[str],
    groupby: str = "",
    split_by: str = "",
    subset_key: str = "",
    subset_value: list[str] | None = None,
    reference: str = "",
) -> ArtifactResult:
    """Render a dot plot; covers Yalu §3C.1–§3C.4 via internal composition."""
    sv = list(subset_value) if subset_value else []
    gene_list = list(genes) if genes else []

    if not gene_list:
        raise ToolExecutionError(
            "genes is required and must be a non-empty list.",
            tool_name="dot_plot",
        )

    if bool(subset_key) != bool(sv):
        raise ToolExecutionError(
            "subset_key and subset_value must be provided together.",
            tool_name="dot_plot",
        )

    if not groupby:
        groupby = _default_groupby(adata)

    if subset_key and sv:
        adata = _apply_subset(adata, subset_key, sv)

    if split_by:
        image_bytes, code = _dot_plot_split(
            adata, gene_list, groupby, split_by, reference or None
        )
    else:
        image_bytes, code = _dot_plot_simple(adata, gene_list, groupby)

    update_viz_state(
        "dotplot",
        genes=gene_list,
        groupby=groupby,
        split_by=split_by or None,
        subset_key=subset_key or None,
        subset_value=sv or None,
        reference=reference or None,
    )

    entities: list[str] = list(gene_list)
    entities.append(groupby)
    if split_by:
        entities.append(split_by)
    entities.extend(sv)

    if len(gene_list) <= 8:
        genes_str = ", ".join(gene_list)
    else:
        genes_str = ", ".join(gene_list[:8]) + f" … (+{len(gene_list) - 8} more)"
    parts = [f"Dot plot of {genes_str} across {groupby}"]
    if split_by:
        parts.append(f", split by {split_by}")
    if sv:
        parts.append(f"; subsetted to {subset_key} ∈ {sv}")
    if reference:
        parts.append(f"; reference condition: {reference}")
    parts.append(".")

    return ArtifactResult(
        text="".join(parts),
        artifact_kind="image",
        tool_name="dot_plot",
        params_used={
            "genes": gene_list,
            "groupby": groupby,
            "split_by": split_by or None,
            "subset_key": subset_key or None,
            "subset_value": sv or None,
            "reference": reference or None,
        },
        entities_acted_on=entities,
        image_bytes=image_bytes,
        code=code,
    )


def _dot_plot_simple(
    adata: "AnnData", genes: list[str], groupby: str
) -> tuple[bytes, str]:
    """Yalu §3C.1: single sc.pl.dotplot call. Scanpy auto-sizes."""
    sc.pl.dotplot(
        adata, var_names=genes, groupby=groupby, show=False,
    )
    code = (
        f"sc.pl.dotplot(adata, var_names={genes!r}, groupby={groupby!r})"
    )
    return _capture_png(), code


def _dot_plot_split(
    adata: "AnnData",
    genes: list[str],
    groupby: str,
    split_by: str,
    reference: str | None,
) -> tuple[bytes, str]:
    """Yalu §3C.2: combined groupby × split_by with condition-colored y-axis labels.

    Defensive copy of cells where both columns are non-null so we can add the
    combined column without mutating upstream state.
    """
    valid_mask = adata.obs[groupby].notna() & adata.obs[split_by].notna()
    if not valid_mask.any():
        raise ToolExecutionError(
            f"No cells have non-null values in both '{groupby}' and "
            f"'{split_by}'. Cannot build combined groupby axis.",
            tool_name="dot_plot",
        )
    sub = adata[valid_mask].copy()
    combined_col = f"{groupby}_{split_by}"
    sub.obs[combined_col] = (
        sub.obs[groupby].astype(str) + "_" + sub.obs[split_by].astype(str)
    )

    sc.pl.dotplot(
        sub, var_names=genes, groupby=combined_col, show=False,
    )

    # Color-code the groupby tick labels by condition. Dotplot puts groupby
    # on the y-axis (no swap_axes); locate the axis whose y-tick labels
    # match the combined-column format ("<group>_<condition>").
    fig = plt.gcf()
    conditions = sorted({str(c) for c in sub.obs[split_by].dropna().unique()})
    color_map = condition_color_map(conditions, reference=reference)

    for ax in fig.get_axes():
        labels = [t.get_text() for t in ax.get_yticklabels()]
        if not labels or not any("_" in lab for lab in labels):
            continue
        for label in ax.get_yticklabels():
            text = label.get_text()
            for cond, color in color_map.items():
                if text.endswith("_" + cond):
                    label.set_color(color)
                    break
            label.set_fontsize(6)
        break

    code = (
        f"# Combined groupby × condition axis with color-coded y-tick labels\n"
        f"adata.obs[{combined_col!r}] = "
        f"adata.obs[{groupby!r}].astype(str) + '_' + adata.obs[{split_by!r}].astype(str)\n"
        f"sc.pl.dotplot(adata, var_names={genes!r}, groupby={combined_col!r})"
    )
    return _capture_png(), code


def _default_groupby(adata: "AnnData") -> str:
    """Yalu §3C.1 implicit default: cell-type column, then cluster column."""
    for cand in ("cell_type", "celltype", "annotation", "label"):
        if cand in adata.obs.columns:
            return cand
    for cand in ("leiden", "louvain", "seurat_clusters"):
        if cand in adata.obs.columns:
            return cand
    raise ToolExecutionError(
        "groupby is required and no cell-type or cluster column was "
        "detected in adata.obs. Please specify groupby explicitly.",
        tool_name="dot_plot",
    )


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
