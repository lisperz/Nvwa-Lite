"""Heatmap visualization — Yalu Layer 1 §3D.

Single ``@register`` tool covering Yalu §3D.1–§3D.5 via internal branching
on ``split_by`` + ``subset_key`` / ``subset_value``. Subset variants compose
server-side via ``_apply_subset``.

split_by mode: per Yalu §3D.2, creates a combined ``<groupby>_<split_by>``
groupby column on a defensive adata copy and color-codes the resulting tick
labels by condition (via ``condition_color_map``). The ``reference`` param
anchors blue.
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
        "Heatmap showing expression of multiple genes across cell groups "
        "(cell type / cluster / condition), with a dendrogram for cluster "
        "similarity. Optional split_by combines groupby × condition into a "
        "single axis with condition-colored tick labels (Yalu §3D.2). "
        "Optional subset_key + subset_value pre-filter cells (e.g. 'heatmap "
        "for cardiomyocytes only'). The reference param anchors the control "
        "condition to blue in split mode."
    ),
    params={
        "genes": {
            "description": "Gene names to display in the heatmap.",
            "field_type": "gene",
        },
        "groupby": {
            "description": (
                "obs column to group cells by (cell type / cluster / condition). "
                "Optional — defaults to the cell-type column, falling back to "
                "the cluster column."
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
def heatmap_plot(
    adata: "AnnData",
    genes: list[str],
    groupby: str = "",
    split_by: str = "",
    subset_key: str = "",
    subset_value: list[str] | None = None,
    reference: str = "",
) -> ArtifactResult:
    """Render a heatmap; covers Yalu §3D.1–§3D.5 via internal composition."""
    sv = list(subset_value) if subset_value else []
    gene_list = list(genes) if genes else []

    if not gene_list:
        raise ToolExecutionError(
            "genes is required and must be a non-empty list.",
            tool_name="heatmap_plot",
        )

    if bool(subset_key) != bool(sv):
        raise ToolExecutionError(
            "subset_key and subset_value must be provided together.",
            tool_name="heatmap_plot",
        )

    if not groupby:
        groupby = _default_groupby(adata)

    if subset_key and sv:
        adata = _apply_subset(adata, subset_key, sv)

    if split_by:
        image_bytes, code = _heatmap_split(
            adata, gene_list, groupby, split_by, reference or None
        )
    else:
        image_bytes, code = _heatmap_simple(adata, gene_list, groupby)

    update_viz_state(
        "heatmap",
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
    parts = [f"Heatmap of {genes_str} across {groupby}"]
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
        tool_name="heatmap_plot",
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


def _heatmap_simple(
    adata: "AnnData", genes: list[str], groupby: str
) -> tuple[bytes, str]:
    """Yalu §3D.1: single sc.pl.heatmap with adaptive figsize + dendrogram."""
    n_groups = int(adata.obs[groupby].nunique())
    n_genes = len(genes)
    # When swap_axes=True: genes on x, groups on y. 0.3in/gene wide, 0.6in/group tall.
    fig_width = max(10, n_genes * 0.3)
    fig_height = max(6, n_groups * 0.6)

    sc.pl.heatmap(
        adata, var_names=genes, groupby=groupby,
        swap_axes=True, figsize=(fig_width, fig_height),
        dendrogram=True, cmap="viridis", show_gene_labels=True,
        show=False,
    )
    code = (
        f"sc.pl.heatmap(adata, var_names={genes!r}, groupby={groupby!r}, "
        f"swap_axes=True, dendrogram=True, cmap='viridis', "
        f"figsize=({fig_width:.1f}, {fig_height:.1f}))"
    )
    return _capture_png(), code


def _heatmap_split(
    adata: "AnnData",
    genes: list[str],
    groupby: str,
    split_by: str,
    reference: str | None,
) -> tuple[bytes, str]:
    """Yalu §3D.2: combined groupby × split_by with condition-colored tick labels.

    Defensive copy of cells where both columns are non-null so we can add the
    combined column without mutating upstream state.
    """
    valid_mask = adata.obs[groupby].notna() & adata.obs[split_by].notna()
    if not valid_mask.any():
        raise ToolExecutionError(
            f"No cells have non-null values in both '{groupby}' and "
            f"'{split_by}'. Cannot build combined groupby axis.",
            tool_name="heatmap_plot",
        )
    sub = adata[valid_mask].copy()
    combined_col = f"{groupby}_{split_by}"
    sub.obs[combined_col] = (
        sub.obs[groupby].astype(str) + "_" + sub.obs[split_by].astype(str)
    )

    n_groups = int(sub.obs[combined_col].nunique())
    n_genes = len(genes)
    fig_width = max(10, n_genes * 0.3)
    fig_height = max(6, n_groups * 0.6)

    sc.pl.heatmap(
        sub, var_names=genes, groupby=combined_col,
        swap_axes=True, figsize=(fig_width, fig_height),
        dendrogram=True, cmap="viridis", show_gene_labels=True,
        show=False,
    )

    # Color-code the groupby tick labels by condition. With swap_axes=True,
    # groupby labels live on the x-axis; locate the axis whose tick labels
    # match the combined-column format ("<group>_<condition>").
    fig = plt.gcf()
    conditions = sorted({str(c) for c in sub.obs[split_by].dropna().unique()})
    color_map = condition_color_map(conditions, reference=reference)

    for ax in fig.get_axes():
        labels = [t.get_text() for t in ax.get_xticklabels()]
        if not labels or not any("_" in lab for lab in labels):
            continue
        for label in ax.get_xticklabels():
            text = label.get_text()
            for cond, color in color_map.items():
                if text.endswith("_" + cond):
                    label.set_color(color)
                    break
            label.set_fontsize(6)
        break

    code = (
        f"# Combined groupby × condition axis with color-coded tick labels\n"
        f"adata.obs[{combined_col!r}] = "
        f"adata.obs[{groupby!r}].astype(str) + '_' + adata.obs[{split_by!r}].astype(str)\n"
        f"sc.pl.heatmap(adata, var_names={genes!r}, groupby={combined_col!r}, "
        f"swap_axes=True, dendrogram=True, cmap='viridis')"
    )
    return _capture_png(), code


def _default_groupby(adata: "AnnData") -> str:
    """Yalu §3D.1 implicit default: cell-type column, then cluster column."""
    for cand in ("cell_type", "celltype", "annotation", "label"):
        if cand in adata.obs.columns:
            return cand
    for cand in ("leiden", "louvain", "seurat_clusters"):
        if cand in adata.obs.columns:
            return cand
    raise ToolExecutionError(
        "groupby is required and no cell-type or cluster column was "
        "detected in adata.obs. Please specify groupby explicitly.",
        tool_name="heatmap_plot",
    )


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
