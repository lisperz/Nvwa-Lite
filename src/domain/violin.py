"""Violin plot — Yalu Layer 1 §3B.

Single ``@register`` tool covering Yalu scenarios 3B.1–3B.6 via internal
branching on ``groupby`` + ``split_by`` + ``subset_key`` / ``subset_value``.
Subset variants (§3B.4, §3B.5, §3B.6) compose server-side per
``local/product/tool_migration_map.md``.

Render path follows Yalu Layer 3 hybrid:
- non-split (3B.1, 3B.3, post-subset 3B.4, 3B.6) → ``sc.pl.violin`` at (8, 4)
- split (3B.2, post-subset 3B.5)                  → ``seaborn.violinplot``
  with hue=split_by at (15, 6)

QC violin (Yalu §2A) is a separate tool living elsewhere.
"""

from __future__ import annotations

import io
import logging
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from src.agent.viz_state import update_viz_state
from src.core.registry import register
from src.core.results import ArtifactResult, ToolExecutionError
from src.domain.subset import _apply_subset

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


@register(
    description=(
        "Violin plot — distribution of a single gene's expression across "
        "groups (cell type / condition / cluster). Optional split_by adds "
        "a second categorical via seaborn hue (Yalu §3B.2 / §3B.5). "
        "Optional subset_key + subset_value pre-filter before plotting "
        "(§3B.4 one cell type, §3B.5 multiple cell types, §3B.6 one "
        "condition). Single-gene only — for QC metric distributions use "
        "the qc violin tool."
    ),
    params={
        "gene": {
            "description": (
                "Single gene symbol to plot (resolved against adata.var_names "
                "— case-insensitive, species-prefix aware)."
            ),
            "field_type": "gene",
        },
        "groupby": {
            "description": (
                "obs column for the x-axis grouping. Empty string defaults "
                "to the cell-type column, falling back to the cluster column."
            ),
            "field_type": "obs_column",
        },
        "split_by": {
            "description": (
                "Optional obs column for split-violin via seaborn hue. Empty "
                "string = single per-group violin."
            ),
            "field_type": "obs_column",
        },
        "subset_key": {
            "description": (
                "Optional obs column for pre-filtering before plotting. Pair "
                "with subset_value. Can be a cell-type column (§3B.4 / §3B.5) "
                "or a condition column (§3B.6)."
            ),
            "field_type": "obs_column",
        },
        "subset_value": {
            "description": (
                "Values to keep on subset_key. Single value (§3B.4 / §3B.6) "
                "or multiple cell types (§3B.5). Resolver dispatches by "
                "subset_key's column role."
            ),
            "field_type": "cell_type",
        },
    },
)
def violin_plot(
    adata: "AnnData",
    gene: str,
    groupby: str = "",
    split_by: str = "",
    subset_key: str = "",
    subset_value: list[str] | None = None,
) -> ArtifactResult:
    """Render a violin plot (Yalu §3B.1–§3B.6).

    Param resolution: ``gene`` is canonicalized against ``adata.var_names``
    by ``gene_lookup``. ``groupby`` / ``split_by`` / ``subset_key`` are
    canonicalized by ``obs_column_lookup``. ``subset_value`` is dispatched
    by the resolver based on ``subset_key``'s column role (cell_type or
    condition lookup; falls back to a direct unique-value match for
    other roles).
    """
    sv = list(subset_value) if subset_value else []

    if bool(subset_key) != bool(sv):
        raise ToolExecutionError(
            "subset_key and subset_value must be provided together.",
            tool_name="violin_plot",
        )

    if not groupby:
        groupby = _default_groupby(adata)

    if subset_key and sv:
        adata = _apply_subset(adata, subset_key, sv, collision_dim=groupby)

    if split_by:
        image_bytes, code = _violin_split(adata, gene, groupby, split_by)
    else:
        image_bytes, code = _violin_simple(adata, gene, groupby)

    update_viz_state(
        "violin",
        genes=[gene],
        groupby=groupby,
        split_by=split_by or None,
        subset_key=subset_key or None,
        subset_value=sv or None,
    )

    entities: list[str] = [gene, groupby]
    if split_by:
        entities.append(split_by)
    entities.extend(sv)

    parts = [f"Violin plot of {gene} across {groupby} for {adata.n_obs:,} cells"]
    if split_by:
        parts.append(f", split by {split_by}")
    if sv:
        parts.append(f"; subsetted to {subset_key} ∈ {sv}")
    parts.append(".")

    return ArtifactResult(
        text="".join(parts),
        artifact_kind="image",
        tool_name="violin_plot",
        params_used={
            "gene": gene,
            "groupby": groupby,
            "split_by": split_by or None,
            "subset_key": subset_key or None,
            "subset_value": sv or None,
        },
        entities_acted_on=entities,
        image_bytes=image_bytes,
        code=code,
    )


def _violin_simple(
    adata: "AnnData", gene: str, groupby: str,
) -> tuple[bytes, str]:
    """sc.pl.violin path — Yalu §3B.1 / §3B.3 / post-subset §3B.4 / §3B.6."""
    if gene not in adata.var_names:
        raise ToolExecutionError(
            f"Gene '{gene}' not found in adata.var_names.",
            tool_name="violin_plot",
        )

    fig, ax = plt.subplots(figsize=(8, 4))
    sc.pl.violin(
        adata, keys=[gene], groupby=groupby,
        rotation=45, ax=ax, show=False,
    )
    ax.set_title(f"{gene} across {groupby}")
    ax.set_xlabel("")
    for label in ax.get_xticklabels():
        label.set_ha("right")
        label.set_rotation_mode("anchor")
    ax.set_ylim(bottom=0)
    plt.tight_layout()

    code = (
        f"sc.pl.violin(adata, keys=[{gene!r}], groupby={groupby!r}, "
        f"rotation=45, show=False)"
    )
    return _capture_png(), code


def _violin_split(
    adata: "AnnData", gene: str, groupby: str, split_by: str,
) -> tuple[bytes, str]:
    """seaborn split-violin path — Yalu §3B.2 / post-subset §3B.5."""
    if gene not in adata.var_names:
        raise ToolExecutionError(
            f"Gene '{gene}' not found in adata.var_names.",
            tool_name="violin_plot",
        )

    expr = adata[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr_flat = expr.flatten() if hasattr(expr, "flatten") else expr

    df = pd.DataFrame({
        "expression": expr_flat,
        groupby: adata.obs[groupby].values,
        split_by: adata.obs[split_by].values,
    })

    fig, ax = plt.subplots(figsize=(15, 6))
    sns.violinplot(
        data=df, x=groupby, y="expression", hue=split_by,
        ax=ax, inner="box", density_norm="width", split=False,
    )
    ax.set_ylim(bottom=0)
    ax.set_xlabel("")
    ax.set_title(f"{gene} across {groupby} by {split_by}")
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("right")
        label.set_rotation_mode("anchor")
    plt.tight_layout()

    code = (
        f"# seaborn split-violin: x={groupby}, hue={split_by}\n"
        f"df = pd.DataFrame({{'expression': adata[:, {gene!r}].X.toarray().flatten(), "
        f"{groupby!r}: adata.obs[{groupby!r}], {split_by!r}: adata.obs[{split_by!r}]}})\n"
        f"sns.violinplot(data=df, x={groupby!r}, y='expression', hue={split_by!r}, "
        f"inner='box', density_norm='width', split=False)"
    )
    return _capture_png(), code


def _default_groupby(adata: "AnnData") -> str:
    """Yalu §3B implicit default: cell-type column, then cluster column."""
    for cand in ("cell_type", "celltype", "annotation", "label"):
        if cand in adata.obs.columns:
            return cand
    for cand in ("leiden", "louvain", "seurat_clusters"):
        if cand in adata.obs.columns:
            return cand
    raise ToolExecutionError(
        "groupby is required and no cell-type or cluster column was "
        "detected in adata.obs. Please specify groupby explicitly.",
        tool_name="violin_plot",
    )


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
