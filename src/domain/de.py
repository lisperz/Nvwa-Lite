"""Differential expression — Yalu Layer 1 §4.2 + §5.

This module hosts the single-target-group DE tools. They share the
``_run_rank_genes`` private helper (same scanpy primitive
``sc.tl.rank_genes_groups`` with ``groups=[X]``) but write to DIFFERENT
uns slots so each tool's output stays schema-clean for its own consumer.

Tools:
- ``run_de`` (§5.1, §5.2, §5.3): pairwise group A vs group B.
- ``find_markers`` (§4.2): one-vs-rest DE for a single cell type.
- ``volcano_plot``: volcano render of the latest pairwise DE.

§5.4 (complex grouping with ``group_map``) deferred until the
``create_group_column`` substrate ships. ``find_all_markers`` (§4.1, multi-
group shape) lives in ``src/domain/markers.py`` — different scanpy call
shape (no ``groups`` arg), different uns slot, no shared helper.

State contract (three slots, content-named — single writer per slot):
- ``run_de`` writes to ``adata.uns[_DE_UNS_KEY]`` (``"nvwa_pairwise_de"``).
  ``volcano_plot`` reads from this slot. Pairwise output only.
- ``find_markers`` writes to ``adata.uns[_ONE_VS_REST_UNS_KEY]``
  (``"nvwa_one_vs_rest_de"``). No current reader; future multi-turn
  workflow (T-055) may add one.
- ``find_all_markers`` writes to ``adata.uns["nvwa_all_markers"]``
  (in markers.py). No current reader; future §4.3 dot_plot_top_markers.

``volcano_plot`` does NOT consume find_markers output — Yalu §5's volcano
follow-up is for pairwise DE only; §4 has no documented volcano scenario.
The structural slot separation IS the guard.

When run_de subsets (§5.3), the pairwise uns entry is bridged back to the
caller's adata so a follow-up ``volcano_plot`` turn finds it.

UX contract (Path D — two-turn explicit):
- ``run_de.text`` ends with an explicit volcano offer (Yalu §5 documents
  the post-DE volcano follow-up). The user sees the DE result first,
  decides whether to follow up. Symmetric to the statistical-disclaimer
  pattern (tool-body owns post-tool messaging, responder paraphrases).
- ``find_markers.text`` does NOT include a volcano offer — Yalu §4 frames
  marker discovery as characterization, not comparison; volcano is a §5
  follow-up.
- ``find_markers.text`` does NOT include the statistical disclaimer —
  Yalu §5 surfaces the Wilcoxon-on-single-cells caveat for DE comparison;
  §4's marker-discovery framing does not.

No imports from ``src/domain/plotting/`` or ``src/domain/analysis/`` —
those are legacy paths slated for M2 cleanup. Volcano render is a
self-contained matplotlib reimplementation per Yalu Layer 3.
"""

from __future__ import annotations

import io
import logging
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from src.agent.viz_state import update_viz_state
from src.core.registry import register
from src.core.results import ArtifactResult, ToolExecutionError
from src.domain.subset import _apply_subset

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


_DE_UNS_KEY = "nvwa_pairwise_de"
_ONE_VS_REST_UNS_KEY = "nvwa_one_vs_rest_de"
_STATISTICAL_DISCLAIMER = (
    "\n\nNote: Wilcoxon rank-sum test on single cells. For publication, "
    "consider pseudo-bulk DESeq2 — single-cell tests can produce overly "
    "significant p-values when multiple samples per condition exist."
)
_VOLCANO_OFFER = (
    "\n\nWould you like me to generate a volcano plot for these results?"
)


def _run_rank_genes(
    adata: "AnnData",
    groupby: str,
    groups: list[str],
    reference: str,
    method: str,
    key_added: str,
) -> pd.DataFrame:
    """Shared scanpy DE invocation for single-target-group calls.

    Used by run_de (pairwise: groups=[g1], reference=g2, key_added=
    _DE_UNS_KEY) and find_markers (one-vs-rest: groups=[celltype],
    reference="rest", key_added=_ONE_VS_REST_UNS_KEY). Same primitive,
    different reference value AND different uns slot — caller passes both.
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        method=method,
        pts=True,
        key_added=key_added,
    )
    return sc.get.rank_genes_groups_df(adata, group=groups[0], key=key_added)


def _bridge_uns_to_caller(
    source: "AnnData", target: "AnnData"
) -> None:
    """Copy the canonical DE uns entry from ``source`` to ``target`` when
    they are different objects (i.e., source is a subset copy). Lets a
    follow-up ``volcano_plot`` on the caller's adata find the DE results.
    """
    if source is not target:
        target.uns[_DE_UNS_KEY] = source.uns[_DE_UNS_KEY]


@register(
    description=(
        "Pairwise differential expression — find genes upregulated in group1 "
        "vs group2 within an obs column (groupby). Optional subset_key + "
        "subset_value pre-filter to one cell type before DE (Yalu §5.3). "
        f"Writes scanpy-native results to adata.uns[{_DE_UNS_KEY!r}] for "
        "follow-up volcano_plot. Always include the statistical disclaimer "
        "in the response."
    ),
    params={
        "groupby": {
            "description": (
                "obs column whose values are compared. Cell-type column for "
                "celltype-vs-celltype (§5.1); condition column for "
                "condition-vs-condition (§5.2 / §5.3 post-subset)."
            ),
            "field_type": "obs_column",
        },
        "group1": {
            "description": (
                "First group value in the groupby column (e.g. 'WT', "
                "'T cell'). Resolver dispatches against groupby's column role."
            ),
            "field_type": "condition",
        },
        "group2": {
            "description": (
                "Second group value in the groupby column (the reference for "
                "the comparison)."
            ),
            "field_type": "condition",
        },
        "subset_key": {
            "description": (
                "Optional obs column for pre-filtering (Yalu §5.3 single-"
                "celltype subset before DE). Pair with subset_value."
            ),
            "field_type": "obs_column",
        },
        "subset_value": {
            "description": (
                "Single cell-type value to keep (§5.3 is single-celltype "
                "subset only)."
            ),
            "field_type": "cell_type",
        },
        "method": {
            "description": (
                "Statistical method: 'wilcoxon' (default), 't-test', or "
                "'logreg'."
            ),
        },
        "n_top_genes": {
            "description": (
                "Top N genes per direction shown in the summary text. "
                "Default 10."
            ),
        },
        "pvals_adj_threshold": {
            "description": (
                "Adj p-value threshold for the 'significant' count and the "
                "exported CSV. Default 0.05."
            ),
        },
    },
)
def run_de(
    adata: "AnnData",
    groupby: str,
    group1: str,
    group2: str,
    subset_key: str = "",
    subset_value: list[str] | None = None,
    method: str = "wilcoxon",
    n_top_genes: int = 10,
    pvals_adj_threshold: float = 0.05,
) -> ArtifactResult:
    """Pairwise DE; covers Yalu §5.1, §5.2, §5.3 via internal subset."""
    sv = list(subset_value) if subset_value else []
    original_adata = adata

    if bool(subset_key) != bool(sv):
        raise ToolExecutionError(
            "subset_key and subset_value must be provided together.",
            tool_name="run_de",
        )
    if sv and len(sv) != 1:
        raise ToolExecutionError(
            f"run_de subset is single cell type only (Yalu §5.3); "
            f"got {len(sv)}: {sv}.",
            tool_name="run_de",
        )

    if subset_key and sv:
        adata = _apply_subset(
            original_adata, subset_key, sv, collision_dim=groupby,
        )

    if groupby not in adata.obs.columns:
        raise ToolExecutionError(
            f"groupby column '{groupby}' not in adata.obs.",
            tool_name="run_de",
        )

    unique_vals = {str(v) for v in adata.obs[groupby].unique()}
    if group1 not in unique_vals:
        raise ToolExecutionError(
            f"group1 '{group1}' not in adata.obs[{groupby!r}].unique() "
            f"({len(unique_vals)} values).",
            tool_name="run_de",
        )
    if group2 not in unique_vals:
        raise ToolExecutionError(
            f"group2 '{group2}' not in adata.obs[{groupby!r}].unique() "
            f"({len(unique_vals)} values).",
            tool_name="run_de",
        )

    de_df = _run_rank_genes(adata, groupby, [group1], group2, method, _DE_UNS_KEY)
    _bridge_uns_to_caller(adata, original_adata)

    sig_df = de_df[de_df["pvals_adj"] < pvals_adj_threshold]
    n_total = int(len(sig_df))
    n_up_g1 = int((sig_df["logfoldchanges"] > 0).sum())
    n_up_g2 = int((sig_df["logfoldchanges"] < 0).sum())

    top_g1 = sig_df[sig_df["logfoldchanges"] > 0].nlargest(n_top_genes, "logfoldchanges")
    top_g2 = sig_df[sig_df["logfoldchanges"] < 0].nsmallest(n_top_genes, "logfoldchanges")

    text_lines = [
        f"DE: {group1} vs {group2} on {groupby}, "
        f"{adata.n_obs:,} cells, method={method}.",
    ]
    if sv:
        text_lines.append(f"Subsetted to {subset_key} ∈ {sv}.")
    text_lines.append(
        f"{n_total} significant (padj<{pvals_adj_threshold}); "
        f"{n_up_g1} up in {group1}, {n_up_g2} up in {group2}."
    )
    text_lines.append("")
    text_lines.append(f"Top {min(n_top_genes, len(top_g1))} up in {group1}:")
    for _, row in top_g1.iterrows():
        text_lines.append(
            f"  {row['names']} (log2FC={row['logfoldchanges']:.2f}, "
            f"padj={row['pvals_adj']:.2g})"
        )
    text_lines.append("")
    text_lines.append(f"Top {min(n_top_genes, len(top_g2))} up in {group2}:")
    for _, row in top_g2.iterrows():
        text_lines.append(
            f"  {row['names']} (log2FC={row['logfoldchanges']:.2f}, "
            f"padj={row['pvals_adj']:.2g})"
        )

    text = "\n".join(text_lines) + _STATISTICAL_DISCLAIMER + _VOLCANO_OFFER

    update_viz_state(
        "de",
        groupby=groupby,
        group1=group1, group2=group2,
        subset_key=subset_key or None,
        subset_value=sv or None,
    )

    entities: list[str] = [groupby, group1, group2]
    entities.extend(sv)

    code = (
        f"sc.tl.rank_genes_groups(adata, groupby={groupby!r}, "
        f"groups=[{group1!r}], reference={group2!r}, method={method!r}, "
        f"key_added={_DE_UNS_KEY!r})\n"
        f"de_df = sc.get.rank_genes_groups_df(adata, group={group1!r}, "
        f"key={_DE_UNS_KEY!r})"
    )

    return ArtifactResult(
        text=text,
        artifact_kind="csv",
        tool_name="run_de",
        params_used={
            "groupby": groupby,
            "group1": group1,
            "group2": group2,
            "subset_key": subset_key or None,
            "subset_value": sv or None,
            "method": method,
            "n_top_genes": n_top_genes,
            "pvals_adj_threshold": pvals_adj_threshold,
        },
        entities_acted_on=entities,
        csv_data=sig_df.to_csv(index=False),
        display_df=sig_df.head(20).to_markdown(index=False),
        code=code,
    )


@register(
    description=(
        "One-vs-rest marker gene analysis for a SINGLE cell type "
        "(Yalu §4.2). Identifies genes specifically expressed in the "
        "target cell type compared to all other cells. Writes scanpy-"
        f"native results to adata.uns[{_ONE_VS_REST_UNS_KEY!r}]. Use "
        "run_de for pairwise comparison of two specific groups; use "
        "find_all_markers for survey across all cell types."
    ),
    params={
        "groupby": {
            "description": (
                "obs column containing cell type labels (typically the "
                "cell-type annotation column)."
            ),
            "field_type": "obs_column",
        },
        "celltype": {
            "description": (
                "Single cell-type value to find markers for "
                "(e.g. 'B cell'). Resolver expands broad terms; "
                "ambiguous expansions are rejected — pick one specific "
                "cell type."
            ),
            "field_type": "cell_type",
        },
        "method": {
            "description": (
                "Statistical method: 'wilcoxon' (default), 't-test', or "
                "'logreg'."
            ),
        },
        "pvals_adj_threshold": {
            "description": (
                "Adj p-value threshold for filtering. Default 0.05."
            ),
        },
        "logfc_threshold": {
            "description": (
                "log2FC threshold for filtering (only-positive markers "
                "when threshold>=0). Default 0 — Yalu §4.2 default returns "
                "positive markers."
            ),
        },
        "n_top_genes": {
            "description": (
                "Top N markers shown in the summary text. Default 10."
            ),
        },
    },
)
def find_markers(
    adata: "AnnData",
    groupby: str,
    celltype: list[str] | None = None,
    method: str = "wilcoxon",
    pvals_adj_threshold: float = 0.05,
    logfc_threshold: float = 0.0,
    n_top_genes: int = 10,
) -> ArtifactResult:
    """One-vs-rest marker analysis for a single cell type; covers Yalu §4.2."""
    ct = list(celltype) if celltype else []
    if len(ct) != 1:
        raise ToolExecutionError(
            f"find_markers requires a single cell type (Yalu §4.2); "
            f"got {len(ct)}: {ct}.",
            tool_name="find_markers",
        )
    target = ct[0]

    if groupby not in adata.obs.columns:
        raise ToolExecutionError(
            f"groupby column '{groupby}' not in adata.obs.",
            tool_name="find_markers",
        )

    unique_vals = {str(v) for v in adata.obs[groupby].unique()}
    if target not in unique_vals:
        raise ToolExecutionError(
            f"celltype '{target}' not in adata.obs[{groupby!r}].unique() "
            f"({len(unique_vals)} values).",
            tool_name="find_markers",
        )

    de_df = _run_rank_genes(
        adata, groupby, [target], "rest", method, _ONE_VS_REST_UNS_KEY,
    )

    filtered = de_df[
        (de_df["pvals_adj"] < pvals_adj_threshold)
        & (de_df["logfoldchanges"] > logfc_threshold)
    ]
    n_total = int(len(filtered))
    top = filtered.nlargest(n_top_genes, "logfoldchanges")

    text_lines = [
        f"Markers for {target} (one-vs-rest on {groupby}, "
        f"{adata.n_obs:,} cells, method={method}).",
        f"{n_total} significant (padj<{pvals_adj_threshold}, "
        f"logFC>{logfc_threshold}).",
        "",
        f"Top {min(n_top_genes, len(top))} markers:",
    ]
    for _, row in top.iterrows():
        text_lines.append(
            f"  {row['names']} (log2FC={row['logfoldchanges']:.2f}, "
            f"padj={row['pvals_adj']:.2g})"
        )

    text = "\n".join(text_lines)

    update_viz_state(
        "find_markers",
        groupby=groupby,
        celltype=target,
    )

    code = (
        f"sc.tl.rank_genes_groups(adata, groupby={groupby!r}, "
        f"groups=[{target!r}], reference='rest', method={method!r}, "
        f"key_added={_ONE_VS_REST_UNS_KEY!r})\n"
        f"de_df = sc.get.rank_genes_groups_df(adata, group={target!r}, "
        f"key={_ONE_VS_REST_UNS_KEY!r})"
    )

    return ArtifactResult(
        text=text,
        artifact_kind="csv",
        tool_name="find_markers",
        params_used={
            "groupby": groupby,
            "celltype": target,
            "method": method,
            "pvals_adj_threshold": pvals_adj_threshold,
            "logfc_threshold": logfc_threshold,
            "n_top_genes": n_top_genes,
        },
        entities_acted_on=[groupby, target],
        csv_data=filtered.to_csv(index=False),
        display_df=filtered.head(20).to_markdown(index=False),
        code=code,
    )


@register(
    description=(
        "Volcano plot of pairwise DE results. Hybrid resolution: when "
        "group1/group2/groupby ALL EMPTY, plots the most recent run_de "
        f"results from adata.uns[{_DE_UNS_KEY!r}]. When ALL THREE are set, "
        "verifies they match the stored DE; raises if they don't. Run "
        "run_de first — volcano_plot does not invoke DE itself."
    ),
    params={
        "group1": {
            "description": (
                "Optional first group; must match the prior run_de's group1. "
                "Empty plots the latest DE."
            ),
            "field_type": "condition",
        },
        "group2": {
            "description": (
                "Optional second group; must match the prior run_de's group2. "
                "Empty plots the latest DE."
            ),
            "field_type": "condition",
        },
        "groupby": {
            "description": (
                "Optional groupby column; must match the prior run_de's "
                "groupby. Empty plots the latest DE."
            ),
            "field_type": "obs_column",
        },
        "lfc_threshold": {
            "description": (
                "log2FC cutoff for significance coloring. Default 0.25."
            ),
        },
        "pvals_adj_threshold": {
            "description": (
                "adj-pvalue cutoff for significance coloring. Default 0.05."
            ),
        },
        "lfc_clip": {
            "description": "Clip log2FC to ±lfc_clip on the x-axis. Default 6.",
        },
        "n_label": {
            "description": "Number of top genes per direction to label. Default 10.",
        },
    },
)
def volcano_plot(
    adata: "AnnData",
    group1: str = "",
    group2: str = "",
    groupby: str = "",
    lfc_threshold: float = 0.25,
    pvals_adj_threshold: float = 0.05,
    lfc_clip: float = 6.0,
    n_label: int = 10,
) -> ArtifactResult:
    """Render a volcano plot for the latest pairwise DE."""
    if _DE_UNS_KEY not in adata.uns:
        raise ToolExecutionError(
            f"No DE results found in adata.uns[{_DE_UNS_KEY!r}]. Run run_de "
            "first.",
            tool_name="volcano_plot",
        )

    n_set = sum(bool(x) for x in (group1, group2, groupby))
    if n_set not in (0, 3):
        raise ToolExecutionError(
            "volcano_plot requires either ALL THREE of (group1, group2, "
            "groupby) set, or ALL THREE empty (plots latest run_de).",
            tool_name="volcano_plot",
        )

    stored = adata.uns[_DE_UNS_KEY]
    stored_params = stored.get("params", {})
    stored_groupby = str(stored_params.get("groupby", ""))
    stored_groups = list(stored_params.get("groups", []))
    stored_g1 = str(stored_groups[0]) if stored_groups else ""
    stored_g2 = str(stored_params.get("reference", ""))

    if n_set == 3:
        if (group1, group2, groupby) != (stored_g1, stored_g2, stored_groupby):
            raise ToolExecutionError(
                f"Requested (group1={group1!r}, group2={group2!r}, "
                f"groupby={groupby!r}) does not match the latest run_de "
                f"({stored_g1!r}, {stored_g2!r}, {stored_groupby!r}). "
                "Re-run run_de or omit these params to plot the latest.",
                tool_name="volcano_plot",
            )
        plot_g1, plot_g2, plot_gb = group1, group2, groupby
    else:
        plot_g1, plot_g2, plot_gb = stored_g1, stored_g2, stored_groupby

    de_df = sc.get.rank_genes_groups_df(adata, group=plot_g1, key=_DE_UNS_KEY)

    image_bytes, code = _render_volcano(
        de_df, plot_g1, plot_g2,
        lfc_threshold=lfc_threshold,
        pvals_adj_threshold=pvals_adj_threshold,
        lfc_clip=lfc_clip,
        n_label=n_label,
    )

    update_viz_state(
        "volcano",
        group1=plot_g1, group2=plot_g2, groupby=plot_gb,
    )

    return ArtifactResult(
        text=f"Volcano plot of {plot_g1} vs {plot_g2} on {plot_gb}.",
        artifact_kind="image",
        tool_name="volcano_plot",
        params_used={
            "group1": plot_g1,
            "group2": plot_g2,
            "groupby": plot_gb,
            "lfc_threshold": lfc_threshold,
            "pvals_adj_threshold": pvals_adj_threshold,
            "lfc_clip": lfc_clip,
            "n_label": n_label,
        },
        entities_acted_on=[plot_g1, plot_g2, plot_gb],
        image_bytes=image_bytes,
        code=code,
    )


def _render_volcano(
    de_df: pd.DataFrame,
    group1: str,
    group2: str,
    *,
    lfc_threshold: float,
    pvals_adj_threshold: float,
    lfc_clip: float,
    n_label: int,
) -> tuple[bytes, str]:
    """Pure-matplotlib volcano. Reimplemented per Yalu Layer 3 to remove
    the legacy ``src/domain/plotting/volcano.py`` dependency.
    """
    df = de_df.copy()
    df["neg_log10_pval"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))
    df["lfc_clipped"] = df["logfoldchanges"].clip(-lfc_clip, lfc_clip)

    sig_up = (df["pvals_adj"] < pvals_adj_threshold) & (
        df["logfoldchanges"] > lfc_threshold
    )
    sig_down = (df["pvals_adj"] < pvals_adj_threshold) & (
        df["logfoldchanges"] < -lfc_threshold
    )
    n_up = int(sig_up.sum())
    n_down = int(sig_down.sum())

    fig, ax = plt.subplots(figsize=(8, 6))

    ns = ~sig_up & ~sig_down
    ax.scatter(
        df.loc[ns, "lfc_clipped"], df.loc[ns, "neg_log10_pval"],
        c="lightgrey", alpha=0.5, s=10, zorder=1,
    )
    ax.scatter(
        df.loc[sig_up, "lfc_clipped"], df.loc[sig_up, "neg_log10_pval"],
        c="#E74C3C", alpha=0.7, s=14,
        label=f"Up in {group1} ({n_up:,})", zorder=2,
    )
    ax.scatter(
        df.loc[sig_down, "lfc_clipped"], df.loc[sig_down, "neg_log10_pval"],
        c="#2E86C1", alpha=0.7, s=14,
        label=f"Up in {group2} ({n_down:,})", zorder=2,
    )

    ax.axhline(
        -np.log10(pvals_adj_threshold),
        color="black", linestyle="--", linewidth=0.8, zorder=0,
    )
    ax.axvline(lfc_threshold, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    ax.axvline(-lfc_threshold, color="grey", linestyle="--", linewidth=0.8, zorder=0)

    top_up = df.loc[sig_up].nlargest(n_label, "neg_log10_pval")
    top_down = df.loc[sig_down].nlargest(n_label, "neg_log10_pval")
    top_genes = pd.concat([top_up, top_down])

    try:
        from adjustText import adjust_text
        texts = [
            ax.text(
                row["lfc_clipped"], row["neg_log10_pval"], row["names"],
                fontsize=7, fontweight="bold", zorder=3,
            )
            for _, row in top_genes.iterrows()
        ]
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.5),
        )
    except ImportError:
        for _, row in top_genes.iterrows():
            ax.annotate(
                row["names"],
                (row["lfc_clipped"], row["neg_log10_pval"]),
                fontsize=7, xytext=(5, 5), textcoords="offset points",
            )

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adj p-value)")
    ax.set_title(f"Volcano: {group1} vs {group2}")
    ax.legend(loc="upper left", fontsize=9)
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)

    code = (
        f"de_df = sc.get.rank_genes_groups_df(adata, group={group1!r}, "
        f"key={_DE_UNS_KEY!r})\n"
        f"# Volcano: {group1} vs {group2}, lfc_threshold={lfc_threshold}, "
        f"pvals_adj<{pvals_adj_threshold}, top {n_label} labels per direction"
    )
    return buf.read(), code
