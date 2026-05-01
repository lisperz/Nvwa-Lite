"""QC tools — Yalu Layer 1 section 2.

Tools:
  - qc_summary_table (Scenario 2B.1) — dataset-level summary table with flagging.
  - qc_violin_plot (Scenarios 2A.1 + 2A.2) — multi-panel QC metric violins
    grouped by cell type, optionally split by condition.
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
from src.domain.analysis.qc_metrics import resolve_qc_metric_column, summarize_qc_metrics

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


@register(
    description=(
        "Dataset-level QC summary table with automated flagging. Computes min, "
        "max, mean, median for QC metrics (nFeature_RNA / n_genes_by_counts, "
        "nCount_RNA / total_counts, pct_counts_mt) and flags elevated mitochondrial "
        "fraction (median > pct_mt_threshold) or low gene detection (median < "
        "n_feature_threshold). Use for prompts like 'give me a QC summary' or "
        "'show QC metrics statistics'. Do NOT use for per-group/per-cell-type QC."
    ),
    params={
        "pct_mt_threshold": {
            "description": "Median pct_counts_mt above this triggers a mito-fraction flag.",
        },
        "n_feature_threshold": {
            "description": "Median nFeature_RNA below this triggers a low-detection flag.",
        },
    },
)
def qc_summary_table(
    adata: "AnnData",
    pct_mt_threshold: float = 10.0,
    n_feature_threshold: int = 200,
) -> ArtifactResult:
    """Compute dataset-level QC summary + automated flagging per Yalu Layer 3.

    Wraps src.domain.analysis.qc_metrics.summarize_qc_metrics (cross-cutting
    primitive that handles column-name resolution across naming conventions:
    nFeature_RNA / n_genes_by_counts, nCount_RNA / total_counts, etc.).
    """
    try:
        df = summarize_qc_metrics(adata)
    except ValueError as e:
        raise ToolExecutionError(str(e), tool_name="qc_summary_table") from e

    # Resolve canonical names once for flag-branch matching
    pct_mt_col = resolve_qc_metric_column(adata, "pct_counts_mt")
    n_feature_col = resolve_qc_metric_column(adata, "n_genes_by_counts")

    # Per-row flag column per Yalu Layer 3 reference
    flags: list[str] = []
    for _, row in df.iterrows():
        metric = str(row["metric"])
        median = float(row["median"])
        flag = ""
        if metric == pct_mt_col and median > pct_mt_threshold:
            flag = (
                f"⚠️ Elevated mitochondrial fraction (median > {pct_mt_threshold}%). "
                f"May indicate cell stress."
            )
        elif metric == n_feature_col and median < n_feature_threshold:
            flag = (
                f"⚠️ Low gene detection (median < {n_feature_threshold}). "
                f"Consider reviewing filtering thresholds."
            )
        flags.append(flag)

    df = df.assign(Flag=flags)

    # Round numeric columns for display (CSV keeps full precision)
    df_display = df.copy()
    numeric_cols = [c for c in df_display.columns if c not in ("metric", "Flag")]
    df_display[numeric_cols] = df_display[numeric_cols].round(2)

    csv_data = df.to_csv(index=False)
    display_df = df_display.to_markdown(index=False)

    n_flagged = sum(1 for f in flags if f)
    if n_flagged:
        text = (
            f"QC summary computed for {adata.n_obs:,} cells across {len(df)} metric(s). "
            f"{n_flagged} metric(s) flagged for review."
        )
    else:
        text = (
            f"QC summary computed for {adata.n_obs:,} cells across {len(df)} metric(s). "
            f"All metrics within typical ranges."
        )

    return ArtifactResult(
        text=text,
        artifact_kind="csv",
        tool_name="qc_summary_table",
        params_used={
            "pct_mt_threshold": pct_mt_threshold,
            "n_feature_threshold": n_feature_threshold,
        },
        entities_acted_on=list(df["metric"]),
        csv_data=csv_data,
        display_df=display_df,
        code="qc_summary_table(adata)",
    )


@register(
    description=(
        "Multi-panel QC metric violin plot grouped by cell type (Yalu §2A). "
        "Renders nFeature_RNA / n_genes_by_counts, nCount_RNA / total_counts, "
        "and pct_counts_mt — auto-detected from adata.obs naming conventions. "
        "Optional split_by (condition column) renders a metric × condition "
        "grid (rows = metrics, columns = conditions) so QC ranges are "
        "comparable across conditions. For dataset-level QC numbers, use "
        "qc_summary_table. For single-gene expression violins, use violin_plot."
    ),
    params={
        "groupby": {
            "description": (
                "obs column for the x-axis grouping (typically the cell-type "
                "column). Empty string defaults to the cell-type column, "
                "falling back to the cluster column."
            ),
            "field_type": "obs_column",
        },
        "split_by": {
            "description": (
                "Optional obs column for per-condition columns in the grid "
                "(typically the condition column). Empty string = no split, "
                "single multi-panel figure (Yalu §2A.1)."
            ),
            "field_type": "obs_column",
        },
        "metrics": {
            "description": (
                "Optional list of QC metric column names. When omitted, "
                "auto-detects nFeature_RNA / nCount_RNA / pct_counts_mt "
                "across naming conventions."
            ),
        },
    },
)
def qc_violin_plot(
    adata: "AnnData",
    groupby: str = "",
    split_by: str = "",
    metrics: list[str] | None = None,
) -> ArtifactResult:
    """Render multi-panel QC violin plots; covers Yalu §2A.1 + §2A.2.

    §2A.1 (no split_by): single multi-panel figure with one violin per QC
    metric, grouped by cell type, mirroring Yalu Layer 3's
    ``sc.pl.violin(..., multi_panel=True)``.

    §2A.2 (split_by set): rows = metrics, columns = conditions. Each cell
    is a single-metric violin grouped by cell type on the condition's
    cell subset. Diverges from Yalu's literal output (N separate PNGs per
    condition) — composes into one matplotlib figure to fit the
    ArtifactResult single-image contract.
    """
    if not groupby:
        groupby = _default_groupby(adata)
    if groupby not in adata.obs.columns:
        raise ToolExecutionError(
            f"groupby '{groupby}' not in adata.obs.",
            tool_name="qc_violin_plot",
        )
    if split_by and split_by not in adata.obs.columns:
        raise ToolExecutionError(
            f"split_by '{split_by}' not in adata.obs.",
            tool_name="qc_violin_plot",
        )

    metric_cols = _resolve_metrics(adata, metrics)
    if not metric_cols:
        if metrics:
            # User asked for specific QC columns that don't exist in this
            # dataset. Surface what they asked for + what's available so the
            # responder can frame the missing-column case honestly without
            # falling back to "couldn't make the violin plot" framing.
            available = _resolve_metrics(adata, None)
            avail_str = (
                f" Available QC metrics in this dataset: {', '.join(available)}."
                if available else ""
            )
            raise ToolExecutionError(
                f"Requested QC column(s) {metrics!r} not found in adata.obs.{avail_str}",
                tool_name="qc_violin_plot",
            )
        raise ToolExecutionError(
            "No QC metrics found in adata.obs (looked for nFeature_RNA / "
            "n_genes_by_counts, nCount_RNA / total_counts, pct_counts_mt). "
            "Provide an explicit metrics list of obs column names.",
            tool_name="qc_violin_plot",
        )

    if split_by:
        conditions = [
            str(c) for c in adata.obs[split_by].dropna().unique()
        ]
        conditions.sort()
        image_bytes, code = _qc_violin_split(
            adata, groupby, split_by, metric_cols, conditions,
        )
    else:
        conditions = []
        image_bytes, code = _qc_violin_simple(adata, groupby, metric_cols)

    update_viz_state(
        "qc_violin",
        groupby=groupby,
        split_by=split_by or None,
    )

    parts = [
        f"QC violin for {len(metric_cols)} metric(s) "
        f"({', '.join(metric_cols)}) across {groupby}"
    ]
    if split_by:
        parts.append(f", split by {split_by} ({len(conditions)} conditions)")
    parts.append(f"; {adata.n_obs:,} cells.")

    return ArtifactResult(
        text="".join(parts),
        artifact_kind="image",
        tool_name="qc_violin_plot",
        params_used={
            "groupby": groupby,
            "split_by": split_by or None,
            "metrics": metric_cols,
        },
        entities_acted_on=list(metric_cols) + list(conditions),
        image_bytes=image_bytes,
        code=code,
    )


def _resolve_metrics(
    adata: "AnnData", metrics: list[str] | None,
) -> list[str]:
    """Return canonical adata.obs column names for QC metrics.

    When ``metrics`` is None: auto-detect via the same priority lists as
    ``summarize_qc_metrics`` (cells / genes / pct_counts_mt). When provided:
    canonicalize each name through ``resolve_qc_metric_column`` and drop
    unresolvable entries with a logged warning.
    """
    if metrics is None:
        cols: list[str] = []
        for cand in ("total_counts", "n_counts", "nCount_RNA"):
            if cand in adata.obs.columns:
                cols.append(cand)
                break
        for cand in ("n_genes_by_counts", "n_genes", "nFeature_RNA"):
            if cand in adata.obs.columns:
                cols.append(cand)
                break
        for cand in ("pct_counts_mt", "percent.mt", "percent_mito"):
            if cand in adata.obs.columns:
                cols.append(cand)
                break
        return cols

    resolved: list[str] = []
    for name in metrics:
        col = resolve_qc_metric_column(adata, name)
        if col is None:
            logger.warning("qc_violin_plot: metric %r not resolvable; skipping.", name)
            continue
        if col not in resolved:
            resolved.append(col)
    return resolved


def _qc_violin_simple(
    adata: "AnnData", groupby: str, metric_cols: list[str],
) -> tuple[bytes, str]:
    """Yalu §2A.1: single multi-panel figure with one violin per metric."""
    sc.pl.violin(
        adata,
        keys=metric_cols,
        groupby=groupby,
        rotation=45,
        show=False,
        multi_panel=True,
    )
    fig = plt.gcf()
    for ax in fig.get_axes():
        for label in ax.get_xticklabels():
            label.set_ha("right")
            label.set_rotation_mode("anchor")
    plt.tight_layout()

    code = (
        f"sc.pl.violin(adata, keys={metric_cols!r}, groupby={groupby!r}, "
        f"rotation=45, multi_panel=True, show=False)"
    )
    return _capture_png(), code


def _qc_violin_split(
    adata: "AnnData",
    groupby: str,
    split_by: str,
    metric_cols: list[str],
    conditions: list[str],
) -> tuple[bytes, str]:
    """Yalu §2A.2: metric × condition grid, single composed figure."""
    n_rows = len(metric_cols)
    n_cols = len(conditions)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(5 * n_cols, 3.5 * n_rows),
        squeeze=False,
    )

    for j, cond in enumerate(conditions):
        sub = adata[adata.obs[split_by].astype(str) == cond]
        for i, metric in enumerate(metric_cols):
            ax = axes[i][j]
            sc.pl.violin(
                sub, keys=[metric], groupby=groupby,
                rotation=45, ax=ax, show=False,
            )
            if i == 0:
                ax.set_title(f"{cond}")
            else:
                ax.set_title("")
            ax.set_xlabel("")
            ax.set_ylabel(metric if j == 0 else "")
            for label in ax.get_xticklabels():
                label.set_ha("right")
                label.set_rotation_mode("anchor")

    plt.tight_layout()

    code = (
        f"# Yalu §2A.2 grid: rows = {metric_cols!r}, cols = {conditions!r}\n"
        f"for cond in adata.obs[{split_by!r}].unique():\n"
        f"    sub = adata[adata.obs[{split_by!r}] == cond]\n"
        f"    sc.pl.violin(sub, keys={metric_cols!r}, "
        f"groupby={groupby!r}, multi_panel=True, show=False)"
    )
    return _capture_png(), code


def _default_groupby(adata: "AnnData") -> str:
    """Yalu §2A implicit default: cell-type column, then cluster column."""
    for cand in ("cell_type", "celltype", "annotation", "label"):
        if cand in adata.obs.columns:
            return cand
    for cand in ("leiden", "louvain", "seurat_clusters"):
        if cand in adata.obs.columns:
            return cand
    raise ToolExecutionError(
        "groupby is required and no cell-type or cluster column was "
        "detected in adata.obs. Please specify groupby explicitly.",
        tool_name="qc_violin_plot",
    )


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
