"""QC tools — Yalu Layer 1 section 2.

Currently:
  - qc_summary_table (Scenario 2B.1)

Future (per local/product/tool_migration_map.md):
  - qc_violin_plot (Scenarios 2A.1 + 2A.2, collapsed with split_by)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from src.core.registry import register
from src.core.results import ArtifactResult, ToolExecutionError
from src.domain.analysis.qc_metrics import resolve_qc_metric_column, summarize_qc_metrics

if TYPE_CHECKING:
    from anndata import AnnData


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
