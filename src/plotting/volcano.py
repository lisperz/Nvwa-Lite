"""Volcano plot for differential expression results.

Pure matplotlib implementation with adjustText for non-overlapping labels.
"""

from __future__ import annotations

import io
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.plotting.executor import PlotResult

logger = logging.getLogger(__name__)


def _try_import_adjusttext():
    """Import adjustText if available, else return None."""
    try:
        from adjustText import adjust_text
        return adjust_text
    except ImportError:
        return None


def plot_volcano(
    de_df: pd.DataFrame,
    group: str,
    *,
    log2fc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    n_labels: int = 10,
) -> PlotResult:
    """Generate a volcano plot from DE results.

    Args:
        de_df: DataFrame with columns: gene, log2fc, pval, pval_adj.
        group: Group name (for title).
        log2fc_threshold: Fold-change cutoff for significance coloring.
        pval_threshold: Adjusted p-value cutoff.
        n_labels: Number of top genes to label.
    """
    df = de_df.copy()
    df["neg_log10_pval"] = -np.log10(df["pval_adj"].clip(lower=1e-300))

    sig_up = (df["log2fc"] > log2fc_threshold) & (df["pval_adj"] < pval_threshold)
    sig_down = (df["log2fc"] < -log2fc_threshold) & (df["pval_adj"] < pval_threshold)
    n_up = sig_up.sum()
    n_down = sig_down.sum()

    fig, ax = plt.subplots(figsize=(12, 8))

    ns = ~sig_up & ~sig_down
    ax.scatter(df.loc[ns, "log2fc"], df.loc[ns, "neg_log10_pval"],
               c="#BBBBBB", alpha=0.4, s=10, label=f"NS ({ns.sum():,})", zorder=1)
    ax.scatter(df.loc[sig_up, "log2fc"], df.loc[sig_up, "neg_log10_pval"],
               c="#e74c3c", alpha=0.7, s=14, label=f"Up ({n_up:,})", zorder=2)
    ax.scatter(df.loc[sig_down, "log2fc"], df.loc[sig_down, "neg_log10_pval"],
               c="#3498db", alpha=0.7, s=14, label=f"Down ({n_down:,})", zorder=2)

    ax.axhline(-np.log10(pval_threshold), ls="--", c="grey", lw=0.8, zorder=0)
    ax.axvline(log2fc_threshold, ls="--", c="grey", lw=0.8, zorder=0)
    ax.axvline(-log2fc_threshold, ls="--", c="grey", lw=0.8, zorder=0)

    # Legend outside the plot area to avoid collision with labels
    ax.legend(
        loc="upper left", bbox_to_anchor=(0, 1),
        frameon=True, framealpha=0.9, edgecolor="#CCCCCC",
        fontsize=9,
    )

    # Label top genes — truncate long names for readability
    df["score"] = df["neg_log10_pval"] * abs(df["log2fc"])
    top = df.nlargest(n_labels, "score")

    def _short_name(name: str, max_len: int = 18) -> str:
        return name if len(name) <= max_len else name[:max_len] + "…"

    adjust_text = _try_import_adjusttext()
    if adjust_text is not None:
        texts = []
        for _, row in top.iterrows():
            texts.append(ax.text(
                row["log2fc"], row["neg_log10_pval"],
                _short_name(row["gene"]),
                fontsize=7, fontweight="bold", alpha=0.85, zorder=3,
            ))
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.5),
            expand=(2.0, 2.0),
            force_text=(1.5, 1.5),
            force_points=(1.5, 1.5),
        )
    else:
        for _, row in top.iterrows():
            ax.annotate(
                _short_name(row["gene"]),
                (row["log2fc"], row["neg_log10_pval"]),
                fontsize=7, alpha=0.8,
                xytext=(5, 5), textcoords="offset points",
            )

    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax.set_title(f"Volcano Plot — {group}", fontsize=14, fontweight="bold")

    # Add some padding so labels at edges aren't clipped
    x_margin = (df["log2fc"].max() - df["log2fc"].min()) * 0.08
    ax.set_xlim(df["log2fc"].min() - x_margin, df["log2fc"].max() + x_margin)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    buf.seek(0)
    image = buf.read()

    code = (
        f'plot_volcano(de_df, group="{group}", '
        f"log2fc_threshold={log2fc_threshold}, pval_threshold={pval_threshold})"
    )
    message = f"Volcano plot for group {group}."
    return PlotResult(image=image, code=code, message=message)
