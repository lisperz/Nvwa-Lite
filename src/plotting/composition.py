"""Composition plots for cell distribution analysis.

Provides stacked bar charts to visualize cell composition across
conditions, samples, or other categorical variables.
"""

from __future__ import annotations

import io
import logging

import matplotlib.pyplot as plt
from anndata import AnnData

from src.analysis.composition import cross_tabulate_metadata
from src.plotting.executor import PlotResult

logger = logging.getLogger(__name__)


def plot_composition(
    adata: AnnData | None = None,
    row_key: str = "",
    col_key: str = "",
    kind: str = "count",
    crosstab: pd.DataFrame | None = None,
) -> PlotResult:
    """Generate a stacked bar chart showing cell composition.

    Args:
        adata: AnnData object with metadata in .obs (if crosstab not provided).
        row_key: Observation key for x-axis groups (e.g., 'condition').
        col_key: Observation key for stacked categories (e.g., 'cell_type').
        kind: Plot type - "count" for absolute counts, "percent" for percentages.
        crosstab: Pre-computed crosstab DataFrame. If provided, skips computation.

    Returns:
        PlotResult with the composition plot.

    Raises:
        ValueError: If keys are not found or kind is invalid.
    """
    if kind not in ("count", "percent"):
        raise ValueError(f"kind must be 'count' or 'percent', got '{kind}'")

    # Get cross-tabulation data
    if crosstab is None:
        if adata is None:
            raise ValueError("Must provide either adata or crosstab")
        normalize = kind == "percent"
        crosstab = cross_tabulate_metadata(adata, row_key, col_key, normalize=normalize)

    logger.info(
        "Plotting composition: %s × %s (kind=%s)",
        row_key, col_key, kind
    )

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 8))

    crosstab.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        colormap="tab20",
        edgecolor="white",
        linewidth=0.5,
    )

    # Formatting
    ax.set_xlabel(row_key, fontsize=12, fontweight="bold")
    if kind == "count":
        ax.set_ylabel("Cell Count", fontsize=12, fontweight="bold")
        title = f"Cell Composition: {col_key} across {row_key}"
    else:
        ax.set_ylabel("Percentage (%)", fontsize=12, fontweight="bold")
        title = f"Cell Composition (%): {col_key} across {row_key}"

    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    ax.legend(
        title=col_key,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        frameon=True,
        framealpha=0.9,
    )
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    # Rotate x-axis labels if needed
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")

    fig.tight_layout()

    # Save to buffer
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    buf.seek(0)
    image = buf.read()

    code = (
        f'plot_composition(adata, row_key="{row_key}", '
        f'col_key="{col_key}", kind="{kind}")'
    )
    message = f"Composition plot: {col_key} distribution across {row_key}."

    return PlotResult(image=image, code=code, message=message)
