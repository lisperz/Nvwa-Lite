"""Plot style helpers for spec-pipeline tools and global matplotlib/scanpy setup.

Lives at ``src/domain/`` root (not under ``src/domain/plotting/`` which is
legacy-path territory). Houses cross-cutting style primitives shared by
plot tools that need consistent visual conventions across the family — most
notably the Yalu condition color standard for split-by-condition modes —
plus the global matplotlib/scanpy configuration applied once at app startup.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc


# Yalu condition color standard (local/from_yalu/.../Nvwa_All_Features_Layer2.md)
_CONDITION_COLOR_PALETTE: list[str] = [
    "#2E86C1",  # Reference / Control — blue
    "#E74C3C",  # Treatment / Experimental — red
    "#27AE60",  # 3rd condition — green
    "#8E44AD",  # 4th condition — purple
]


def condition_color_map(
    conditions: list[str],
    reference: str | None = None,
) -> dict[str, str]:
    """Map each condition to a hex color per Yalu's standard palette.

    When ``reference`` is provided and present in ``conditions``, it gets the
    blue slot and the rest take red/green/purple in alphabetical order. When
    ``reference`` is None or not in conditions, the alphabetically-first
    condition gets blue. Beyond 4 conditions the palette wraps via modulo
    (rare; >4 conditions is an edge case in scRNA-seq comparisons).
    """
    if not conditions:
        return {}

    sorted_conds = sorted(str(c) for c in conditions)
    if reference and str(reference) in sorted_conds:
        ref_str = str(reference)
        sorted_conds = [ref_str] + [c for c in sorted_conds if c != ref_str]

    return {
        cond: _CONDITION_COLOR_PALETTE[i % len(_CONDITION_COLOR_PALETTE)]
        for i, cond in enumerate(sorted_conds)
    }


def configure_plot_style() -> None:
    """Configure scanpy and matplotlib for clean, legible plots."""
    matplotlib.use("Agg")  # Headless backend for Docker

    sc.set_figure_params(dpi=150, fontsize=12, frameon=False)
    sc.settings.verbosity = 0

    plt.rcParams.update({
        "figure.figsize": (10, 6),
        "font.family": "sans-serif",
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "figure.autolayout": True,
    })
