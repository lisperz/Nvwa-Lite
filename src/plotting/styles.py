"""Global matplotlib/scanpy style configuration for publication-ready figures."""

import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc


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
