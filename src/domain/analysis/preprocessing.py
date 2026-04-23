"""Preprocessing pipeline for raw scRNA-seq data.

Runs the standard scanpy workflow:
QC → filter → normalize → HVG → PCA → UMAP → clustering.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


@dataclass
class PreprocessingResult:
    """Summary of a preprocessing run."""

    message: str
    steps_performed: list[str] = field(default_factory=list)
    n_cells_before: int = 0
    n_cells_after: int = 0
    n_genes_before: int = 0
    n_genes_after: int = 0


def run_preprocessing(
    adata: AnnData,
    *,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mito: float = 20.0,
    n_top_genes: int = 2000,
    resolution: float = 0.5,
) -> tuple[AnnData, PreprocessingResult]:
    """Run the full preprocessing pipeline on raw count data.

    Returns a new AnnData (subsetting creates copies) plus a summary.
    """
    result = PreprocessingResult(
        message="",
        n_cells_before=adata.n_obs,
        n_genes_before=adata.n_vars,
    )
    steps: list[str] = []
    adata = adata.copy()

    # QC metrics
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True,
    )
    steps.append("QC metrics calculated")
    logger.info("QC done: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    steps.append(f"Filtered (min_genes={min_genes}, min_cells={min_cells})")

    # Mito filter
    if "pct_counts_mt" in adata.obs.columns:
        before = adata.n_obs
        adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito, :].copy()
        removed = before - adata.n_obs
        steps.append(f"Removed {removed} high-mito cells (>{max_pct_mito}%)")
    logger.info("After filter: %d cells, %d genes", adata.n_obs, adata.n_vars)

    # Normalize + log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    steps.append("Normalized and log-transformed")

    # Save raw (log-normalized, full gene set) for DE analysis
    adata.raw = adata.copy()
    steps.append("Raw layer saved (log-normalized)")

    # Highly variable genes (seurat flavor works on log-normalized data)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var["highly_variable"]].copy()
    steps.append(f"Selected {adata.n_vars} highly variable genes")

    # Scale
    sc.pp.scale(adata, max_value=10)
    steps.append("Scaled")

    # PCA
    sc.tl.pca(adata, svd_solver="arpack")
    steps.append("PCA computed")

    # Neighbors + UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    steps.append("UMAP computed")

    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)
    steps.append(f"Leiden clustering (resolution={resolution})")
    logger.info("Done: %d cells, %d genes", adata.n_obs, adata.n_vars)

    result.steps_performed = steps
    result.n_cells_after = adata.n_obs
    result.n_genes_after = adata.n_vars
    result.message = (
        f"Preprocessing complete.\n"
        f"Cells: {result.n_cells_before:,} → {result.n_cells_after:,}\n"
        f"Genes: {result.n_genes_before:,} → {result.n_genes_after:,}\n"
        f"Steps: {'; '.join(steps)}"
    )
    return adata, result
