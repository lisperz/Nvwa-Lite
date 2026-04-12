"""Load `.h5ad` files with Seurat → Scanpy obsm key compatibility.

Files exported via Seurat's `SaveH5Seurat` + `Convert` pipeline store dimensional
reductions under uppercase keys (`UMAP`, `PCA`, `TSNE`, `DIFFMAP`). Scanpy and
every tool in this codebase expects the lowercase Scanpy convention
(`X_umap`, `X_pca`, `X_tsne`, `X_diffmap`). Renaming once at load time means the
30+ downstream tools and their key-presence checks work without modification.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
from anndata import AnnData

logger = logging.getLogger(__name__)

SEURAT_TO_SCANPY_OBSM: dict[str, str] = {
    "UMAP": "X_umap",
    "PCA": "X_pca",
    "TSNE": "X_tsne",
    "DIFFMAP": "X_diffmap",
}


def rename_seurat_obsm_keys(adata: AnnData) -> list[tuple[str, str]]:
    """Rename Seurat-style obsm keys in place to Scanpy convention.

    If the Scanpy-convention key already exists, the Seurat key is left alone
    (existing data wins — no clobber). Uses `obsm.pop` so the underlying matrix
    is not copied.

    Returns the list of `(old_key, new_key)` pairs that were renamed.
    """
    renames: list[tuple[str, str]] = []
    for seurat_key, scanpy_key in SEURAT_TO_SCANPY_OBSM.items():
        if seurat_key not in adata.obsm:
            continue
        if scanpy_key in adata.obsm:
            logger.warning(
                "Both '%s' and '%s' present in obsm; keeping existing '%s', "
                "leaving '%s' untouched",
                seurat_key, scanpy_key, scanpy_key, seurat_key,
            )
            continue
        adata.obsm[scanpy_key] = adata.obsm.pop(seurat_key)
        renames.append((seurat_key, scanpy_key))
        logger.info("Renamed obsm key '%s' -> '%s'", seurat_key, scanpy_key)
    return renames


def load_h5ad(path: str | Path) -> AnnData:
    """Read an `.h5ad` and apply Seurat → Scanpy obsm key renaming."""
    adata = ad.read_h5ad(path)
    rename_seurat_obsm_keys(adata)
    return adata
