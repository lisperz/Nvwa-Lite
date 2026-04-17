"""Unit tests for src/analysis/h5ad_loader.py — Seurat → Scanpy obsm key rename.

Self-contained: builds tiny synthetic AnnData fixtures rather than relying on
the project dataset. Run with:
  PYTHONPATH=. uv run pytest tests/unit/test_h5ad_loader.py -v
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from src.analysis.h5ad_loader import (
    SEURAT_TO_SCANPY_OBSM,
    rename_seurat_obsm_keys,
)


def _make_adata(obsm: dict[str, np.ndarray] | None = None) -> AnnData:
    """Build a minimal AnnData with 10 cells × 5 genes and the given obsm dict."""
    X = np.zeros((10, 5), dtype=np.float32)
    a = AnnData(X=X)
    for key, arr in (obsm or {}).items():
        a.obsm[key] = arr
    return a


class TestRenameSeuratObsmKeys:

    def test_renames_all_four_seurat_keys(self):
        umap = np.random.rand(10, 2).astype(np.float32)
        pca = np.random.rand(10, 50).astype(np.float32)
        tsne = np.random.rand(10, 2).astype(np.float32)
        diffmap = np.random.rand(10, 5).astype(np.float32)
        a = _make_adata({"UMAP": umap, "PCA": pca, "TSNE": tsne, "DIFFMAP": diffmap})

        renames = rename_seurat_obsm_keys(a)

        assert set(a.obsm.keys()) == {"X_umap", "X_pca", "X_tsne", "X_diffmap"}
        assert set(renames) == set(SEURAT_TO_SCANPY_OBSM.items())
        np.testing.assert_array_equal(a.obsm["X_umap"], umap)
        np.testing.assert_array_equal(a.obsm["X_pca"], pca)

    def test_pure_scanpy_adata_is_noop(self):
        umap = np.random.rand(10, 2).astype(np.float32)
        a = _make_adata({"X_umap": umap})

        renames = rename_seurat_obsm_keys(a)

        assert renames == []
        assert list(a.obsm.keys()) == ["X_umap"]
        np.testing.assert_array_equal(a.obsm["X_umap"], umap)

    def test_conflict_keeps_existing_scanpy_key_and_does_not_clobber(self):
        scanpy_umap = np.ones((10, 2), dtype=np.float32)
        seurat_umap = np.zeros((10, 2), dtype=np.float32)
        a = _make_adata({"X_umap": scanpy_umap, "UMAP": seurat_umap})

        renames = rename_seurat_obsm_keys(a)

        assert ("UMAP", "X_umap") not in renames
        np.testing.assert_array_equal(a.obsm["X_umap"], scanpy_umap)
        assert "UMAP" in a.obsm  # left untouched
        np.testing.assert_array_equal(a.obsm["UMAP"], seurat_umap)

    def test_partial_seurat_keys(self):
        """Renames the Seurat keys that exist; leaves Scanpy/unrelated keys alone."""
        umap = np.random.rand(10, 2).astype(np.float32)
        existing_pca = np.random.rand(10, 50).astype(np.float32)
        a = _make_adata({"UMAP": umap, "X_pca": existing_pca})

        renames = rename_seurat_obsm_keys(a)

        assert renames == [("UMAP", "X_umap")]
        assert set(a.obsm.keys()) == {"X_umap", "X_pca"}
        np.testing.assert_array_equal(a.obsm["X_umap"], umap)
        np.testing.assert_array_equal(a.obsm["X_pca"], existing_pca)

    def test_empty_obsm_is_noop(self):
        a = _make_adata({})
        assert rename_seurat_obsm_keys(a) == []
        assert len(a.obsm) == 0

    def test_rename_is_reference_not_copy(self):
        """`pop` + assign should reuse the same array, not copy it."""
        umap = np.random.rand(10, 2).astype(np.float32)
        a = _make_adata({"UMAP": umap})
        original_id = id(a.obsm["UMAP"])

        rename_seurat_obsm_keys(a)

        # np.asarray is a no-op for ndarrays — same underlying buffer, no copy.
        assert id(a.obsm["X_umap"]) == original_id

    def test_dataframe_obsm_is_coerced_to_ndarray(self):
        """Seurat-converted .h5ad stores obsm as DataFrame — must become ndarray
        so scanpy's `basis_values[:, dims]` slicing works."""
        a = _make_adata()
        umap_df = pd.DataFrame(
            np.random.rand(10, 2).astype(np.float32),
            index=a.obs_names,
            columns=["umap_1", "umap_2"],
        )
        pca_df = pd.DataFrame(
            np.random.rand(10, 50).astype(np.float32),
            index=a.obs_names,
            columns=[f"PC_{i}" for i in range(1, 51)],
        )
        a.obsm["UMAP"] = umap_df
        a.obsm["PCA"] = pca_df

        rename_seurat_obsm_keys(a)

        assert isinstance(a.obsm["X_umap"], np.ndarray)
        assert isinstance(a.obsm["X_pca"], np.ndarray)
        # Scanpy-style slicing must work on the coerced value.
        sliced = a.obsm["X_umap"][:, (0, 1)]
        assert sliced.shape == (10, 2)
        np.testing.assert_array_equal(a.obsm["X_umap"], umap_df.to_numpy())
