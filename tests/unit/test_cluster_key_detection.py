"""Unit tests for _get_cluster_key() auto-detection logic."""

import pandas as pd
import pytest
from anndata import AnnData

from src.agent.tools import bind_dataset, _get_cluster_key


class TestClusterKeyDetection:
    """Test auto-detection of clustering columns in adata.obs."""

    def test_seurat_clusters_only(self):
        """When only seurat_clusters exists, it should be detected."""
        obs = pd.DataFrame({"seurat_clusters": ["0", "1", "0", "1"]})
        adata = AnnData(obs=obs)

        bind_dataset(adata)
        result = _get_cluster_key()

        assert result == "seurat_clusters"

    def test_leiden_priority(self):
        """When leiden exists, it should be detected first."""
        obs = pd.DataFrame({
            "leiden": ["0", "1", "0", "1"],
            "seurat_clusters": ["A", "B", "A", "B"]
        })
        adata = AnnData(obs=obs)

        bind_dataset(adata)
        result = _get_cluster_key()

        assert result == "leiden"

    def test_fallback_to_default(self):
        """When no clustering column exists, fallback to leiden."""
        obs = pd.DataFrame({"cell_type": ["T", "B", "T", "B"]})
        adata = AnnData(obs=obs)

        bind_dataset(adata)
        result = _get_cluster_key()

        assert result == "leiden"
