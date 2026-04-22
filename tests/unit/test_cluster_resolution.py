"""Unit tests for cluster resolution utilities."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from src.analysis.cluster_resolution import create_cluster_index_mapping


class TestCreateClusterIndexMapping:
    """Tests for create_cluster_index_mapping function."""

    def test_integer_cluster_ids_sort_numerically(self):
        """Integer cluster IDs should sort numerically, not lexicographically."""
        # Create test data with integer cluster IDs 0-20
        n_cells = 100
        cluster_ids = np.random.choice(range(21), size=n_cells).astype(str)
        adata = AnnData(
            X=np.random.rand(n_cells, 10),
            obs=pd.DataFrame({"cluster": cluster_ids}),
        )

        mapping = create_cluster_index_mapping(adata, "cluster")

        # Extract the cluster IDs in order
        ordered_clusters = [mapping[i] for i in sorted(mapping.keys())]

        # Should be numeric order: 0, 1, 2, ..., 20
        # NOT lexicographic: 0, 1, 10, 11, ..., 19, 2, 20, 3, ...
        expected = [str(i) for i in range(21)]
        assert ordered_clusters == expected, (
            f"Expected numeric sort {expected}, got {ordered_clusters}"
        )

    def test_string_labels_sort_alphabetically(self):
        """String cell type labels should sort alphabetically."""
        cell_types = ["T cells", "B cells", "NK cells", "Monocytes"]
        n_cells = 100
        cluster_ids = np.random.choice(cell_types, size=n_cells)
        adata = AnnData(
            X=np.random.rand(n_cells, 10),
            obs=pd.DataFrame({"celltype": cluster_ids}),
        )

        mapping = create_cluster_index_mapping(adata, "celltype")

        ordered_clusters = [mapping[i] for i in sorted(mapping.keys())]

        # Should be alphabetical
        expected = sorted(cell_types)
        assert ordered_clusters == expected

    def test_mixed_values_fall_back_to_string_sort(self):
        """Mixed integer and string values should fall back to string sort."""
        mixed_values = ["0", "1", "B cells", "10", "T cells"]
        n_cells = 100
        cluster_ids = np.random.choice(mixed_values, size=n_cells)
        adata = AnnData(
            X=np.random.rand(n_cells, 10),
            obs=pd.DataFrame({"cluster": cluster_ids}),
        )

        mapping = create_cluster_index_mapping(adata, "cluster")

        ordered_clusters = [mapping[i] for i in sorted(mapping.keys())]

        # Should fall back to lexicographic sort
        expected = sorted(mixed_values)
        assert ordered_clusters == expected
