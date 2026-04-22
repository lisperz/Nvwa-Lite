#!/usr/bin/env python3
"""Manual verification script for Issue #25 - cluster numeric sort fix."""

import sys
from pathlib import Path

import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.analysis.cluster_resolution import create_cluster_index_mapping


def main():
    """Verify that integer cluster IDs sort numerically."""
    dataset_path = "data/uploads/GSE223414_test_subset.h5ad"

    print(f"Loading dataset: {dataset_path}")
    adata = sc.read_h5ad(dataset_path)

    print(f"Dataset shape: {adata.shape}")
    print(f"Available obs columns: {list(adata.obs.columns)}")

    # Check if seurat_clusters exists
    if "seurat_clusters" in adata.obs.columns:
        print("\nTesting with 'seurat_clusters' column:")
        mapping = create_cluster_index_mapping(adata, "seurat_clusters")

        print("\nCluster mapping (should be numeric order 0, 1, 2, ..., not 0, 1, 10, 11, ...):")
        for idx in sorted(mapping.keys()):
            print(f"  {idx}: {mapping[idx]}")

        # Verify numeric order
        cluster_ids = [mapping[i] for i in sorted(mapping.keys())]
        try:
            int_ids = [int(cid) for cid in cluster_ids]
            if int_ids == sorted(int_ids):
                print("\n✅ SUCCESS: Cluster IDs are in numeric order!")
            else:
                print("\n❌ FAIL: Cluster IDs are NOT in numeric order!")
                print(f"Expected: {sorted(int_ids)}")
                print(f"Got: {int_ids}")
        except ValueError:
            print("\n⚠️  Cluster IDs are not all integers (expected for some datasets)")
    else:
        print("\n⚠️  'seurat_clusters' column not found in this dataset")
        print("Available columns:", list(adata.obs.columns))


if __name__ == "__main__":
    main()
