#!/usr/bin/env python3
"""Test script for cluster resolution with numeric indices."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import scanpy as sc
from src.analysis.cluster_resolution import (
    create_cluster_index_mapping,
    parse_cluster_identifier,
    resolve_group_identifier,
    resolve_pairwise_groups,
)


def test_cluster_resolution():
    """Test cluster resolution with pbmc3k dataset."""
    print("Loading pbmc3k_processed dataset...")
    adata = sc.datasets.pbmc3k_processed()

    print(f"\nDataset info:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Obs columns: {list(adata.obs.columns)}")

    # Check louvain column
    print(f"\nLouvain column unique values:")
    louvain_values = sorted(adata.obs['louvain'].unique())
    print(f"  {louvain_values}")

    # Test create_cluster_index_mapping
    print("\n" + "="*60)
    print("TEST 1: Create cluster index mapping")
    print("="*60)
    mapping = create_cluster_index_mapping(adata, 'louvain')
    print("Numeric index -> Cell type mapping:")
    for idx, name in mapping.items():
        print(f"  {idx}: {name}")

    # Test parse_cluster_identifier
    print("\n" + "="*60)
    print("TEST 2: Parse cluster identifiers")
    print("="*60)
    test_cases = ["1", "Cluster 1", "cluster_5", "B cells", "CD4 T cells"]
    for identifier in test_cases:
        cleaned, is_numeric = parse_cluster_identifier(identifier)
        print(f"  '{identifier}' -> cleaned='{cleaned}', is_numeric={is_numeric}")

    # Test resolve_group_identifier with numeric IDs
    print("\n" + "="*60)
    print("TEST 3: Resolve numeric cluster IDs")
    print("="*60)
    numeric_tests = ["0", "1", "5", "Cluster 1", "cluster_3"]
    for identifier in numeric_tests:
        try:
            resolved, col = resolve_group_identifier(adata, identifier, 'louvain')
            print(f"  '{identifier}' -> '{resolved}' (column: {col})")
        except ValueError as e:
            print(f"  '{identifier}' -> ERROR: {e}")

    # Test resolve_group_identifier with cell type names
    print("\n" + "="*60)
    print("TEST 4: Resolve cell type names")
    print("="*60)
    celltype_tests = ["B cells", "CD4 T cells", "NK cells"]
    for identifier in celltype_tests:
        try:
            resolved, col = resolve_group_identifier(adata, identifier, 'louvain')
            print(f"  '{identifier}' -> '{resolved}' (column: {col})")
        except ValueError as e:
            print(f"  '{identifier}' -> ERROR: {e}")

    # Test resolve_pairwise_groups
    print("\n" + "="*60)
    print("TEST 5: Resolve pairwise comparisons")
    print("="*60)
    pairwise_tests = [
        ("1", "5"),
        ("Cluster 1", "Cluster 5"),
        ("B cells", "CD4 T cells"),
        ("1", "B cells"),  # Mixed: numeric and cell type
    ]
    for group1, group2 in pairwise_tests:
        try:
            resolved1, resolved2, col = resolve_pairwise_groups(
                adata, group1, group2, 'louvain'
            )
            print(f"  '{group1}' vs '{group2}' -> '{resolved1}' vs '{resolved2}' (column: {col})")
        except ValueError as e:
            print(f"  '{group1}' vs '{group2}' -> ERROR: {e}")

    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)


if __name__ == "__main__":
    test_cluster_resolution()
