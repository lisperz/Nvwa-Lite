#!/usr/bin/env python3
"""Test script for DE intent disambiguation and scope handling."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import scanpy as sc
from src.analysis.differential import run_differential_expression
from src.analysis.cluster_resolution import resolve_analysis_scope


def test_intent_disambiguation():
    """Test DE intent disambiguation with pbmc3k dataset."""
    print("Loading pbmc3k_processed dataset...")
    adata = sc.datasets.pbmc3k_processed()

    print(f"\nDataset info:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Louvain clusters: {sorted(adata.obs['louvain'].unique())}")

    # Test 1: All clusters (Mode 1)
    print("\n" + "="*60)
    print("TEST 1: One-vs-rest for ALL clusters")
    print("="*60)
    result = run_differential_expression(adata, groupby='louvain', target_group=None)
    print(f"Result group: {result.group}")
    print(f"Result reference: {result.reference}")
    print(f"Message:\n{result.message}")

    # Test 2: Single cluster by cell type name (Mode 2)
    print("\n" + "="*60)
    print("TEST 2: One-vs-rest for SINGLE cluster (B cells)")
    print("="*60)
    result = run_differential_expression(adata, groupby='louvain', target_group="B cells")
    print(f"Result group: {result.group}")
    print(f"Result reference: {result.reference}")
    print(f"Message:\n{result.message}")
    print(f"Top 5 genes:\n{result.results_df.head()}")

    # Test 3: Single cluster by numeric ID (Mode 2)
    print("\n" + "="*60)
    print("TEST 3: One-vs-rest for SINGLE cluster (Cluster 3 = CD8 T cells)")
    print("="*60)
    # First resolve cluster 3
    resolved_target, resolved_col = resolve_analysis_scope(adata, "3", "louvain")
    print(f"Cluster 3 resolved to: {resolved_target}")

    result = run_differential_expression(adata, groupby='louvain', target_group=resolved_target)
    print(f"Result group: {result.group}")
    print(f"Result reference: {result.reference}")
    print(f"Message:\n{result.message}")
    print(f"Top 5 genes:\n{result.results_df.head()}")

    # Test 4: Resolve analysis scope
    print("\n" + "="*60)
    print("TEST 4: Resolve analysis scope")
    print("="*60)
    test_cases = [
        (None, "All clusters"),
        ("all", "All clusters"),
        ("B cells", "B cells"),
        ("3", "Cluster 3"),
        ("6", "Cluster 6"),
    ]
    for target, description in test_cases:
        try:
            resolved, col = resolve_analysis_scope(adata, target, "louvain")
            print(f"  {description:20s} -> resolved={resolved}, col={col}")
        except ValueError as e:
            print(f"  {description:20s} -> ERROR: {e}")

    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)


if __name__ == "__main__":
    test_intent_disambiguation()
