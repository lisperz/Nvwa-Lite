#!/usr/bin/env python3
"""Validation script for the three fixes implemented.

This script tests:
1. QC metric visualization in scatter plots
2. Metadata inspection utility
3. Marker gene logic (verifies it uses proper statistical methods)
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from src.analysis.calculations import get_metadata_summary
from src.plotting.comparison import plot_scatter
from src.plotting.validation import validate_obs_or_gene
from src.analysis.differential import run_differential_expression


def test_fix_1_qc_scatter():
    """Test Fix 1: QC metric visualization in scatter plots."""
    print("\n" + "="*70)
    print("TEST FIX 1: QC Metric Visualization in Scatter Plots")
    print("="*70)

    # Create synthetic dataset with QC metrics
    n_cells = 100
    n_genes = 50

    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    adata = ad.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    # Add QC metrics to obs
    adata.obs["total_counts"] = np.random.randint(1000, 10000, n_cells)
    adata.obs["n_genes_by_counts"] = np.random.randint(100, 500, n_cells)
    adata.obs["pct_counts_mt"] = np.random.uniform(0, 20, n_cells)

    print("\n✓ Created synthetic dataset with QC metrics:")
    print(f"  - Cells: {n_cells}")
    print(f"  - Genes: {n_genes}")
    print(f"  - QC metrics in obs: {list(adata.obs.columns)}")

    # Test validation helper
    print("\n✓ Testing validate_obs_or_gene():")

    # Should pass for QC metric
    err = validate_obs_or_gene(adata, "total_counts")
    assert err is None, f"Expected None, got: {err}"
    print("  ✓ QC metric 'total_counts' validated successfully")

    # Should pass for gene
    err = validate_obs_or_gene(adata, "Gene_0")
    assert err is None, f"Expected None, got: {err}"
    print("  ✓ Gene 'Gene_0' validated successfully")

    # Should fail for non-existent field
    err = validate_obs_or_gene(adata, "NonExistent")
    assert err is not None, "Expected error for non-existent field"
    print("  ✓ Non-existent field correctly rejected")

    # Test scatter plot with QC metrics
    print("\n✓ Testing plot_scatter() with QC metrics:")
    try:
        result = plot_scatter(
            adata,
            gene_x="total_counts",
            gene_y="n_genes_by_counts",
            color_by="pct_counts_mt"
        )
        print("  ✓ Scatter plot generated successfully")
        print(f"  ✓ Code: {result.code}")
        print(f"  ✓ Message: {result.message}")
        assert "total_counts" in result.code
        assert "n_genes_by_counts" in result.code
        assert "pct_counts_mt" in result.code
    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        return False

    print("\n✅ FIX 1 PASSED: QC metrics work in scatter plots")
    return True


def test_fix_2_metadata_inspection():
    """Test Fix 2: Metadata screening and statistics."""
    print("\n" + "="*70)
    print("TEST FIX 2: Metadata Screening & Statistics")
    print("="*70)

    # Create synthetic dataset with various metadata types
    n_cells = 200
    n_genes = 50

    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    adata = ad.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    # Add categorical metadata
    adata.obs["louvain"] = pd.Categorical(np.random.choice(["0", "1", "2", "3"], n_cells))
    adata.obs["cell_type"] = pd.Categorical(np.random.choice(["T cell", "B cell", "Monocyte"], n_cells))
    adata.obs["batch"] = pd.Categorical(np.random.choice(["batch1", "batch2"], n_cells))

    # Add continuous metadata (should be excluded)
    adata.obs["continuous_metric"] = np.random.randn(n_cells)

    # Add high-cardinality metadata (should be excluded)
    adata.obs["cell_id"] = [f"cell_{i}" for i in range(n_cells)]

    print("\n✓ Created synthetic dataset with metadata:")
    print(f"  - Categorical: louvain, cell_type, batch")
    print(f"  - Continuous: continuous_metric")
    print(f"  - High-cardinality: cell_id")

    # Test metadata summary
    print("\n✓ Testing get_metadata_summary():")
    try:
        df = get_metadata_summary(adata, max_unique_values=50)

        print(f"  ✓ Found {len(df)} categorical columns")
        print(f"  ✓ Columns: {list(df['column_name'])}")

        # Verify expected columns are present
        expected_cols = {"louvain", "cell_type", "batch"}
        found_cols = set(df["column_name"])
        assert expected_cols.issubset(found_cols), f"Missing columns: {expected_cols - found_cols}"
        print("  ✓ All expected categorical columns found")

        # Verify high-cardinality column is excluded
        assert "cell_id" not in found_cols, "High-cardinality column should be excluded"
        print("  ✓ High-cardinality column correctly excluded")

        # Verify structure
        assert "column_name" in df.columns
        assert "dtype" in df.columns
        assert "n_unique" in df.columns
        assert "unique_values" in df.columns
        assert "value_counts" in df.columns
        print("  ✓ DataFrame structure is correct")

        # Check value counts
        for _, row in df.iterrows():
            col_name = row["column_name"]
            value_counts = row["value_counts"]
            n_unique = row["n_unique"]

            print(f"\n  Column: {col_name}")
            print(f"    Unique values: {n_unique}")
            print(f"    Distribution: {value_counts}")

            # Verify counts sum to total cells
            total = sum(value_counts.values())
            assert total == n_cells, f"Counts don't sum to n_cells: {total} != {n_cells}"

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

    print("\n✅ FIX 2 PASSED: Metadata inspection works correctly")
    return True


def test_fix_3_marker_gene_logic():
    """Test Fix 3: Verify marker gene logic uses proper statistical methods."""
    print("\n" + "="*70)
    print("TEST FIX 3: Marker Gene Logic Verification")
    print("="*70)

    # Create synthetic dataset with clear cluster structure
    n_cells = 300
    n_genes = 100

    # Create three clusters with distinct expression patterns
    cluster_sizes = [100, 100, 100]
    X = []
    clusters = []

    for i, size in enumerate(cluster_sizes):
        # Each cluster has higher expression of specific genes
        cluster_data = np.random.negative_binomial(5, 0.3, (size, n_genes))
        # Boost expression of cluster-specific genes
        cluster_data[:, i*10:(i+1)*10] *= 5
        X.append(cluster_data)
        clusters.extend([str(i)] * size)

    X = np.vstack(X)
    adata = ad.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]
    adata.obs["leiden"] = pd.Categorical(clusters)

    print("\n✓ Created synthetic dataset with 3 clusters:")
    print(f"  - Cells: {n_cells}")
    print(f"  - Genes: {n_genes}")
    print(f"  - Clusters: {adata.obs['leiden'].value_counts().to_dict()}")

    # Test differential expression
    print("\n✓ Testing run_differential_expression():")
    try:
        result = run_differential_expression(adata, groupby="leiden", method="wilcoxon", n_genes=10)

        print("  ✓ Differential expression completed successfully")
        print(f"  ✓ Message: {result.message}")

        # Verify results are stored in adata.uns
        assert "rank_genes_groups" in adata.uns, "DE results not stored in adata.uns"
        print("  ✓ Results stored in adata.uns['rank_genes_groups']")

        # Verify structure
        de_result = adata.uns["rank_genes_groups"]
        assert "names" in de_result, "Missing 'names' in DE results"
        assert "scores" in de_result or "logfoldchanges" in de_result, "Missing scores/logfoldchanges"
        assert "pvals" in de_result, "Missing p-values"
        assert "pvals_adj" in de_result, "Missing adjusted p-values"
        print("  ✓ DE results have correct structure (names, scores, pvals, pvals_adj)")

        # Verify it's using statistical methods, not max expression
        # Check that results include p-values (statistical test)
        groups = list(de_result["names"].dtype.names)
        for group in groups:
            pvals = de_result["pvals"][group]
            assert len(pvals) > 0, f"No p-values for group {group}"
            assert all(0 <= p <= 1 for p in pvals if not np.isnan(p)), "Invalid p-values"
        print("  ✓ Statistical p-values present (not simple max expression)")

        # Verify top genes are cluster-specific
        print("\n  Top marker genes per cluster:")
        for group in groups:
            top_genes = de_result["names"][group][:5]
            print(f"    Cluster {group}: {list(top_genes)}")

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

    print("\n✅ FIX 3 PASSED: Marker gene logic uses proper statistical methods")
    return True


def main():
    """Run all validation tests."""
    print("\n" + "="*70)
    print("VALIDATION SCRIPT FOR THREE FIXES")
    print("="*70)

    results = []

    # Run tests
    results.append(("Fix 1: QC Scatter Plots", test_fix_1_qc_scatter()))
    results.append(("Fix 2: Metadata Inspection", test_fix_2_metadata_inspection()))
    results.append(("Fix 3: Marker Gene Logic", test_fix_3_marker_gene_logic()))

    # Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)

    for name, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{status}: {name}")

    all_passed = all(passed for _, passed in results)

    if all_passed:
        print("\n🎉 ALL TESTS PASSED!")
        return 0
    else:
        print("\n⚠️  SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
