# Unit Testing Guide for Layer 1 Tools

## What Are Layer 1 Unit Tests?

**Layer 1** = Core analysis/plotting functions (the actual Scanpy/matplotlib code)

**Unit tests** = Test each individual function in isolation to verify it works correctly

**Why?** Catch bugs early, prevent regressions, ensure consistent behavior

---

## What Yuxin Wants You to Do

Complete unit tests for all Layer 1 functions in:
- `src/analysis/*.py` - Analysis functions (DE, clustering, QC, composition, etc.)
- `src/plotting/*.py` - Plotting functions (UMAP, violin, dot plot, etc.)

---

## Example: test_differential_simple.py

I've created a working example at `tests/unit/test_differential_simple.py` with 9 tests:

✅ **Basic functionality tests** - Does the function work?
✅ **Error handling tests** - Does it raise errors for invalid inputs?
✅ **Edge case tests** - Does it handle n_genes=0, empty groups, etc.?
✅ **State contamination test** - Does running twice give same results?

**All 9 tests pass!**

---

## How to Run Unit Tests

```bash
# Run all unit tests
cd nvwa-lite
PYTHONPATH=. uv run pytest tests/unit/ -v

# Run a specific test file
PYTHONPATH=. uv run pytest tests/unit/test_differential_simple.py -v

# Run a specific test function
PYTHONPATH=. uv run pytest tests/unit/test_differential_simple.py::test_run_differential_expression_all_groups -v
```

---

## What You Need to Create

### Priority 1: Core Analysis Functions

1. **`test_composition.py`** - Test `src/analysis/composition.py`
   - `cross_tabulate_metadata()` - cell count tables
   - Test: returns correct counts, handles missing groups, no state contamination

2. **`test_qc_metrics.py`** - Test `src/analysis/qc_metrics.py`
   - `summarize_qc_metrics()` - QC summary statistics
   - `get_obs_column_statistics()` - column stats
   - Test: correct calculations, handles missing columns

3. **`test_preprocessing.py`** - Test `src/analysis/preprocessing.py`
   - `run_preprocessing()` - normalization + clustering
   - Test: creates expected columns (X_pca, X_umap, leiden), doesn't crash

4. **`test_calculations.py`** - Test `src/analysis/calculations.py`
   - `calculate_mito_percentage()` - mito% calculation
   - `get_metadata_summary()` - metadata summaries
   - Test: correct math, handles edge cases

### Priority 2: Plotting Functions

5. **`test_plotting_executor.py`** - Test `src/plotting/executor.py`
   - `plot_umap()`, `plot_violin()`, `plot_feature()` - core plots
   - Test: returns PlotResult, PNG bytes are valid, no crashes

6. **`test_volcano.py`** - Test `src/plotting/volcano.py`
   - `plot_volcano()` - volcano plot generation
   - Test: returns valid PNG, handles empty DE results

---

## Test Template

Use this template for each new test file:

```python
"""Unit tests for [module name] (Layer 1).

Run with: PYTHONPATH=. pytest tests/unit/test_[module].py -v
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

from src.analysis.[module] import [function_to_test]


@pytest.fixture
def small_adata():
    """Create a small synthetic dataset for testing."""
    np.random.seed(42)
    n_obs = 100
    n_vars = 50

    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(float)
    adata = AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_vars)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_obs)]

    # Add metadata
    adata.obs["leiden"] = pd.Categorical(
        ["cluster_0"] * 40 + ["cluster_1"] * 35 + ["cluster_2"] * 25
    )

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata


def test_function_basic(small_adata):
    """Test basic functionality."""
    result = function_to_test(small_adata, param1="value1")

    assert result is not None
    assert isinstance(result, ExpectedType)
    # Add more assertions


def test_function_invalid_input(small_adata):
    """Test error handling."""
    with pytest.raises(ValueError):
        function_to_test(small_adata, param1="invalid")
```

---

## Tips

1. **Start simple** - Test that functions don't crash first
2. **Use small data** - 100 cells, 50 genes is enough
3. **Check types** - `isinstance(result, ExpectedType)`
4. **Check structure** - Does DataFrame have expected columns?
5. **Test errors** - Use `pytest.raises(ValueError)` for invalid inputs
6. **No state contamination** - Run function twice, compare results

---

## Known Bugs to Test For

From the CEO test (April 1), these bugs must be caught by unit tests:

1. **State contamination** (A4) - Same query returns different cell counts
   - Test: Run function twice with same params, assert results are identical

2. **Pairwise DE loop** (D2) - Cell type filter broken
   - Test: `run_pairwise_de()` with subset, verify it completes

3. **Dot plot top N bug** (C6) - Returns global genes instead of per-cell-type
   - Test: Request top 3 per cell type, verify each cell type has 3 genes

4. **Empty group handling** - Functions crash on empty groups
   - Test: Create adata with empty cluster, verify graceful error

---

## Success Criteria

✅ All test files created for Priority 1 modules
✅ At least 5 tests per module (basic + error + edge cases)
✅ All tests pass with `PYTHONPATH=. pytest tests/unit/ -v`
✅ Known bugs from CEO test are covered

---

## Questions?

- Look at `test_differential_simple.py` as a working example
- Look at `test_router.py` to see Yuxin's style
- Run tests frequently: `PYTHONPATH=. pytest tests/unit/test_[your_file].py -v`
