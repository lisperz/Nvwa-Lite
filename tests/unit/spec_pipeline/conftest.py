"""Shared fixtures for spec-pipeline unit tests.

Loads `local/data/pbmc_test.h5ad` once per session, subsets to ~200 cells × 500
genes for speed, and returns a per-test copy so mutations don't leak across
tests. Augmentation (cell_type columns, uns["nvwa_meta"], condition columns) is
done inline in each test.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

try:
    import anndata as ad
except ImportError:
    ad = None

_REPO_ROOT = Path(__file__).resolve().parents[3]
_PBMC_PATH = _REPO_ROOT / "local" / "data" / "pbmc_test.h5ad"
_SUBSET_CELLS = 200
_SUBSET_GENES = 500
_RNG_SEED = 42


@pytest.fixture(scope="session")
def _pbmc_base():
    """Load + subset once per session. Never returned directly; use `adata`."""
    if ad is None:
        pytest.skip("anndata not installed")
    if not _PBMC_PATH.exists():
        pytest.skip(f"pbmc_test.h5ad not found at {_PBMC_PATH}")
    a = ad.read_h5ad(_PBMC_PATH)
    rng = np.random.default_rng(_RNG_SEED)
    n_cells = min(_SUBSET_CELLS, a.n_obs)
    n_genes = min(_SUBSET_GENES, a.n_vars)
    cell_idx = sorted(rng.choice(a.n_obs, size=n_cells, replace=False))
    gene_idx = sorted(rng.choice(a.n_vars, size=n_genes, replace=False))
    return a[cell_idx, gene_idx].copy()


@pytest.fixture
def adata(_pbmc_base):
    """Per-test copy of the subsetted PBMC AnnData. Mutate freely."""
    return _pbmc_base.copy()
