"""Shared test fixtures and dataset profile discovery for Nvwa-Lite test suite.

Environment variables:
  DATASET_PATH     — path to any .h5ad file (required by Layer 1 tests)
  DATASET_NAME     — "temple" enables Temple-specific regression tests
  SAVE_TEST_PLOTS  — "1" writes plot PNGs to test-records/ for visual review
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from anndata import AnnData

# ---------------------------------------------------------------------------
# Repo root on sys.path so all tests can `from src.xxx import yyy`
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# ---------------------------------------------------------------------------
# Environment config
# ---------------------------------------------------------------------------
PNG_MAGIC = b"\x89PNG"
MIN_PLOT_BYTES = 5_000

_SAVE_PLOTS: bool = os.environ.get("SAVE_TEST_PLOTS", "0") == "1"

# Default to Temple path for backwards compatibility
DATASET_PATH = Path(
    os.environ.get("DATASET_PATH")
    or os.environ.get("TEMPLE_DATASET_PATH", "/Users/zhuchen/Downloads/GSE223414_slim.h5ad")
)
DATASET_NAME: str = os.environ.get("DATASET_NAME", "").lower()

# Per-dataset plot subdirectory — prevents files from different datasets overwriting each other
# e.g. test-records/2026-04-10_layer1-plots/temple/
_PLOT_DIR = (
    _REPO_ROOT
    / "test-records"
    / f"{datetime.now().strftime('%Y-%m-%d')}_layer1-plots"
    / (DATASET_NAME or DATASET_PATH.stem)
)


# ---------------------------------------------------------------------------
# DatasetProfile — describes any dataset's structure without biology assumptions
# ---------------------------------------------------------------------------

@dataclass
class DatasetProfile:
    """Portable description of a loaded AnnData dataset's structure.

    Discovered automatically from adata by discover_profile().
    Tests use this instead of hardcoding column names or biology.
    """
    name: str                            # "temple" | "pbmc" | "unknown"
    n_cells: int
    n_genes: int

    # Metadata columns (None if not found)
    cell_type_column: str | None         # e.g. "cell_type", "louvain"
    condition_column: str | None         # e.g. "orig.ident" — None for single-sample datasets
    cluster_column: str | None           # e.g. "seurat_clusters", "louvain"

    # Grouping column for DE/marker/violin tests (cell_type preferred, else cluster)
    primary_groupby: str                 # never None — falls back to first obs column

    # Embedding availability
    has_umap: bool
    has_pca: bool

    # QC columns (None if not found)
    qc_mito_column: str | None
    qc_counts_column: str | None
    qc_genes_column: str | None

    # A few genes confirmed to exist in var_names with non-trivial expression.
    # Used by tests that need a real gene name without hardcoding biology.
    sample_genes: list[str] = field(default_factory=list)

    # Whether this is the Temple cardiac dataset (enables regression tests)
    @property
    def is_temple(self) -> bool:
        return self.name == "temple"


def _pick_sample_genes(adata: "AnnData", n: int = 5) -> list[str]:
    """Return up to n gene names confirmed to have non-zero expression in >5% of cells."""
    import scipy.sparse as sp

    threshold = max(1, int(adata.n_obs * 0.05))
    genes: list[str] = []

    X = adata.X
    if sp.issparse(X):
        # sum non-zeros per gene
        nonzero_per_gene = np.diff(X.tocsc().indptr)
    else:
        nonzero_per_gene = np.count_nonzero(X, axis=0)

    # Sort genes by non-zero cell count descending, pick top n
    candidate_indices = np.argsort(nonzero_per_gene)[::-1]
    for idx in candidate_indices:
        if nonzero_per_gene[idx] >= threshold:
            genes.append(str(adata.var_names[idx]))
        if len(genes) >= n:
            break

    return genes


def _detect_condition_column(adata: "AnnData") -> str | None:
    """Heuristic detection of a multi-condition/sample column in adata.obs."""
    candidates = ["orig.ident", "condition", "sample", "batch", "group", "donor", "timepoint"]
    for col in candidates:
        if col in adata.obs.columns:
            n_unique = adata.obs[col].nunique()
            # A condition column has ≥2 unique values but isn't high-cardinality
            if 2 <= n_unique <= 100:
                return col
    return None


def discover_profile(adata: "AnnData", name: str = "unknown") -> DatasetProfile:
    """Auto-discover a DatasetProfile from any AnnData object.

    Uses existing detection functions from src.analysis to avoid duplicating logic.
    """
    from src.analysis.cluster_resolution import detect_grouping_columns
    from src.analysis.qc_metrics import resolve_qc_metric_column
    from src.types import detect_dataset_state

    state = detect_dataset_state(adata, filename=name)
    grouping = detect_grouping_columns(adata)

    annotation_cols = grouping["annotation_columns"]
    cluster_cols = grouping["cluster_columns"]

    cell_type_col = annotation_cols[0] if annotation_cols else None
    cluster_col = cluster_cols[0] if cluster_cols else None
    condition_col = _detect_condition_column(adata)

    # Best groupby for DE/violin: prefer cell type annotation, fall back to cluster,
    # then first available obs column with <100 unique values
    if cell_type_col:
        primary_groupby = cell_type_col
    elif cluster_col:
        primary_groupby = cluster_col
    else:
        primary_groupby = next(
            (c for c in adata.obs.columns if adata.obs[c].nunique() < 100),
            adata.obs.columns[0],
        )

    return DatasetProfile(
        name=name,
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        cell_type_column=cell_type_col,
        condition_column=condition_col,
        cluster_column=cluster_col,
        primary_groupby=primary_groupby,
        has_umap=state.has_umap,
        has_pca=state.has_pca,
        qc_mito_column=resolve_qc_metric_column(adata, "mito"),
        qc_counts_column=resolve_qc_metric_column(adata, "total_counts"),
        qc_genes_column=resolve_qc_metric_column(adata, "n_genes"),
        sample_genes=_pick_sample_genes(adata),
    )


def _prepare_adata(adata: "AnnData") -> "AnnData":
    """Apply dataset-agnostic fixes to any loaded AnnData.

    Converts DataFrame embeddings (Seurat export quirk) to numpy arrays
    for any dataset — not just Temple.
    """
    for src_key, dst_key in [("UMAP", "X_umap"), ("PCA", "X_pca")]:
        if src_key in adata.obsm and dst_key not in adata.obsm:
            import pandas as pd
            if isinstance(adata.obsm[src_key], pd.DataFrame):
                adata.obsm[dst_key] = np.array(adata.obsm[src_key])
    return adata


# ---------------------------------------------------------------------------
# Module-scoped fixtures shared across all test files
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def adata():
    """Load and prepare the dataset specified by DATASET_PATH.

    Session-scoped so the h5ad file is read only once across all test files.
    """
    import anndata as ad

    if not DATASET_PATH.exists():
        pytest.skip(f"Dataset not found: {DATASET_PATH}")

    raw = ad.read_h5ad(DATASET_PATH)
    return _prepare_adata(raw)


@pytest.fixture(scope="session")
def profile(adata):
    """Discover the DatasetProfile for the loaded dataset."""
    name = DATASET_NAME or (
        "temple" if "cell_type" in adata.obs.columns
        and "Early cardiomyocyte" in adata.obs.get("cell_type", []).unique()
        else "unknown"
    )
    return discover_profile(adata, name=name)


@pytest.fixture(scope="session")
def adata_with_de(adata, profile):
    """AnnData with one-vs-rest DE pre-computed on profile.primary_groupby.

    Session-scoped so DE is only computed once regardless of how many test
    files import this fixture.
    """
    from src.analysis.differential import run_differential_expression

    adata_copy = adata.copy()
    run_differential_expression(adata_copy, groupby=profile.primary_groupby, n_genes=20)
    return adata_copy


@pytest.fixture(scope="session")
def plot_dir() -> Path:
    """Return dated plot output directory. Only created when SAVE_TEST_PLOTS=1."""
    if _SAVE_PLOTS:
        _PLOT_DIR.mkdir(parents=True, exist_ok=True)
    return _PLOT_DIR


# ---------------------------------------------------------------------------
# Helpers exported for use in test files
# ---------------------------------------------------------------------------

def valid_png(image: bytes) -> bool:
    """Return True if image is a real PNG above the minimum size threshold."""
    return (
        isinstance(image, (bytes, bytearray))
        and image[:4] == PNG_MAGIC
        and len(image) >= MIN_PLOT_BYTES
    )
