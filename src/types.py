"""Dataset state tracking for the Nvwa-Lite pipeline.

Tracks what preprocessing steps have been applied to an AnnData object,
enabling the agent to make informed decisions about what tools to suggest.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


@dataclass
class DatasetState:
    """Tracks the processing state of a loaded dataset."""

    source: str = "unknown"
    filename: str = ""
    has_raw_counts: bool = False
    is_normalized: bool = False
    has_hvg: bool = False
    has_pca: bool = False
    has_umap: bool = False
    has_clustering: bool = False
    has_de_results: bool = False
    cluster_key: str = ""
    n_cells: int = 0
    n_genes: int = 0

    def summary(self) -> str:
        """Human-readable summary of what processing has been done."""
        steps: list[str] = []
        if self.is_normalized:
            steps.append("normalized")
        if self.has_hvg:
            steps.append("HVG selected")
        if self.has_pca:
            steps.append("PCA computed")
        if self.has_umap:
            steps.append("UMAP computed")
        if self.has_clustering:
            steps.append(f"clustered (key={self.cluster_key})")
        if self.has_de_results:
            steps.append("DE results available")
        done = ", ".join(steps) if steps else "none"
        return (
            f"Source: {self.source} | File: {self.filename}\n"
            f"Cells: {self.n_cells:,} | Genes: {self.n_genes:,}\n"
            f"Processing: {done}"
        )


def detect_dataset_state(
    adata: AnnData,
    source: str = "unknown",
    filename: str = "",
) -> DatasetState:
    """Inspect an AnnData object and determine its processing state.

    Uses a hierarchical detection strategy:
    1. Check for X_umap -> ANALYSIS COMPLETE (ready for all visualizations)
    2. Check for X_pca -> PREPROCESSED (needs UMAP/clustering)
    3. Check for clustering labels in obs -> PARTIALLY PREPROCESSED (can still plot)
    4. Otherwise -> RAW (needs full preprocessing)
    """
    state = DatasetState(
        source=source,
        filename=filename,
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
    )
    state.has_raw_counts = adata.raw is not None
    state.is_normalized = "log1p" in getattr(adata.uns, "keys", lambda: [])()
    state.has_hvg = "highly_variable" in adata.var.columns
    state.has_pca = "X_pca" in adata.obsm
    state.has_umap = "X_umap" in adata.obsm

    # Detect clustering key (self-healing logic)
    for key in ("leiden", "louvain"):
        if key in adata.obs.columns:
            state.has_clustering = True
            state.cluster_key = key
            break

    state.has_de_results = "rank_genes_groups" in adata.uns

    # Self-healing: If clustering labels exist but is_normalized is False,
    # assume the data is actually preprocessed (just missing log1p marker)
    if state.has_clustering and not state.is_normalized:
        # Check if data looks normalized (values in reasonable range)
        import numpy as np
        if hasattr(adata.X, 'max'):
            max_val = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
            # If max value is < 20, likely log-normalized
            if max_val < 20:
                state.is_normalized = True

    return state
