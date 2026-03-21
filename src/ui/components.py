"""Reusable Streamlit sidebar widgets."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st

if TYPE_CHECKING:
    from src.types import DatasetState

UPLOAD_DIR = Path("data/uploads")
MAX_UPLOAD_MB = 2000


def file_upload_widget() -> Path | None:
    """Sidebar widget for .h5ad file upload. Returns path or None."""
    with st.sidebar:
        st.subheader("Upload Data")
        st.caption(f"Max file size: {MAX_UPLOAD_MB} MB")
        uploaded = st.file_uploader(
            "Upload a .h5ad file",
            type=["h5ad"],
            help="Upload your own scRNA-seq dataset in .h5ad format.",
        )
        if uploaded is not None:
            file_bytes = uploaded.getvalue()
            file_size_mb = len(file_bytes) / (1024 * 1024)
            if file_size_mb > MAX_UPLOAD_MB:
                st.error(
                    f"File too large ({file_size_mb:.0f} MB). "
                    f"Maximum allowed: {MAX_UPLOAD_MB} MB."
                )
                return None
            UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
            dest = UPLOAD_DIR / uploaded.name
            dest.write_bytes(file_bytes)
            st.success(f"Uploaded: {uploaded.name} ({file_size_mb:.1f} MB)")
            return dest
    return None


def pipeline_panel(state: DatasetState | None) -> None:
    """Sidebar widget displaying the current processing pipeline status."""
    with st.sidebar:
        st.subheader("Pipeline")
        if state is not None:
            steps: list[str] = []
            if state.is_normalized:
                steps.append("Normalized")
            if state.has_pca:
                steps.append("PCA")
            if state.has_umap:
                steps.append("UMAP")
            if state.has_clustering:
                steps.append(f"Clustered ({state.cluster_key})")
            if state.has_de_results:
                steps.append("DE")
            label = " → ".join(steps) if steps else "Raw / unprocessed"
            st.caption(label)
        else:
            st.caption("Raw / unprocessed")


def example_queries() -> None:
    """Sidebar widget with example queries for quick start."""
    with st.sidebar:
        st.subheader("Example Queries")
        examples = [
            "What's in this dataset?",
            "Preprocess the data.",
            "Show me the UMAP plot colored by cell type.",
            "Show a violin plot for MS4A1 across all clusters.",
            "Run differential expression.",
            "Show volcano plot for cluster 0.",
            "Show a heatmap for CD3E, MS4A1, NKG7, LYZ.",
        ]
        for ex in examples:
            st.caption(f"__{ex}__")
