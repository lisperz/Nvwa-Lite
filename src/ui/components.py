"""Reusable Streamlit sidebar widgets."""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st

if TYPE_CHECKING:
    from anndata import AnnData
    from src.types import DatasetState

UPLOAD_DIR = Path("data/uploads")


def api_key_input() -> str:
    """Sidebar widget for OpenAI API key input."""
    env_key = os.environ.get("OPENAI_API_KEY", "")

    with st.sidebar:
        st.subheader("OpenAI API Key")
        if env_key:
            st.success("API key loaded from environment.")
            return env_key

        key = st.text_input(
            "Enter your OpenAI API key:",
            type="password",
            placeholder="sk-...",
        )
        if not key:
            st.warning("Please provide an API key to continue.")
        return key


def file_upload_widget() -> Path | None:
    """Sidebar widget for .h5ad file upload. Returns path or None."""
    with st.sidebar:
        st.subheader("Upload Data")
        uploaded = st.file_uploader(
            "Upload a .h5ad file",
            type=["h5ad"],
            help="Upload your own scRNA-seq dataset in .h5ad format.",
        )
        if uploaded is not None:
            UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
            dest = UPLOAD_DIR / uploaded.name
            dest.write_bytes(uploaded.getvalue())
            st.success(f"Uploaded: {uploaded.name}")
            return dest
    return None


def dataset_selector() -> str:
    """Sidebar widget for dataset selection (built-in or uploaded)."""
    with st.sidebar:
        st.subheader("Dataset")
        choice = st.selectbox(
            "Select dataset:",
            options=["PBMC3k (built-in)", "Upload your own"],
            index=0,
        )
    return choice


def dataset_info_panel(adata: AnnData, state: DatasetState | None = None) -> None:
    """Sidebar widget displaying dataset summary statistics."""
    total_genes = adata.n_vars
    if adata.raw is not None:
        total_genes = adata.raw.n_vars

    with st.sidebar:
        st.subheader("Dataset Info")
        col1, col2 = st.columns(2)
        col1.metric("Cells", f"{adata.n_obs:,}")
        col2.metric("Genes", f"{total_genes:,}")

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
            st.caption(f"Pipeline: {label}")

        with st.expander("Observation keys"):
            for key in sorted(adata.obs.columns):
                st.text(f"• {key}")

        with st.expander("Common marker genes"):
            markers = ["CD3E", "CD3D", "MS4A1", "CD79A", "NKG7",
                        "CST3", "LYZ", "GNLY", "FCER1A", "PPBP"]
            raw_names = set(adata.raw.var_names) if adata.raw else set()
            all_names = set(adata.var_names) | raw_names
            available = [m for m in markers if m in all_names]
            st.text("\n".join(f"• {g}" for g in available))


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
