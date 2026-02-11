"""Reusable Streamlit sidebar widgets."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import streamlit as st

if TYPE_CHECKING:
    from anndata import AnnData


def api_key_input() -> str:
    """Sidebar widget for OpenAI API key input.

    Falls back to OPENAI_API_KEY environment variable.
    """
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


def dataset_selector() -> str:
    """Sidebar widget for dataset selection."""
    with st.sidebar:
        st.subheader("Dataset")
        choice = st.selectbox(
            "Select dataset:",
            options=["PBMC3k (built-in)"],
            index=0,
        )
    return choice


def dataset_info_panel(adata: AnnData) -> None:
    """Sidebar widget displaying dataset summary statistics."""
    total_genes = adata.n_vars
    if adata.raw is not None:
        total_genes = adata.raw.n_vars

    with st.sidebar:
        st.subheader("Dataset Info")
        col1, col2 = st.columns(2)
        col1.metric("Cells", f"{adata.n_obs:,}")
        col2.metric("Genes", f"{total_genes:,}")

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
            "Show me the UMAP plot colored by cell type.",
            "Visualize the expression of CD3E on the UMAP.",
            "Show a violin plot for MS4A1 across all clusters.",
            "Create a dot plot for CD3E, MS4A1, NKG7.",
            "What genes are available in this dataset?",
        ]
        for ex in examples:
            st.caption(f"💬 _{ex}_")
