"""Reusable Streamlit sidebar widgets."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st

if TYPE_CHECKING:
    from src.types import DatasetState

logger = logging.getLogger(__name__)

UPLOAD_DIR = Path("data/uploads")
PRELOADED_DIR = Path("data/preloaded")
MAX_UPLOAD_MB = 2000
PRELOADED_DATASETS_FILE = Path("preloaded_datasets.json")


@st.cache_resource
def _load_preloaded_config() -> dict:
    """Load and cache preloaded_datasets.json. Returns empty dict if missing."""
    if not PRELOADED_DATASETS_FILE.exists():
        return {}
    try:
        return json.loads(PRELOADED_DATASETS_FILE.read_text())
    except Exception as e:
        logger.warning("Failed to load preloaded_datasets.json: %s", e)
        return {}


def preloaded_dataset_widget(user_id: str) -> Path | None:
    """Sidebar widget for pre-registered datasets.

    Reads preloaded_datasets.json and renders an 'Available Datasets' section
    if the user has any registered datasets. Returns the selected local Path
    when the user clicks 'Use This Dataset', otherwise None.

    The local path is never shown in the UI — only display_name is shown.

    Once a dataset is selected, it persists in session_state until the user
    explicitly selects a different dataset or uploads a file.
    """
    config = _load_preloaded_config()
    datasets = config.get(user_id, [])
    if not datasets:
        return None

    # Check if we already have a preloaded dataset active
    if "preloaded_dataset_path" in st.session_state:
        return st.session_state.preloaded_dataset_path

    with st.sidebar:
        st.subheader("Available Datasets")

        options = {d["display_name"]: d for d in datasets}
        selected_name = st.radio(
            "Select a dataset",
            options=list(options.keys()),
            label_visibility="collapsed",
        )
        selected = options[selected_name]

        if selected.get("description"):
            st.caption(selected["description"])

        if st.button("Use This Dataset", key="use_preloaded_btn"):
            path = PRELOADED_DIR / selected["filename"]

            if not path.exists():
                st.error(
                    f"Dataset file not found on server. "
                    f"Please contact your administrator."
                )
                logger.error("Preloaded dataset missing: %s", path)
                return None

            min_mb = selected.get("file_size_min_mb", 0)
            actual_mb = path.stat().st_size / (1024 * 1024)
            if min_mb and actual_mb < min_mb:
                st.error(
                    f"Dataset file appears to be the wrong version "
                    f"({actual_mb:.0f} MB, expected >= {min_mb} MB). "
                    f"Please contact your administrator."
                )
                logger.error(
                    "Preloaded dataset size mismatch: %s is %.0f MB, expected >= %.0f MB",
                    path, actual_mb, min_mb,
                )
                return None

            # Persist in session state so it survives reruns
            st.session_state.preloaded_dataset_path = path
            return path

    return None


def file_upload_widget(user_id: str | None = None, session_id: str | None = None) -> tuple[Path | None, str | None]:
    """Sidebar widget for .h5ad file upload. Returns (local_path, s3_key) or (None, None)."""
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
                return None, None

            # Save locally first (for immediate loading)
            UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
            dest = UPLOAD_DIR / uploaded.name
            dest.write_bytes(file_bytes)

            # Upload to S3 if configured and user/session provided
            s3_key = None
            if user_id and session_id and os.getenv("S3_BUCKET_NAME"):
                try:
                    from src.storage.service import S3StorageService
                    s3_service = S3StorageService(
                        bucket_name=os.getenv("S3_BUCKET_NAME"),
                        region=os.getenv("AWS_REGION", "us-east-2")
                    )
                    s3_key = s3_service.upload_file(user_id, session_id, file_bytes, uploaded.name, "upload")
                except Exception as e:
                    st.warning(f"S3 upload failed (using local storage): {e}")

            st.success(f"Uploaded: {uploaded.name} ({file_size_mb:.1f} MB)")
            # Clear preloaded selection so the upload takes priority
            st.session_state.pop("preloaded_dataset_path", None)
            return dest, s3_key
    return None, None


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
