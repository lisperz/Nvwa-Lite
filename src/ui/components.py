"""Reusable Streamlit sidebar widgets."""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st

if TYPE_CHECKING:
    from src.types import DatasetState

UPLOAD_DIR = Path("data/uploads")
MAX_UPLOAD_MB = 2000


def file_upload_widget(
    user_id: str | None = None,
    session_id: str | None = None,
    use_direct_s3: bool = True,
) -> tuple[Path | None, str | None, Path | None]:
    """Sidebar widget for .h5ad file upload.

    Returns:
        (local_path, s3_key, temp_path)
        - Direct S3 path: (None, s3_key, temp_path) — load from temp_path
        - Proxy path: (local_path, s3_key, None) — load from local_path
        - No file: (None, None, None)
    """
    # Direct S3 path: browser uploads directly, Streamlit never touches file bytes
    if use_direct_s3 and user_id and session_id and os.getenv("S3_BUCKET_NAME"):
        try:
            from src.storage.service import S3StorageService
            from src.ui.s3_uploader import s3_direct_upload_widget

            s3_service = S3StorageService(
                bucket_name=os.getenv("S3_BUCKET_NAME"),
                region=os.getenv("AWS_REGION", "us-east-2"),
            )

            with st.sidebar:
                filename, s3_key, temp_path = s3_direct_upload_widget(
                    user_id=user_id,
                    session_id=session_id,
                    s3_service=s3_service,
                    max_size_mb=MAX_UPLOAD_MB,
                )

            if filename and s3_key and temp_path:
                return None, s3_key, temp_path

            return None, None, None

        except Exception as e:
            st.warning(f"Direct S3 upload unavailable: {e}. Falling back to proxy upload.")

    # Proxy upload fallback (file goes through Streamlit server)
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
                return None, None, None

            UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
            dest = UPLOAD_DIR / uploaded.name
            dest.write_bytes(file_bytes)

            s3_key = None
            if user_id and session_id and os.getenv("S3_BUCKET_NAME"):
                try:
                    from src.storage.service import S3StorageService
                    s3_service = S3StorageService(
                        bucket_name=os.getenv("S3_BUCKET_NAME"),
                        region=os.getenv("AWS_REGION", "us-east-2"),
                    )
                    s3_key = s3_service.upload_file(
                        user_id, session_id, file_bytes, uploaded.name, "upload"
                    )
                except Exception as e:
                    st.warning(f"S3 upload failed (using local storage): {e}")

            st.success(f"Uploaded: {uploaded.name} ({file_size_mb:.1f} MB)")
            return dest, s3_key, None

    return None, None, None


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
