"""Streamlit chat interface for Nvwa-Lite.

Provides a chat-based UI where biologists can ask natural language
questions about their scRNA-seq data and receive visualizations.
"""

from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path

# Ensure project root is on sys.path so `src` is importable
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import anndata as ad
import streamlit as st

from src.agent.core import create_agent
from src.agent.tools import clear_plot_results, get_plot_results, set_adata_replaced_callback
from src.plotting.styles import configure_plot_style
from src.types import DatasetState, detect_dataset_state
from src.ui.components import (
    dataset_info_panel,
    example_queries,
    file_upload_widget,
)

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "nvwa-lite.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(page_title="Nvwa-Lite", page_icon="🧬", layout="wide")
configure_plot_style()

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

st.title("🧬 Nvwa-Lite")
st.caption("Single-cell RNA-seq visualization powered by natural language.")

api_key = os.environ.get("OPENAI_API_KEY", "")
uploaded_path = file_upload_widget()


# ---------------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------------

@st.cache_resource
def load_uploaded(path: str):
    """Load and cache an uploaded .h5ad file."""
    logger.info("Loading uploaded file: %s", path)
    adata = ad.read_h5ad(path)
    logger.info("Uploaded loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)
    return adata


# Track current dataset source to detect changes
current_filename = uploaded_path.name if uploaded_path else ""

# Check if we need to reload (new file uploaded)
need_reload = (
    "adata" not in st.session_state
    or st.session_state.get("_dataset_filename") != current_filename
)

if uploaded_path is None:
    st.info("Please upload a .h5ad file in the sidebar to get started.")
    st.stop()

if need_reload:
    adata = load_uploaded(str(uploaded_path))
    ds_state = detect_dataset_state(adata, source="upload", filename=uploaded_path.name)

    st.session_state.adata = adata
    st.session_state.ds_state = ds_state
    st.session_state._dataset_filename = current_filename
    st.session_state.messages = []
    st.session_state.chat_history = []

# Use session state adata (may have been replaced by preprocessing)
adata = st.session_state.adata
ds_state = st.session_state.ds_state

dataset_info_panel(adata, state=ds_state)
example_queries()


def _on_adata_replaced(new_adata) -> None:
    """Callback when preprocessing replaces the dataset."""
    st.session_state.adata = new_adata
    st.session_state.ds_state = detect_dataset_state(
        new_adata, ds_state.source, ds_state.filename,
    )


set_adata_replaced_callback(_on_adata_replaced)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clean_response_text(text: str) -> str:
    """Strip base64 image markdown and data URIs that the LLM may leak."""
    text = re.sub(r"!\[.*?\]\(data:image/[^)]+\)", "", text)
    text = re.sub(r"data:image/\S+", "", text)
    return text.strip()


def _render_message(msg: dict) -> None:
    """Render a single chat message (code blocks, images, text)."""
    plots = msg.get("plots", [])
    for i, plot in enumerate(plots):
        with st.expander(f"📝 Generated Code ({i + 1})", expanded=False):
            st.code(plot["code"], language="python")
        st.image(plot["image"], use_container_width=True)
        st.download_button(
            label=f"Download plot {i + 1}",
            data=plot["image"],
            file_name=f"nvwa_plot_{msg.get('idx', 0)}_{i}.png",
            mime="image/png",
            key=f"dl_{msg.get('idx', 0)}_{i}",
        )
    # Legacy single-plot support
    if not plots and msg.get("image"):
        if msg.get("code"):
            with st.expander("📝 Generated Code", expanded=False):
                st.code(msg["code"], language="python")
        st.image(msg["image"], use_container_width=True)
    if msg.get("content"):
        st.markdown(msg["content"])


# ---------------------------------------------------------------------------
# Chat state
# ---------------------------------------------------------------------------

if "messages" not in st.session_state:
    st.session_state.messages = []
if "chat_history" not in st.session_state:
    st.session_state.chat_history = []

# ---------------------------------------------------------------------------
# Chat display
# ---------------------------------------------------------------------------

for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        _render_message(msg)

# ---------------------------------------------------------------------------
# Chat input
# ---------------------------------------------------------------------------

if not api_key:
    st.error("OpenAI API key not configured. Please set OPENAI_API_KEY in the environment.")
    st.stop()

if prompt := st.chat_input("Ask about your data... (e.g., 'Show me the UMAP plot')"):
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    with st.chat_message("assistant"):
        with st.spinner("Analyzing your request..."):
            clear_plot_results()
            agent = create_agent(
                adata=adata, api_key=api_key, dataset_state=ds_state,
            )
            response = agent.invoke(
                user_input=prompt,
                chat_history=st.session_state.chat_history,
            )
            plot_results = get_plot_results()
            clean_text = _clean_response_text(response.text)

            msg_idx = len(st.session_state.messages)
            msg_record: dict = {
                "role": "assistant",
                "content": clean_text,
                "idx": msg_idx,
            }

            if plot_results:
                plots_data: list[dict] = []
                for i, pr in enumerate(plot_results):
                    plots_data.append({"code": pr.code, "image": pr.image})
                    with st.expander(f"📝 Generated Code ({i + 1})", expanded=False):
                        st.code(pr.code, language="python")
                    st.image(pr.image, use_container_width=True)
                    st.download_button(
                        label=f"Download plot {i + 1}",
                        data=pr.image,
                        file_name=f"nvwa_plot_{msg_idx}_{i}.png",
                        mime="image/png",
                        key=f"dl_{msg_idx}_{i}",
                    )
                msg_record["plots"] = plots_data

            st.markdown(clean_text)

    st.session_state.messages.append(msg_record)
    st.session_state.chat_history.append(("user", prompt))
    st.session_state.chat_history.append(("assistant", clean_text))
