"""Streamlit chat interface for Nvwa-Lite.

Provides a chat-based UI where biologists can ask natural language
questions about their scRNA-seq data and receive visualizations.
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path

# Ensure project root is on sys.path so `src` is importable
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import scanpy as sc
import streamlit as st

from src.agent.core import create_agent
from src.agent.tools import clear_last_plot_result, get_last_plot_result
from src.plotting.styles import configure_plot_style
from src.ui.components import (
    api_key_input,
    dataset_info_panel,
    dataset_selector,
    example_queries,
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

st.set_page_config(
    page_title="Nvwa-Lite",
    page_icon="🧬",
    layout="wide",
)

# Apply publication-ready plot styles
configure_plot_style()

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

st.title("🧬 Nvwa-Lite")
st.caption("Single-cell RNA-seq visualization powered by natural language.")

api_key = api_key_input()
dataset_choice = dataset_selector()


@st.cache_resource
def load_pbmc3k():
    """Load and cache the PBMC3k processed dataset."""
    logger.info("Loading PBMC3k dataset...")
    adata = sc.datasets.pbmc3k_processed()
    logger.info("Dataset loaded: %d cells, %d genes", adata.n_obs, adata.n_vars)
    return adata


adata = load_pbmc3k()
dataset_info_panel(adata)
example_queries()


def _clean_response_text(text: str) -> str:
    """Strip base64 image markdown and data URIs that the LLM may leak."""
    text = re.sub(r"!\[.*?\]\(data:image/[^)]+\)", "", text)
    text = re.sub(r"data:image/\S+", "", text)
    return text.strip()


def _render_message(msg: dict) -> None:
    """Render a single chat message (code block, image, text)."""
    if msg.get("code"):
        with st.expander("📝 Generated Code", expanded=False):
            st.code(msg["code"], language="python")
    if msg.get("image"):
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
    st.info("Please enter your OpenAI API key in the sidebar to get started.")
    st.stop()

if prompt := st.chat_input("Ask about your data... (e.g., 'Show me the UMAP plot')"):
    # Display user message
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    # Run agent
    with st.chat_message("assistant"):
        with st.spinner("Analyzing your request..."):
            clear_last_plot_result()
            agent = create_agent(adata=adata, api_key=api_key)
            response = agent.invoke(
                user_input=prompt,
                chat_history=st.session_state.chat_history,
            )
            plot_result = get_last_plot_result()
            clean_text = _clean_response_text(response.text)

            msg_record: dict = {"role": "assistant", "content": clean_text}
            if plot_result:
                msg_record["code"] = plot_result.code
                msg_record["image"] = plot_result.image
                with st.expander("📝 Generated Code", expanded=False):
                    st.code(plot_result.code, language="python")
                st.image(plot_result.image, use_container_width=True)

            st.markdown(clean_text)

    st.session_state.messages.append(msg_record)
    st.session_state.chat_history.append(("user", prompt))
    st.session_state.chat_history.append(("assistant", clean_text))
