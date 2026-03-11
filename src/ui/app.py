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
from src.agent.tools import clear_plot_results, clear_table_results, get_plot_results, get_table_results, set_adata_replaced_callback
from src.auth.service import AuthService
from src.plotting.styles import configure_plot_style
from src.session.manager import SessionManager
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
# Authentication
# ---------------------------------------------------------------------------

# Initialize services
# Load tokens from pilot_tokens.json file (fallback to env vars if file doesn't exist)
tokens_file = Path("pilot_tokens.json")
auth_service = AuthService(tokens_file=tokens_file if tokens_file.exists() else None)
session_manager = SessionManager()

# Check authentication
if "user" not in st.session_state:
    st.title("🧬 Nvwa-Lite")
    st.caption("Single-cell RNA-seq visualization powered by natural language.")

    with st.form("auth_form"):
        st.subheader("🔐 Authentication Required")
        st.caption("Enter your pilot access token to continue.")
        token_input = st.text_input("Access Token", type="password", help="Enter the token provided by your administrator")
        submit = st.form_submit_button("Authenticate")

        if submit:
            if not token_input:
                st.error("Please enter a token.")
            else:
                user = auth_service.validate_token(token_input)
                if user:
                    st.session_state.user = user
                    st.session_state.token = token_input
                    logger.info(f"User authenticated: {user.user_id}")
                    st.success(f"Welcome, {user.email}!")
                    st.rerun()
                else:
                    st.error("Invalid token. Please check your credentials and try again.")
                    logger.warning(f"Failed authentication attempt with token: {token_input[:8]}...")

    st.stop()

# User is authenticated
user = st.session_state.user
logger.info(f"Session active for user: {user.user_id}")

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

st.title("🧬 Nvwa-Lite")
st.caption("Single-cell RNA-seq visualization powered by natural language.")

with st.sidebar:
    st.caption(f"👤 Logged in as: {user.email}")
    if st.button("🚪 Logout"):
        # Clean up session
        if "session_id" in st.session_state:
            session_manager.end_session(st.session_state.session_id)
        # Clear session state
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()

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

    # Create new session for this dataset
    import uuid
    session_id = str(uuid.uuid4())

    # Try to create session (respects concurrency limits)
    session = session_manager.create_session(
        user_id=user.user_id,
        session_id=session_id,
        dataset_s3_key=f"users/{user.user_id}/sessions/{session_id}/uploads/{uploaded_path.name}"
    )

    if session is None:
        st.error("Unable to create session. You may have an active session or the system is at capacity. Please try again later.")
        st.stop()

    st.session_state.adata = adata
    st.session_state.ds_state = ds_state
    st.session_state._dataset_filename = current_filename
    st.session_state.session_id = session_id
    st.session_state.messages = []
    st.session_state.chat_history = []
    logger.info(f"New session created: {session_id} for user: {user.user_id}")

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
    """Render a single chat message (code blocks, images, tables, text)."""
    # Render plots
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

    # Render tables
    tables = msg.get("tables", [])
    for i, table in enumerate(tables):
        with st.expander(f"📝 Generated Code ({i + 1})", expanded=False):
            st.code(table["code"], language="python")
        st.markdown(table["display_df"])
        st.download_button(
            label=f"📥 Download table {i + 1} (CSV)",
            data=table["csv_data"],
            file_name=f"nvwa_de_results_{msg.get('idx', 0)}_{i}.csv",
            mime="text/csv",
            key=f"dl_table_{msg.get('idx', 0)}_{i}",
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
    # Ensure we have an active session
    if "session_id" not in st.session_state:
        st.error("No active session. Please upload a dataset first.")
        st.stop()

    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    with st.chat_message("assistant"):
        # Create a placeholder for progressive updates
        status_placeholder = st.empty()
        status_placeholder.info("🤔 Thinking...")

        try:
            clear_plot_results()
            clear_table_results()

            status_placeholder.info("🔧 Initializing agent...")
            agent = create_agent(
                adata=adata,
                api_key=api_key,
                dataset_state=ds_state,
                user_id=user.user_id,
                session_id=st.session_state.session_id,
            )

            status_placeholder.info("💬 Processing your request...")
            response = agent.invoke(
                user_input=prompt,
                chat_history=st.session_state.chat_history,
            )

            # Clear status message
            status_placeholder.empty()

            plot_results = get_plot_results()
            table_results = get_table_results()
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

            if table_results:
                tables_data: list[dict] = []
                for i, tr in enumerate(table_results):
                    tables_data.append({
                        "code": tr.code,
                        "csv_data": tr.csv_data,
                        "display_df": tr.display_df
                    })
                    with st.expander(f"📝 Generated Code ({i + 1})", expanded=False):
                        st.code(tr.code, language="python")
                    st.markdown(tr.display_df)
                    st.download_button(
                        label=f"📥 Download table {i + 1} (CSV)",
                        data=tr.csv_data,
                        file_name=f"nvwa_de_results_{msg_idx}_{i}.csv",
                        mime="text/csv",
                        key=f"dl_table_{msg_idx}_{i}",
                    )
                msg_record["tables"] = tables_data

            st.markdown(clean_text)

        except Exception as e:
            status_placeholder.empty()
            st.error(f"An error occurred: {str(e)}")
            logger.exception("Error processing user request")
            msg_record = {
                "role": "assistant",
                "content": f"I encountered an error: {str(e)}",
                "idx": len(st.session_state.messages),
            }

    st.session_state.messages.append(msg_record)
    st.session_state.chat_history.append(("user", prompt))
    if "clean_text" in locals():
        st.session_state.chat_history.append(("assistant", clean_text))
