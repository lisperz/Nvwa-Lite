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
from src.agent.tools import clear_plot_results, clear_table_results, get_plot_results, get_table_results, set_adata_replaced_callback, set_plot_generated_callback
from src.agent.viz_state import VisualizationState, get_viz_state
from src.auth.service import AuthService
from src.plotting.styles import configure_plot_style
from src.session.manager import SessionManager
from src.types import DatasetState, detect_dataset_state
from src.ui.components import (
    example_queries,
    file_upload_widget,
    pipeline_panel,
)
from src.ui.feedback_dialog import show_feedback_dialog
from src.ui.feedback_trigger import init_feedback_timer, mark_plot_generated, check_feedback_timer
from src.db.logger import DatabaseLogger

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

# Initialize services (cached to avoid recreation on every rerun)
@st.cache_resource
def get_auth_service():
    """Get cached authentication service."""
    tokens_file = Path("pilot_tokens.json")
    return AuthService(tokens_file=tokens_file if tokens_file.exists() else None)

@st.cache_resource
def get_session_manager():
    """Get cached session manager."""
    import redis

    redis_host = os.environ.get("REDIS_HOST", "localhost")
    redis_port = int(os.environ.get("REDIS_PORT", "6379"))
    max_concurrent = int(os.environ.get("MAX_CONCURRENT_SESSIONS", "20"))
    max_per_user = int(os.environ.get("MAX_SESSIONS_PER_USER", "2"))
    timeout_mins = int(os.environ.get("SESSION_TIMEOUT_MINUTES", "30"))

    try:
        redis_client = redis.Redis(host=redis_host, port=redis_port, decode_responses=True)
        redis_client.ping()
    except (redis.ConnectionError, redis.TimeoutError):
        redis_client = None

    return SessionManager(
        redis_client=redis_client,
        max_concurrent_sessions=max_concurrent,
        max_sessions_per_user=max_per_user,
        session_timeout_minutes=timeout_mins,
    )

auth_service = get_auth_service()
session_manager = get_session_manager()

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

# Log session info only once per session (not on every rerun)
if "session_logged" not in st.session_state:
    logger.info(f"Session active for user: {user.user_id}")
    st.session_state.session_logged = True

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

# Generate session_id early for S3 upload
if "session_id" not in st.session_state:
    import uuid
    st.session_state.session_id = str(uuid.uuid4())

uploaded_path, s3_key = file_upload_widget(user.user_id, st.session_state.session_id)


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
    try:
        adata = load_uploaded(str(uploaded_path))
    except Exception:
        logger.exception("Failed to load uploaded file: %s", uploaded_path.name)
        st.error(
            f"Could not read **{uploaded_path.name}**. The file may be corrupt or "
            "not a valid `.h5ad` (AnnData) file. Please upload a different file."
        )
        st.stop()
    ds_state = detect_dataset_state(adata, source="upload", filename=uploaded_path.name)

    # Try to create session (respects concurrency limits)
    session = session_manager.create_session(
        user_id=user.user_id,
        session_id=st.session_state.session_id,
        dataset_s3_key=s3_key or f"users/{user.user_id}/sessions/{st.session_state.session_id}/uploads/{uploaded_path.name}"
    )

    if session is None:
        st.error("Unable to create session. You may have an active session or the system is at capacity. Please try again later.")
        st.stop()

    st.session_state.adata = adata
    st.session_state.ds_state = ds_state
    st.session_state.viz_state = VisualizationState()
    st.session_state._dataset_filename = current_filename
    st.session_state.messages = []
    st.session_state.chat_history = []
    st.session_state.plot_generated = False
    st.session_state.feedback_submitted = False
    logger.info(f"New session created: {st.session_state.session_id} for user: {user.user_id}")

# Use session state adata (may have been replaced by preprocessing)
adata = st.session_state.adata
ds_state = st.session_state.ds_state

pipeline_panel(ds_state)
example_queries()


def _on_adata_replaced(new_adata) -> None:
    """Callback when preprocessing replaces the dataset."""
    st.session_state.adata = new_adata
    st.session_state.ds_state = detect_dataset_state(
        new_adata, ds_state.source, ds_state.filename,
    )


def _on_plot_generated() -> None:
    """Callback when any plot is generated."""
    st.session_state.plot_generated = True
    mark_plot_generated()


set_adata_replaced_callback(_on_adata_replaced)
set_plot_generated_callback(_on_plot_generated)


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

    # Add user message and update chat history
    st.session_state.messages.append({"role": "user", "content": prompt})
    st.session_state.chat_history.append(("user", prompt))

    # Display user message immediately
    with st.chat_message("user"):
        st.markdown(prompt)

    # Use spinner to show processing state without causing shadow effect
    with st.spinner("🤔 Thinking..."):
        try:
            clear_plot_results()
            clear_table_results()

            agent = create_agent(
                adata=adata,
                api_key=api_key,
                dataset_state=ds_state,
                viz_state=st.session_state.get("viz_state"),
                user_id=user.user_id,
                session_id=st.session_state.session_id,
            )

            response = agent.invoke(
                user_input=prompt,
                chat_history=st.session_state.chat_history[:-1],
                filename=ds_state.filename,
            )

            # Read back updated viz state after tool calls
            updated_viz = get_viz_state()
            if updated_viz is not None:
                st.session_state.viz_state = updated_viz

            plot_results = get_plot_results()
            table_results = get_table_results()
            clean_text = _clean_response_text(response.text)

            msg_idx = len(st.session_state.messages)
            msg_record: dict = {
                "role": "assistant",
                "content": clean_text,
                "idx": msg_idx,
            }

            # Build message record with plots and tables
            if plot_results:
                plots_data: list[dict] = []
                for i, pr in enumerate(plot_results):
                    plots_data.append({"code": pr.code, "image": pr.image})
                msg_record["plots"] = plots_data

            if table_results:
                tables_data: list[dict] = []
                for i, tr in enumerate(table_results):
                    tables_data.append({
                        "code": tr.code,
                        "csv_data": tr.csv_data,
                        "display_df": tr.display_df
                    })
                msg_record["tables"] = tables_data

            # Add to messages and chat history
            st.session_state.messages.append(msg_record)
            st.session_state.chat_history.append(("assistant", clean_text))

        except Exception as e:
            logger.exception("Error processing user request")
            msg_record = {
                "role": "assistant",
                "content": f"I encountered an error: {str(e)}",
                "idx": len(st.session_state.messages),
            }
            st.session_state.messages.append(msg_record)
            st.session_state.chat_history.append(("assistant", f"Error: {str(e)}"))

    # Trigger rerun to display the new messages
    st.rerun()


# ---------------------------------------------------------------------------
# Feedback Widget
# ---------------------------------------------------------------------------

# NEW-01(a): disabled for first-customer demo. Flip to True to re-enable.
# NEW-01(b) — clicking the widget mid-response kills the in-flight agent turn —
# is still open; root-cause trace deferred to a post-demo debug session.
FEEDBACK_WIDGET_ENABLED = False

# Initialize feedback state
init_feedback_timer()
if "plot_generated" not in st.session_state:
    st.session_state.plot_generated = False
if "feedback_submitted" not in st.session_state:
    st.session_state.feedback_submitted = False
if "feedback_skipped" not in st.session_state:
    st.session_state.feedback_skipped = False
if "show_feedback_dialog" not in st.session_state:
    st.session_state.show_feedback_dialog = False

# Auto-trigger feedback after 10 seconds if plot generated
if (FEEDBACK_WIDGET_ENABLED and
    st.session_state.plot_generated and
    not st.session_state.feedback_submitted and
    "session_id" in st.session_state):

    # Check timer every 2 seconds
    check_feedback_timer(delay_seconds=300)

# Show dialog if triggered
if FEEDBACK_WIDGET_ENABLED and st.session_state.get("show_feedback_dialog", False):
    show_feedback_dialog()

# Feedback icon button (bottom left corner)
if (FEEDBACK_WIDGET_ENABLED and
    st.session_state.plot_generated and
    not st.session_state.feedback_submitted and
    "session_id" in st.session_state):

    st.markdown("""
    <style>
    .feedback-icon {
        position: fixed;
        bottom: 80px;
        left: 20px;
        z-index: 999;
    }
    </style>
    """, unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="feedback-icon">', unsafe_allow_html=True)
        if st.button("💬 Feedback", key="feedback_icon_btn"):
            st.session_state.show_feedback_dialog = True
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)

# Handle feedback submission
if st.session_state.get("feedback_submitted", False) and st.session_state.get("feedback_data"):
    db_logger = DatabaseLogger()
    feedback_data = st.session_state.feedback_data
    db_logger.log_feedback(
        session_id=st.session_state.session_id,
        user_id=user.user_id,
        q1_score=feedback_data["q1_score"],
        q2_time_saved=feedback_data["q2_time_saved"],
        q3_open_text=feedback_data["q3_open_text"]
    )
    st.session_state.feedback_data = None
    st.success("Thanks! Your feedback helps us improve.")
    logger.info(f"Feedback submitted for session {st.session_state.session_id}")
    st.rerun()


