"""Timer-based feedback trigger using Streamlit fragments."""

import streamlit as st
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


def init_feedback_timer():
    """Initialize feedback timer state."""
    if "feedback_plot_time" not in st.session_state:
        st.session_state.feedback_plot_time = None
    if "feedback_timer_checked" not in st.session_state:
        st.session_state.feedback_timer_checked = False


def mark_plot_generated():
    """Mark that a plot was just generated."""
    st.session_state.feedback_plot_time = datetime.now()
    st.session_state.feedback_timer_checked = False
    logger.info(f"[FEEDBACK] Plot generated at {st.session_state.feedback_plot_time}")


@st.fragment(run_every=2)
def check_feedback_timer(delay_seconds: int = 10):
    """Check if feedback should be shown (runs every 2 seconds)."""

    logger.info(f"[FEEDBACK] Fragment running - plot_time={st.session_state.feedback_plot_time}, checked={st.session_state.feedback_timer_checked}")

    # Don't check if already triggered or no plot generated
    if st.session_state.feedback_timer_checked or st.session_state.feedback_plot_time is None:
        logger.info(f"[FEEDBACK] Skipping check - already checked or no plot")
        return

    # Calculate elapsed time
    elapsed = (datetime.now() - st.session_state.feedback_plot_time).total_seconds()
    logger.info(f"[FEEDBACK] Elapsed time: {elapsed:.1f}s / {delay_seconds}s")

    # Trigger feedback if enough time has passed
    if elapsed >= delay_seconds:
        logger.info(f"[FEEDBACK] Triggering feedback dialog!")
        st.session_state.feedback_timer_checked = True
        st.session_state.show_feedback_dialog = True
        st.rerun()
