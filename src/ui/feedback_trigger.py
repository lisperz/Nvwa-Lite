"""Timer-based feedback trigger using JavaScript timer."""

import streamlit as st
import streamlit.components.v1 as components
from datetime import datetime


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


def check_feedback_timer(delay_seconds: int = 10):
    """Check if feedback should be shown using JavaScript timer.

    This approach uses a JavaScript setTimeout to trigger a page reload
    after the specified delay, which works reliably in production.
    """
    # Don't check if already triggered or no plot generated
    if st.session_state.feedback_timer_checked or st.session_state.feedback_plot_time is None:
        return

    # Calculate elapsed time
    elapsed = (datetime.now() - st.session_state.feedback_plot_time).total_seconds()

    # If enough time has passed, trigger feedback dialog
    if elapsed >= delay_seconds:
        st.session_state.feedback_timer_checked = True
        st.session_state.show_feedback_dialog = True
        st.rerun()
    else:
        # Not enough time yet - inject JavaScript timer to reload page after remaining time
        remaining_ms = int((delay_seconds - elapsed) * 1000)

        # Use JavaScript setTimeout to trigger page reload
        components.html(
            f"""
            <script>
                setTimeout(function() {{
                    window.parent.location.reload();
                }}, {remaining_ms});
            </script>
            """,
            height=0,
        )
