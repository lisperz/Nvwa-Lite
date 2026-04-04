"""Streamlit-native feedback widget using st.dialog."""

import streamlit as st
from datetime import datetime


def show_feedback_dialog():
    """Show feedback dialog as a modal overlay."""

    @st.dialog("How did this analysis go?", width="large")
    def feedback_form():
        st.markdown("### How accurate was the AI's biological interpretation? *")
        rating = st.radio(
            "Rating",
            options=[1, 2, 3, 4, 5],
            format_func=lambda x: "⭐" * x,
            horizontal=True,
            label_visibility="collapsed"
        )

        st.markdown("### Compared to doing this yourself, how much time did Nvwa save you? *")
        time_saved = st.radio(
            "Time saved",
            options=[
                "Less than 30 min",
                "30 min – 2 hours",
                "2 – 8 hours",
                "More than 1 day"
            ],
            label_visibility="collapsed"
        )

        st.markdown("### Did anything go wrong or surprise you?")
        open_text = st.text_area(
            "Optional feedback",
            placeholder="Optional feedback...",
            label_visibility="collapsed"
        )

        col1, col2 = st.columns([1, 1])
        with col1:
            if st.button("Skip", use_container_width=True):
                st.session_state.show_feedback_dialog = False
                st.session_state.feedback_timer_checked = False
                st.session_state.feedback_plot_time = datetime.now()  # Reset timer
                st.rerun()
        with col2:
            if st.button("Submit Feedback", type="primary", use_container_width=True):
                st.session_state.feedback_data = {
                    "q1_score": rating,
                    "q2_time_saved": time_saved,
                    "q3_open_text": open_text or None
                }
                st.session_state.feedback_submitted = True
                st.session_state.show_feedback_dialog = False
                st.rerun()

        return None

    return feedback_form()
