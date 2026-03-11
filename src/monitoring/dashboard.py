"""Monitoring dashboard for Nvwa-Lite MVP.

Real-time metrics and analytics for pilot deployment monitoring.
Displays KPIs, user activity, tool usage, and error tracking.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Ensure project root is on sys.path
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import streamlit as st

from src.monitoring.analytics import AnalyticsService

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "dashboard.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Nvwa-Lite Dashboard",
    page_icon="📊",
    layout="wide",
)

st.title("📊 Nvwa-Lite Monitoring Dashboard")
st.caption("Real-time metrics for MVP pilot deployment")

# ---------------------------------------------------------------------------
# Analytics service
# ---------------------------------------------------------------------------

analytics = AnalyticsService(log_dir=LOG_DIR)

# ---------------------------------------------------------------------------
# Time range selector
# ---------------------------------------------------------------------------

time_range = st.selectbox(
    "Time Range",
    options=[1, 6, 12, 24, 48, 168],  # hours
    format_func=lambda x: f"Last {x} hours" if x < 24 else f"Last {x // 24} days",
    index=3,  # Default to 24 hours
)

st.divider()

# ---------------------------------------------------------------------------
# Key metrics
# ---------------------------------------------------------------------------

col1, col2, col3, col4 = st.columns(4)

active_users = analytics.get_active_users(hours=time_range)
total_sessions = analytics.get_total_sessions(hours=time_range)
error_count = analytics.get_error_count(hours=time_range)
token_usage = analytics.get_total_token_usage(hours=time_range)

with col1:
    st.metric("Active Users", len(active_users))

with col2:
    st.metric("Total Sessions", total_sessions)

with col3:
    st.metric("Errors", error_count, delta_color="inverse")

with col4:
    st.metric("Total Tokens", f"{token_usage['total']:,}")

st.divider()

# ---------------------------------------------------------------------------
# Tool usage and performance
# ---------------------------------------------------------------------------

col1, col2 = st.columns(2)

with col1:
    st.subheader("🔧 Tool Usage")
    tool_stats = analytics.get_tool_usage_stats(hours=time_range)

    if tool_stats:
        # Sort by usage count
        sorted_tools = sorted(tool_stats.items(), key=lambda x: x[1], reverse=True)

        for tool_name, count in sorted_tools:
            st.metric(tool_name, count)
    else:
        st.info("No tool usage data available for this time range.")

with col2:
    st.subheader("⚡ Performance")

    avg_response_time = analytics.get_average_response_time(hours=time_range)
    st.metric("Avg Response Time", f"{avg_response_time:.2f}s")

    st.caption("Token Breakdown")
    st.metric("Prompt Tokens", f"{token_usage['prompt']:,}")
    st.metric("Completion Tokens", f"{token_usage['completion']:,}")

st.divider()

# ---------------------------------------------------------------------------
# User activity
# ---------------------------------------------------------------------------

st.subheader("👥 User Activity")

user_activity = analytics.get_user_activity(hours=time_range)

if user_activity:
    # Display as table
    st.dataframe(
        user_activity,
        use_container_width=True,
        hide_index=True,
        column_config={
            "user_id": "User ID",
            "messages": st.column_config.NumberColumn("Messages", format="%d"),
            "tools": st.column_config.NumberColumn("Tool Calls", format="%d"),
            "sessions": st.column_config.NumberColumn("Sessions", format="%d"),
        },
    )
else:
    st.info("No user activity data available for this time range.")

st.divider()

# ---------------------------------------------------------------------------
# Recent errors
# ---------------------------------------------------------------------------

st.subheader("🚨 Recent Errors")

recent_errors = analytics.get_recent_errors(limit=10)

if recent_errors:
    for error in recent_errors:
        with st.expander(
            f"⚠️ {error['tool_name']} - {error['timestamp'][:19]}",
            expanded=False,
        ):
            st.text(f"User: {error['user_id']}")
            st.text(f"Session: {error['session_id']}")
            st.code(error['error'], language="text")
else:
    st.success("No errors in this time range! 🎉")

# ---------------------------------------------------------------------------
# Auto-refresh
# ---------------------------------------------------------------------------

st.divider()
st.caption("Dashboard auto-refreshes every 30 seconds")

# Auto-refresh every 30 seconds
import time
time.sleep(30)
st.rerun()

