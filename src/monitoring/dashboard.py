"""Admin monitoring dashboard — powered by RDS PostgreSQL.

Tabs: Overview · Users · Tools · Sessions
Run via docker-compose (port 8502) or: uv run streamlit run src/monitoring/dashboard.py
"""

from __future__ import annotations

import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

import pandas as pd
import streamlit as st

import src.monitoring.analytics as analytics

# ---------------------------------------------------------------------------
# Logging
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

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Nvwa Admin Dashboard",
    page_icon="🧬",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Admin auth gate
# ---------------------------------------------------------------------------

_ADMIN_PASSWORD = os.getenv("ADMIN_PASSWORD", "")

if not _ADMIN_PASSWORD:
    st.error("ADMIN_PASSWORD environment variable is not set.")
    st.stop()

if not st.session_state.get("admin_authed"):
    st.title("🔒 Admin Login")
    pw = st.text_input("Password", type="password")
    if st.button("Login"):
        if pw == _ADMIN_PASSWORD:
            st.session_state["admin_authed"] = True
            st.rerun()
        else:
            st.error("Incorrect password.")
    st.stop()

st.title("🧬 Nvwa Admin Dashboard")

# ---------------------------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------------------------

with st.sidebar:
    st.header("Filters")
    hours = st.selectbox(
        "Time window",
        options=[1, 6, 12, 24, 48, 168],
        format_func=lambda x: f"Last {x}h" if x < 24 else f"Last {x // 24}d",
        index=3,
    )
    auto_refresh = st.toggle("Auto-refresh (30 s)", value=True)
    st.divider()
    st.caption(f"DB time: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------

tab_overview, tab_users, tab_tools, tab_sessions = st.tabs(
    ["📊 Overview", "👥 Users", "🔧 Tools", "📁 Sessions"]
)

# ============================================================
# TAB 1 — OVERVIEW
# ============================================================

with tab_overview:
    kpi = analytics.get_overview(hours)

    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Active Users", kpi["active_users"])
    c2.metric("Sessions", kpi["active_sessions"])
    c3.metric("Messages", kpi["total_messages"])
    c4.metric("Tool Calls", kpi["total_tool_calls"])
    c5.metric("Errors", kpi["total_errors"], delta_color="inverse")
    c6.metric("Est. Cost", f"${kpi['estimated_cost_usd']:.4f}")

    st.divider()

    col_left, col_right = st.columns([2, 1])

    with col_left:
        st.subheader("Hourly message activity")
        hourly = analytics.get_hourly_activity(hours)
        if hourly:
            df_h = pd.DataFrame(hourly)
            df_h["hour"] = pd.to_datetime(df_h["hour"])
            df_h = df_h.set_index("hour")
            st.bar_chart(df_h["messages"])
        else:
            st.info("No message data for this window.")

    with col_right:
        st.subheader("Token usage")
        st.metric("Input tokens", f"{kpi['input_tokens']:,}")
        st.metric("Output tokens", f"{kpi['output_tokens']:,}")
        total = kpi["input_tokens"] + kpi["output_tokens"]
        st.metric("Total tokens", f"{total:,}")
        st.caption(
            "Pricing: $0.15 / 1M input · $0.60 / 1M output (gpt-4o-mini)"
        )

# ============================================================
# TAB 2 — USERS
# ============================================================

with tab_users:
    users = analytics.get_user_breakdown(hours)

    if not users:
        st.info("No user activity in this time window.")
    else:
        df_u = pd.DataFrame(users)
        df_u["cost_usd"] = df_u.apply(
            lambda r: analytics._cost(r["input_tokens"], r["output_tokens"]), axis=1
        )
        df_u["last_seen"] = pd.to_datetime(df_u["last_seen"]).dt.strftime(
            "%Y-%m-%d %H:%M"
        )

        st.dataframe(
            df_u[
                [
                    "user_id", "sessions", "messages",
                    "tool_calls", "total_tokens", "cost_usd", "last_seen",
                ]
            ],
            use_container_width=True,
            hide_index=True,
            column_config={
                "user_id": "User",
                "sessions": st.column_config.NumberColumn("Sessions"),
                "messages": st.column_config.NumberColumn("Messages"),
                "tool_calls": st.column_config.NumberColumn("Tool Calls"),
                "total_tokens": st.column_config.NumberColumn("Tokens"),
                "cost_usd": st.column_config.NumberColumn("Cost (USD)", format="$%.4f"),
                "last_seen": "Last Seen",
            },
        )

        st.divider()
        st.subheader("Session drill-down")

        user_ids = [r["user_id"] for r in users]
        selected_user = st.selectbox("Select user", user_ids)

        if selected_user:
            sessions = analytics.get_recent_sessions(hours=hours * 4, limit=20)
            user_sessions = [s for s in sessions if s["user_id"] == selected_user]

            if not user_sessions:
                st.info("No sessions found for this user.")
            else:
                session_labels = {
                    s["session_id"]: (
                        f"{s['session_id'][:12]}…  |  {s['filename']}  |  "
                        f"{s['message_count']} msgs  |  {s['status']}"
                    )
                    for s in user_sessions
                }
                chosen = st.selectbox(
                    "Select session",
                    options=list(session_labels.keys()),
                    format_func=lambda k: session_labels[k],
                )

                if chosen:
                    msgs = analytics.get_session_messages(chosen)
                    if msgs:
                        for m in msgs:
                            role = m["role"]
                            icon = "🧑" if role == "user" else "🤖"
                            with st.chat_message(role):
                                st.markdown(f"{icon} **{role}**")
                                st.write(m["content"])
                                if m.get("response_ms"):
                                    st.caption(f"{m['response_ms']:.0f} ms")
                    else:
                        st.info("No messages recorded for this session.")

# ============================================================
# TAB 3 — TOOLS
# ============================================================

with tab_tools:
    tool_stats = analytics.get_tool_stats(hours)

    if not tool_stats:
        st.info("No tool execution data for this time window.")
    else:
        df_t = pd.DataFrame(tool_stats)
        df_t["error_rate"] = (df_t["errors"] / df_t["total_calls"] * 100).round(1)

        col_left, col_right = st.columns([2, 1])

        with col_left:
            st.subheader("Call count by tool")
            df_chart = df_t.set_index("tool_name")[["total_calls"]].sort_values(
                "total_calls"
            )
            st.bar_chart(df_chart)

        with col_right:
            st.subheader("Summary")
            st.metric("Distinct tools used", len(df_t))
            most_used = df_t.iloc[0]["tool_name"] if len(df_t) else "—"
            st.metric("Most-called tool", most_used)
            error_tools = df_t[df_t["errors"] > 0]
            st.metric("Tools with errors", len(error_tools))

        st.divider()
        st.subheader("Full tool table")
        st.dataframe(
            df_t[["tool_name", "total_calls", "errors", "error_rate", "avg_ms", "max_ms"]],
            use_container_width=True,
            hide_index=True,
            column_config={
                "tool_name": "Tool",
                "total_calls": st.column_config.NumberColumn("Calls"),
                "errors": st.column_config.NumberColumn("Errors"),
                "error_rate": st.column_config.NumberColumn("Error %", format="%.1f%%"),
                "avg_ms": st.column_config.NumberColumn("Avg ms"),
                "max_ms": st.column_config.NumberColumn("Max ms"),
            },
        )

        st.divider()
        st.subheader("🚨 Recent errors")
        errors = analytics.get_recent_errors(hours)
        if errors:
            for err in errors:
                ts = str(err["created_at"])[:19]
                label = f"{err['tool_name']}  ·  {err['user_id']}  ·  {ts}"
                with st.expander(label, expanded=False):
                    st.caption(f"Session: {err['session_id']}")
                    st.code(err["error_preview"] or "", language="text")
        else:
            st.success("No errors in this window! 🎉")

# ============================================================
# TAB 4 — SESSIONS
# ============================================================

with tab_sessions:
    sessions = analytics.get_recent_sessions(hours, limit=50)

    if not sessions:
        st.info("No sessions started in this time window.")
    else:
        df_s = pd.DataFrame(sessions)
        df_s["started_at"] = pd.to_datetime(df_s["started_at"]).dt.strftime(
            "%Y-%m-%d %H:%M"
        )
        df_s["duration"] = df_s["duration_seconds"].apply(
            lambda s: f"{s // 60}m {s % 60}s" if s is not None else "—"
        )

        st.dataframe(
            df_s[
                [
                    "started_at", "user_id", "filename",
                    "message_count", "duration", "status",
                ]
            ],
            use_container_width=True,
            hide_index=True,
            column_config={
                "started_at": "Started",
                "user_id": "User",
                "filename": "Dataset",
                "message_count": st.column_config.NumberColumn("Messages"),
                "duration": "Duration",
                "status": "Status",
            },
        )

# ---------------------------------------------------------------------------
# Auto-refresh
# ---------------------------------------------------------------------------

if auto_refresh:
    import time
    st.divider()
    st.caption("Auto-refreshing every 30 seconds…")
    time.sleep(30)
    st.rerun()
