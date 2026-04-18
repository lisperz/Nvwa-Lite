"""Snapshot of `adata.uns` after DE-mutating tool calls (T-035.b).

Two pure functions, zero coupling to agent loop:
- `snapshot_uns(adata, tool_name, prior_keys)` — build the event payload.
- `emit_uns_snapshot(logger, user_id, session_id, tool_name, adata, prior_keys)` —
  emit via `EventLogger.log_session_event` with `event_type="uns_snapshot"`.

Payload preview rules (plan Q2):
- Every value: `{"type": type(value).__name__}`.
- dict values: + `"subkeys"` — top-level keys only, no recursion.
- DataFrame / ndarray: + `"shape"`.
- Never dump contents.
"""

from __future__ import annotations

from typing import Any


def _preview(value: Any) -> dict[str, Any]:
    preview: dict[str, Any] = {"type": type(value).__name__}
    if isinstance(value, dict):
        preview["subkeys"] = list(value.keys())
    elif hasattr(value, "shape"):
        preview["shape"] = tuple(value.shape)
    return preview


def snapshot_uns(adata: Any, tool_name: str, prior_keys: list[str]) -> dict[str, Any]:
    """Return the uns_snapshot event payload for `adata` post-mutation."""
    keys_after = list(adata.uns.keys())
    prior = set(prior_keys)
    new_keys = [k for k in keys_after if k not in prior]
    payload_preview = {k: _preview(adata.uns[k]) for k in new_keys}
    return {
        "tool_name": tool_name,
        "keys_after": keys_after,
        "new_keys": new_keys,
        "payload_preview": payload_preview,
    }


def emit_uns_snapshot(
    logger: Any,
    user_id: str,
    session_id: str,
    tool_name: str,
    adata: Any,
    prior_keys: list[str],
) -> None:
    """Log an uns_snapshot event. No-ops when `logger` is None."""
    if logger is None:
        return
    metadata = snapshot_uns(adata, tool_name, prior_keys)
    logger.log_session_event(
        user_id=user_id,
        session_id=session_id,
        event_type="uns_snapshot",
        metadata=metadata,
    )
