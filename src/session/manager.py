"""Session management with concurrency control.

Manages user sessions with Redis-backed state tracking and enforces
concurrency limits (max 1 session per user, max 20 concurrent system-wide).

Session count is derived from live Redis keys (SCAN session:*) rather than
a separate counter, so TTL expiry automatically reduces the count — no drift.
"""

from __future__ import annotations

import json
import logging
import threading
from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

try:
    import redis
    REDIS_AVAILABLE = True
except ImportError:
    REDIS_AVAILABLE = False

logger = logging.getLogger(__name__)


class SessionStatus(Enum):
    """Session lifecycle states."""

    ACTIVE = "active"
    IDLE = "idle"
    COMPLETED = "completed"
    EXPIRED = "expired"


@dataclass
class Session:
    """Represents a user analysis session."""

    session_id: str
    user_id: str
    dataset_s3_key: str
    status: SessionStatus
    created_at: datetime
    last_activity: datetime
    message_count: int = 0
    tool_calls: list[str] = field(default_factory=list)


class SessionManager:
    """Manage user sessions with concurrency control.

    For the pilot, enforces:
    - Max 1 active session per user
    - Max 2 concurrent sessions system-wide
    - 30-minute session timeout (auto-cleanup)
    """

    def __init__(
        self,
        redis_client: redis.Redis | None = None,
        max_concurrent_sessions: int = 20,
        max_sessions_per_user: int = 2,
        session_timeout_minutes: int = 30,
    ):
        """Initialize the session manager.

        Args:
            redis_client: Redis client for state storage. If None, uses in-memory fallback.
            max_concurrent_sessions: Maximum concurrent sessions system-wide (default 20).
            max_sessions_per_user: Maximum sessions per user (default 2).
            session_timeout_minutes: Session timeout in minutes (default 30).
        """
        if redis_client is None and REDIS_AVAILABLE:
            # Try to connect to local Redis
            try:
                redis_client = redis.Redis(
                    host="localhost",
                    port=6379,
                    decode_responses=True,
                )
                redis_client.ping()  # Test connection
            except (redis.ConnectionError, redis.TimeoutError):
                redis_client = None

        self.redis = redis_client
        self.max_concurrent_sessions = max_concurrent_sessions
        self.max_sessions_per_user = max_sessions_per_user
        self.session_timeout_minutes = session_timeout_minutes

        # Fallback in-memory storage if Redis unavailable
        if self.redis is None:
            self._memory_sessions: dict[str, dict] = {}
            self._memory_user_sessions: dict[str, set[str]] = {}

        # Reentrant lock — cleanup endpoint thread + main Streamlit thread both call end_session
        self._lock = threading.RLock()

    def create_session(
        self,
        user_id: str,
        session_id: str,
        dataset_s3_key: str,
    ) -> Session | None:
        """Create a new session if concurrency limits allow.

        Before checking limits, cleans orphaned set entries (BUG #22 fix):
        a hard browser reload loses st.session_state but leaves the old Redis
        entry alive. Orphaned entries (key expired but set reference remains)
        are cleaned automatically. Live sessions are NOT evicted — the user
        is hard-capped at max_sessions_per_user.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            dataset_s3_key: S3 key of the uploaded dataset.

        Returns:
            Session object if created, None if limits exceeded.
        """
        with self._lock:
            # Clean orphaned user session entries (BUG #22 fix)
            self._cleanup_expired_user_sessions(user_id, current_session_id=session_id)

            # Clean stale sessions (heartbeat fallback — last_activity > 45s)
            self._cleanup_stale_sessions(user_id)

            # Check user's active sessions (after cleanup)
            user_session_count = self._get_user_session_count(user_id)
            if user_session_count >= self.max_sessions_per_user:
                return None  # User at session limit

            # Check system-wide concurrency (derived from live keys — BUG #21 fix)
            active_count = self._get_active_count()
            if active_count >= self.max_concurrent_sessions:
                return None  # System at capacity

            # Create session
            session = Session(
                session_id=session_id,
                user_id=user_id,
                dataset_s3_key=dataset_s3_key,
                status=SessionStatus.ACTIVE,
                created_at=datetime.utcnow(),
                last_activity=datetime.utcnow(),
            )

            self._store_session(session)
            return session

    def get_session(self, session_id: str) -> Session | None:
        """Get session by ID.

        Args:
            session_id: Session identifier.

        Returns:
            Session object if found, None otherwise.
        """
        if self.redis:
            session_key = f"session:{session_id}"
            data = self.redis.hgetall(session_key)
            if not data:
                return None

            return Session(
                session_id=session_id,
                user_id=data["user_id"],
                dataset_s3_key=data["dataset_s3_key"],
                status=SessionStatus(data["status"]),
                created_at=datetime.fromisoformat(data["created_at"]),
                last_activity=datetime.fromisoformat(data["last_activity"]),
                message_count=int(data.get("message_count", 0)),
                tool_calls=json.loads(data.get("tool_calls", "[]")),
            )
        else:
            data = self._memory_sessions.get(session_id)
            if not data:
                return None

            return Session(
                session_id=session_id,
                user_id=data["user_id"],
                dataset_s3_key=data["dataset_s3_key"],
                status=SessionStatus(data["status"]),
                created_at=datetime.fromisoformat(data["created_at"]),
                last_activity=datetime.fromisoformat(data["last_activity"]),
                message_count=data.get("message_count", 0),
                tool_calls=data.get("tool_calls", []),
            )

    def update_activity(self, session_id: str, tool_name: str | None = None) -> None:
        """Update session activity timestamp.

        Args:
            session_id: Session identifier.
            tool_name: Optional tool name to record.
        """
        if self.redis:
            session_key = f"session:{session_id}"
            self.redis.hset(session_key, "last_activity", datetime.utcnow().isoformat())
            self.redis.hincrby(session_key, "message_count", 1)

            if tool_name:
                tool_calls = json.loads(self.redis.hget(session_key, "tool_calls") or "[]")
                tool_calls.append(tool_name)
                self.redis.hset(session_key, "tool_calls", json.dumps(tool_calls))

            # Refresh expiration
            self.redis.expire(session_key, self.session_timeout_minutes * 60)
        else:
            if session_id in self._memory_sessions:
                self._memory_sessions[session_id]["last_activity"] = datetime.utcnow().isoformat()
                self._memory_sessions[session_id]["message_count"] += 1

                if tool_name:
                    self._memory_sessions[session_id]["tool_calls"].append(tool_name)

    def end_session(self, session_id: str) -> None:
        """End a session and free up resources.

        Args:
            session_id: Session identifier.
        """
        with self._lock:
            session = self.get_session(session_id)
            if not session:
                return

            if self.redis:
                session_key = f"session:{session_id}"
                user_sessions_key = f"user:{session.user_id}:sessions"

                # Remove from user's sessions
                self.redis.srem(user_sessions_key, session_id)

                # Mark as completed
                self.redis.hset(session_key, "status", SessionStatus.COMPLETED.value)

                # Keep for 24 hours for analytics
                self.redis.expire(session_key, 86400)
            else:
                if session_id in self._memory_sessions:
                    self._memory_sessions[session_id]["status"] = SessionStatus.COMPLETED.value
                    user_id = self._memory_sessions[session_id]["user_id"]
                    if user_id in self._memory_user_sessions:
                        self._memory_user_sessions[user_id].discard(session_id)

    def heartbeat(self, session_id: str) -> None:
        """Update last_activity for a session (heartbeat fallback).

        Called every ~15s by the Streamlit fragment. If the tab dies and
        sendBeacon fails, _cleanup_stale_sessions will end this session
        when last_activity is older than 45 seconds.

        Args:
            session_id: Session identifier.
        """
        with self._lock:
            if self.redis:
                session_key = f"session:{session_id}"
                if self.redis.exists(session_key):
                    self.redis.hset(session_key, "last_activity", datetime.utcnow().isoformat())
            else:
                if session_id in self._memory_sessions:
                    self._memory_sessions[session_id]["last_activity"] = datetime.utcnow().isoformat()

    def _cleanup_stale_sessions(self, user_id: str) -> None:
        """End ACTIVE sessions with no heartbeat for > 45 seconds (fallback only).

        Called inside create_session() before limit checks. Only fires when
        sendBeacon cleanup failed — e.g. browser crash or network drop.
        Never evicts sessions that are still receiving heartbeats.

        Args:
            user_id: User whose stale sessions to clean.
        """
        stale_threshold = datetime.utcnow() - timedelta(seconds=45)

        if self.redis:
            user_sessions_key = f"user:{user_id}:sessions"
            session_ids = self.redis.smembers(user_sessions_key)
            for sid in session_ids:
                data = self.redis.hgetall(f"session:{sid}")
                if not data:
                    continue
                if data.get("status") != SessionStatus.ACTIVE.value:
                    continue
                last_activity_str = data.get("last_activity")
                if last_activity_str:
                    try:
                        last_activity = datetime.fromisoformat(last_activity_str)
                        if last_activity < stale_threshold:
                            logger.info("Cleaning stale session %s (last_activity=%s)", sid, last_activity_str)
                            self.end_session(sid)
                    except ValueError:
                        pass
        else:
            for sid, data in list(self._memory_sessions.items()):
                if data.get("user_id") != user_id:
                    continue
                if data.get("status") != SessionStatus.ACTIVE.value:
                    continue
                last_activity_str = data.get("last_activity")
                if last_activity_str:
                    try:
                        last_activity = datetime.fromisoformat(last_activity_str)
                        if last_activity < stale_threshold:
                            logger.info("Cleaning stale session %s (last_activity=%s)", sid, last_activity_str)
                            self.end_session(sid)
                    except ValueError:
                        pass

    def get_queue_position(self, user_id: str) -> int | None:
        """Get user's position in the queue (if waiting).

        Args:
            user_id: User identifier.

        Returns:
            Queue position (1-indexed) or None if not queued.
        """
        # For Phase 1, we don't implement a queue - just return None
        # Phase 2 will add proper job queue with positions
        return None

    def _has_active_session(self, user_id: str) -> bool:
        """Check if user has an active session."""
        return self._get_user_session_count(user_id) > 0

    def _get_user_session_count(self, user_id: str) -> int:
        """Get number of active sessions for a user."""
        if self.redis:
            user_sessions_key = f"user:{user_id}:sessions"
            # Clean up orphaned session IDs
            session_ids = self.redis.smembers(user_sessions_key)
            for session_id in session_ids:
                if not self.redis.exists(f"session:{session_id}"):
                    self.redis.srem(user_sessions_key, session_id)
            return self.redis.scard(user_sessions_key)
        else:
            return len(self._memory_user_sessions.get(user_id, set()))

    def _get_active_count(self) -> int:
        """Get current number of active sessions.

        BUG #21 fix: count live session:* keys via SCAN instead of reading a
        separate counter. When Redis TTL expires a session key the count drops
        automatically — no drift possible.
        """
        if self.redis:
            count = 0
            cursor = 0
            while True:
                cursor, keys = self.redis.scan(cursor, match="session:*", count=100)
                count += len(keys)
                if cursor == 0:
                    break
            return count
        else:
            # In-memory fallback: count sessions not yet marked completed
            return sum(
                1 for s in self._memory_sessions.values()
                if s.get("status") == SessionStatus.ACTIVE.value
            )

    def _cleanup_expired_user_sessions(self, user_id: str, current_session_id: str) -> None:
        """Remove orphaned set entries for sessions whose Redis key has expired.

        BUG #22 fix: a hard browser reload clears st.session_state so the old
        session_id is lost, but if its TTL hasn't elapsed the Redis key still
        exists and counts against the user's slot limit. This method only
        removes set entries whose Redis key is already gone — it never evicts
        live sessions. The user is hard-capped at max_sessions_per_user.

        Args:
            user_id: User whose expired set entries to clean.
            current_session_id: The new session_id being created (skip if found).
        """
        if not self.redis:
            return

        user_sessions_key = f"user:{user_id}:sessions"
        session_ids = self.redis.smembers(user_sessions_key)
        for sid in session_ids:
            if sid == current_session_id:
                continue
            if not self.redis.exists(f"session:{sid}"):
                # Key already TTL-expired — remove the stale set reference
                self.redis.srem(user_sessions_key, sid)

    def _store_session(self, session: Session) -> None:
        """Store session in Redis or memory."""
        if self.redis:
            session_key = f"session:{session.session_id}"
            user_sessions_key = f"user:{session.user_id}:sessions"

            self.redis.hset(
                session_key,
                mapping={
                    "user_id": session.user_id,
                    "dataset_s3_key": session.dataset_s3_key,
                    "status": session.status.value,
                    "created_at": session.created_at.isoformat(),
                    "last_activity": session.last_activity.isoformat(),
                    "message_count": session.message_count,
                    "tool_calls": json.dumps(session.tool_calls),
                },
            )

            # Add to user's sessions
            self.redis.sadd(user_sessions_key, session.session_id)

            # Set expiration
            self.redis.expire(session_key, self.session_timeout_minutes * 60)
        else:
            self._memory_sessions[session.session_id] = {
                "user_id": session.user_id,
                "dataset_s3_key": session.dataset_s3_key,
                "status": session.status.value,
                "created_at": session.created_at.isoformat(),
                "last_activity": session.last_activity.isoformat(),
                "message_count": session.message_count,
                "tool_calls": session.tool_calls,
            }

            if session.user_id not in self._memory_user_sessions:
                self._memory_user_sessions[session.user_id] = set()
            self._memory_user_sessions[session.user_id].add(session.session_id)

