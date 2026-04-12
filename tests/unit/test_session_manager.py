"""Tests for session manager BUG #21 and #22 fixes.

Note: These tests use the in-memory fallback (no Redis). The bugs are Redis-specific:
- BUG #21 is fixed by SCAN-based counting in _get_active_count() for Redis
- BUG #22 is fixed by _evict_stale_user_sessions() which is a no-op in in-memory mode

For Redis-specific behavior, manual testing or integration tests are required.
"""

import pytest
from src.session.manager import SessionManager, SessionStatus


class TestSessionCounterBug21:
    """Test BUG #21 fix: SCAN-based active count (Redis) vs in-memory fallback."""

    def test_active_count_empty(self):
        """When no sessions exist, count should be 0."""
        mgr = SessionManager(redis_client=None, max_concurrent_sessions=20)
        assert mgr._get_active_count() == 0

    def test_active_count_after_create(self):
        """Active count should increment after session creation."""
        mgr = SessionManager(redis_client=None, max_concurrent_sessions=20)
        session = mgr.create_session(
            user_id="user1",
            session_id="sid1",
            dataset_s3_key="s3://bucket/file.h5ad"
        )
        assert session is not None
        assert mgr._get_active_count() == 1

    def test_active_count_after_end_session(self):
        """Active count should decrement after ending session."""
        mgr = SessionManager(redis_client=None, max_concurrent_sessions=20)
        session = mgr.create_session(
            user_id="user1",
            session_id="sid1",
            dataset_s3_key="s3://bucket/file.h5ad"
        )
        assert mgr._get_active_count() == 1
        mgr.end_session("sid1")
        assert mgr._get_active_count() == 0

    def test_active_count_multiple_sessions(self):
        """Count tracks all active sessions independently."""
        mgr = SessionManager(redis_client=None, max_concurrent_sessions=20)
        for i in range(3):
            mgr.create_session(
                user_id=f"user{i}",
                session_id=f"sid{i}",
                dataset_s3_key=f"s3://bucket/file{i}.h5ad"
            )
        assert mgr._get_active_count() == 3
        mgr.end_session("sid1")
        assert mgr._get_active_count() == 2
        mgr.end_session("sid0")
        assert mgr._get_active_count() == 1

    def test_system_capacity_limit(self):
        """System-wide max_concurrent_sessions limit is enforced."""
        mgr = SessionManager(redis_client=None, max_concurrent_sessions=2)
        # Create 2 sessions (at capacity)
        s1 = mgr.create_session("user1", "sid1", "s3://bucket/f1.h5ad")
        s2 = mgr.create_session("user2", "sid2", "s3://bucket/f2.h5ad")
        assert s1 is not None and s2 is not None
        assert mgr._get_active_count() == 2

        # Try to create 3rd (should fail)
        s3 = mgr.create_session("user3", "sid3", "s3://bucket/f3.h5ad")
        assert s3 is None
        assert mgr._get_active_count() == 2


class TestUserSessionLimit:
    """Test per-user session limits — max_sessions_per_user=2 is preserved."""

    def test_user_can_have_two_concurrent_sessions(self):
        """User can hold up to max_sessions_per_user=2 valid sessions."""
        mgr = SessionManager(redis_client=None, max_sessions_per_user=2)
        s1 = mgr.create_session("user1", "sid1", "s3://bucket/f1.h5ad")
        s2 = mgr.create_session("user1", "sid2", "s3://bucket/f2.h5ad")
        assert s1 is not None
        assert s2 is not None
        assert mgr._get_user_session_count("user1") == 2

    def test_user_third_session_rejected(self):
        """A third session for the same user is rejected when at limit."""
        mgr = SessionManager(redis_client=None, max_sessions_per_user=2)
        mgr.create_session("user1", "sid1", "s3://bucket/f1.h5ad")
        mgr.create_session("user1", "sid2", "s3://bucket/f2.h5ad")
        s3 = mgr.create_session("user1", "sid3", "s3://bucket/f3.h5ad")
        # Hard cap: third session rejected, no auto-eviction
        assert s3 is None
        assert mgr._get_user_session_count("user1") == 2

    def test_existing_sessions_not_evicted_on_third_attempt(self):
        """When a third session is rejected, existing sessions remain ACTIVE."""
        mgr = SessionManager(redis_client=None, max_sessions_per_user=2)
        s1 = mgr.create_session("user1", "sid1", "s3://bucket/f1.h5ad")
        s2 = mgr.create_session("user1", "sid2", "s3://bucket/f2.h5ad")

        # Attempt third session (should be rejected)
        s3 = mgr.create_session("user1", "sid3", "s3://bucket/f3.h5ad")
        assert s3 is None

        # Verify first two sessions are still ACTIVE (not COMPLETED)
        s1_check = mgr.get_session("sid1")
        s2_check = mgr.get_session("sid2")
        assert s1_check is not None and s1_check.status == SessionStatus.ACTIVE
        assert s2_check is not None and s2_check.status == SessionStatus.ACTIVE

    def test_user_session_count_after_end(self):
        """Ending a session frees the user slot."""
        mgr = SessionManager(redis_client=None, max_sessions_per_user=2)
        mgr.create_session("user1", "sid1", "s3://bucket/f1.h5ad")
        mgr.create_session("user1", "sid2", "s3://bucket/f2.h5ad")
        assert mgr._get_user_session_count("user1") == 2

        mgr.end_session("sid1")
        assert mgr._get_user_session_count("user1") == 1

        # Now can create a new one
        s3 = mgr.create_session("user1", "sid3", "s3://bucket/f3.h5ad")
        assert s3 is not None
        assert mgr._get_user_session_count("user1") == 2

    def test_different_users_independent(self):
        """Session limits are per-user, not shared."""
        mgr = SessionManager(redis_client=None, max_sessions_per_user=2, max_concurrent_sessions=20)
        mgr.create_session("user1", "sid1", "s3://bucket/f1.h5ad")
        mgr.create_session("user1", "sid2", "s3://bucket/f2.h5ad")
        # user2 should still be able to create sessions
        s = mgr.create_session("user2", "sid3", "s3://bucket/f3.h5ad")
        assert s is not None
        assert mgr._get_user_session_count("user2") == 1


class TestSessionStatus:
    """Test session lifecycle and status transitions."""

    def test_session_created_as_active(self):
        """Newly created sessions are ACTIVE."""
        mgr = SessionManager(redis_client=None)
        s = mgr.create_session("user1", "sid1", "s3://bucket/file.h5ad")
        assert s is not None
        assert s.status == SessionStatus.ACTIVE

    def test_session_marked_completed_on_end(self):
        """Ending a session marks it COMPLETED."""
        mgr = SessionManager(redis_client=None)
        s = mgr.create_session("user1", "sid1", "s3://bucket/file.h5ad")
        mgr.end_session("sid1")
        s_ended = mgr.get_session("sid1")
        assert s_ended is not None
        assert s_ended.status == SessionStatus.COMPLETED

    def test_end_nonexistent_session_noop(self):
        """Ending a nonexistent session is safe (no-op)."""
        mgr = SessionManager(redis_client=None)
        mgr.end_session("nonexistent")
        assert mgr._get_active_count() == 0

