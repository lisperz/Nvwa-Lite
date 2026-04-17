"""Tests for cleanup server HTTP endpoint."""

import json
from unittest.mock import Mock, patch
from src.session.cleanup_server import CleanupHandler
from src.session.manager import SessionManager


def _make_handler(session_manager, auth_service, path: str, body: str) -> CleanupHandler:
    """Instantiate CleanupHandler without triggering BaseHTTPRequestHandler.__init__."""
    with patch.object(CleanupHandler, "__init__", lambda self, *a, **kw: None):
        handler = CleanupHandler()
    handler.server = Mock()
    handler.server.session_manager = session_manager
    handler.server.auth_service = auth_service
    handler.path = path
    handler.headers = {"Content-Length": str(len(body))}
    handler.rfile = Mock()
    handler.rfile.read.return_value = body.encode()
    handler.send_response = Mock()
    handler.send_error = Mock()
    handler.end_headers = Mock()
    return handler


class TestCleanupEndpoint:
    """Test cleanup endpoint validates and ends sessions."""

    def test_valid_cleanup_request(self):
        """Valid cleanup request should end the session."""
        mgr = SessionManager(redis_client=None)
        mgr.create_session("user1", "sid1", "s3://bucket/file.h5ad")
        assert mgr._get_active_count() == 1

        auth_service = Mock()
        user_mock = Mock()
        user_mock.user_id = "user1"
        auth_service.validate_token.return_value = user_mock

        payload = json.dumps({"session_id": "sid1", "token": "valid_token"})
        handler = _make_handler(mgr, auth_service, "/cleanup", payload)
        handler.do_POST()

        handler.send_response.assert_called_once_with(200)
        assert mgr._get_active_count() == 0

    def test_invalid_token(self):
        """Invalid token should return 401."""
        mgr = SessionManager(redis_client=None)
        mgr.create_session("user1", "sid1", "s3://bucket/file.h5ad")

        auth_service = Mock()
        auth_service.validate_token.return_value = None

        payload = json.dumps({"session_id": "sid1", "token": "bad_token"})
        handler = _make_handler(mgr, auth_service, "/cleanup", payload)
        handler.do_POST()

        handler.send_error.assert_called_once()
        assert handler.send_error.call_args[0][0] == 401
        # Session should still be active
        assert mgr._get_active_count() == 1

    def test_missing_session_id(self):
        """Missing session_id should return 400."""
        auth_service = Mock()
        payload = json.dumps({"token": "valid_token"})
        handler = _make_handler(SessionManager(redis_client=None), auth_service, "/cleanup", payload)
        handler.do_POST()

        handler.send_error.assert_called_once()
        assert handler.send_error.call_args[0][0] == 400

    def test_wrong_path(self):
        """Request to wrong path should return 404."""
        handler = _make_handler(SessionManager(redis_client=None), Mock(), "/other", "{}")
        handler.do_POST()

        handler.send_error.assert_called_once()
        assert handler.send_error.call_args[0][0] == 404

    def test_idempotent_cleanup(self):
        """Calling cleanup twice should return 200 both times (idempotent)."""
        mgr = SessionManager(redis_client=None)
        mgr.create_session("user1", "sid1", "s3://bucket/file.h5ad")

        auth_service = Mock()
        user_mock = Mock()
        user_mock.user_id = "user1"
        auth_service.validate_token.return_value = user_mock

        payload = json.dumps({"session_id": "sid1", "token": "valid_token"})

        # First call
        h1 = _make_handler(mgr, auth_service, "/cleanup", payload)
        h1.do_POST()
        h1.send_response.assert_called_once_with(200)

        # Second call (session already ended)
        h2 = _make_handler(mgr, auth_service, "/cleanup", payload)
        h2.do_POST()
        h2.send_response.assert_called_once_with(200)

        assert mgr._get_active_count() == 0
