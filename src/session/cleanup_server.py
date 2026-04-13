"""HTTP cleanup endpoint for browser sendBeacon on tab close/navigate.

Runs a lightweight HTTP server in a background thread to accept cleanup
requests from the browser when the tab closes or navigates away. This
enables immediate session cleanup without waiting for heartbeat timeout.
"""

import json
import logging
from http.server import BaseHTTPRequestHandler, HTTPServer
from threading import Thread
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.auth.service import AuthService
    from src.session.manager import SessionManager

logger = logging.getLogger(__name__)


class CleanupHandler(BaseHTTPRequestHandler):
    """Handles POST /cleanup requests from browser sendBeacon."""

    def do_POST(self) -> None:
        """Handle cleanup request."""
        if self.path != "/cleanup":
            self.send_error(404)
            return

        try:
            content_length = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(content_length)
            data = json.loads(body)

            session_id = data.get("session_id")
            token = data.get("token")

            if not session_id or not token:
                self.send_error(400, "Missing session_id or token")
                return

            # Validate token
            user = self.server.auth_service.validate_token(token)
            if not user:
                self.send_error(401, "Invalid token")
                return

            # Verify session belongs to this user
            session = self.server.session_manager.get_session(session_id)
            if session and session.user_id == user.user_id:
                self.server.session_manager.end_session(session_id)
                logger.info("Cleaned session %s via sendBeacon for user %s", session_id, user.user_id)

            # Always return 200 (idempotent — session may already be ended)
            self.send_response(200)
            self.end_headers()

        except Exception as e:
            logger.error("Cleanup endpoint error: %s", e, exc_info=True)
            self.send_error(500, str(e))

    def log_message(self, format: str, *args) -> None:
        """Suppress default HTTP server logs (use logger instead)."""
        pass


def start_cleanup_server(
    session_manager: "SessionManager",
    auth_service: "AuthService",
    port: int = 8503,
) -> None:
    """Start cleanup HTTP server in a background daemon thread.

    Args:
        session_manager: SessionManager instance.
        auth_service: AuthService instance for token validation.
        port: Port to listen on (default 8503).
    """
    server = HTTPServer(("", port), CleanupHandler)
    server.session_manager = session_manager
    server.auth_service = auth_service

    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()

    logger.info("Cleanup server started on port %d", port)
