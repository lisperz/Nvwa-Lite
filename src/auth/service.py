"""Token-based authentication service for pilot users.

Provides simple token validation for the MVP pilot phase (5-10 users).
Tokens are pre-generated and distributed to pilot users via email.
"""

from __future__ import annotations

import os
import secrets
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


@dataclass
class User:
    """Represents an authenticated pilot user."""

    user_id: str
    token: str
    email: str
    created_at: datetime
    is_active: bool = True


class AuthService:
    """Token-based authentication for pilot users.

    For the MVP pilot, we use pre-generated tokens stored in environment
    variables or a simple JSON file. This is sufficient for 5-10 users.

    For production, this would be replaced with JWT tokens, OAuth, or
    a proper identity provider.
    """

    def __init__(self, tokens_file: Path | None = None):
        """Initialize the auth service.

        Args:
            tokens_file: Optional path to JSON file with pilot tokens.
                        If not provided, loads from environment variables.
        """
        self.tokens: dict[str, User] = {}
        self._load_pilot_tokens(tokens_file)

    def _load_pilot_tokens(self, tokens_file: Path | None = None) -> None:
        """Load pilot user tokens from file or environment."""
        if tokens_file and tokens_file.exists():
            # Load from JSON file
            import json
            with open(tokens_file) as f:
                data = json.load(f)
                for user_id, info in data.items():
                    self.tokens[info["token"]] = User(
                        user_id=user_id,
                        token=info["token"],
                        email=info["email"],
                        created_at=datetime.fromisoformat(info.get("created_at", datetime.utcnow().isoformat())),
                        is_active=info.get("is_active", True)
                    )
        else:
            # Load from environment variables (format: PILOT_TOKEN_001=user_id:email:token)
            for key, value in os.environ.items():
                if key.startswith("PILOT_TOKEN_"):
                    try:
                        user_id, email, token = value.split(":")
                        self.tokens[token] = User(
                            user_id=user_id,
                            token=token,
                            email=email,
                            created_at=datetime.utcnow(),
                            is_active=True
                        )
                    except ValueError:
                        continue

    def validate_token(self, token: str) -> User | None:
        """Validate a token and return the associated user.

        Args:
            token: The authentication token to validate.

        Returns:
            User object if token is valid and active, None otherwise.
        """
        user = self.tokens.get(token)
        if user and user.is_active:
            return user
        return None

    def generate_pilot_token(self, user_id: str, email: str) -> str:
        """Generate a secure token for a pilot user.

        This is used during pilot setup to create tokens for users.

        Args:
            user_id: Unique identifier for the user (e.g., "beta_user_001").
            email: User's email address.

        Returns:
            The generated token string.
        """
        token = secrets.token_urlsafe(32)
        self.tokens[token] = User(
            user_id=user_id,
            token=token,
            email=email,
            created_at=datetime.utcnow(),
            is_active=True
        )
        return token

    def revoke_token(self, token: str) -> bool:
        """Revoke a user's token.

        Args:
            token: The token to revoke.

        Returns:
            True if token was found and revoked, False otherwise.
        """
        user = self.tokens.get(token)
        if user:
            user.is_active = False
            return True
        return False

    def get_all_users(self) -> list[User]:
        """Get all registered pilot users.

        Returns:
            List of all User objects.
        """
        return list(self.tokens.values())
