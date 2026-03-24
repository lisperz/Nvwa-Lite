"""Seed pilot users from pilot_tokens.json into the RDS users table.

Usage (from nvwa-lite/ directory):
    DATABASE_URL="postgresql://..." uv run python scripts/seed_users.py

Or set DATABASE_URL in .env and run:
    uv run python scripts/seed_users.py
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

# Load .env if present
env_path = Path(__file__).parent.parent / ".env"
if env_path.exists():
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            key, _, val = line.partition("=")
            os.environ.setdefault(key.strip(), val.strip())

try:
    import psycopg2
except ImportError:
    print("ERROR: psycopg2 not installed. Run: uv add psycopg2-binary")
    sys.exit(1)

dsn = os.getenv("DATABASE_URL")
if not dsn:
    print("ERROR: DATABASE_URL not set")
    sys.exit(1)

tokens_path = Path(__file__).parent.parent / "pilot_tokens.json"
if not tokens_path.exists():
    print(f"ERROR: {tokens_path} not found")
    sys.exit(1)

tokens: dict = json.loads(tokens_path.read_text())

conn = psycopg2.connect(dsn)
inserted = 0
skipped = 0

try:
    with conn:
        with conn.cursor() as cur:
            for user_id, info in tokens.items():
                cur.execute(
                    """
                    INSERT INTO users (user_id, email, token)
                    VALUES (%s, %s, %s)
                    ON CONFLICT (user_id) DO NOTHING
                    """,
                    (user_id, info["email"], info["token"]),
                )
                if cur.rowcount:
                    inserted += 1
                    print(f"  Inserted: {user_id} ({info['email']})")
                else:
                    skipped += 1
                    print(f"  Skipped (exists): {user_id}")
finally:
    conn.close()

print(f"\nDone. {inserted} inserted, {skipped} skipped.")
