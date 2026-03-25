#!/usr/bin/env python3
"""Run database migrations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.db.client import get_conn

def run_migration(sql_file: Path) -> None:
    """Execute a SQL migration file."""
    print(f"Running migration: {sql_file.name}")

    with get_conn() as conn:
        if conn is None:
            print("ERROR: Database connection failed")
            sys.exit(1)

        with conn.cursor() as cur:
            sql = sql_file.read_text()
            cur.execute(sql)

        print(f"✓ Migration {sql_file.name} completed successfully")

if __name__ == "__main__":
    migrations_dir = Path(__file__).parent.parent / "migrations"
    migration_file = migrations_dir / "001_add_s3_key_to_artifacts.sql"

    if not migration_file.exists():
        print(f"ERROR: Migration file not found: {migration_file}")
        sys.exit(1)

    run_migration(migration_file)
