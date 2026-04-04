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

    # Run all migrations in order
    migrations = sorted(migrations_dir.glob("*.sql"))

    if not migrations:
        print("No migration files found")
        sys.exit(1)

    for migration_file in migrations:
        run_migration(migration_file)

    print("\n✓ All migrations completed successfully")

