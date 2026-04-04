#!/bin/bash
# Run database migration to add S3 support

set -e

cd "$(dirname "$0")/.."

echo "Running database migration..."
uv run python scripts/run_migration.py

echo "✓ Database migration completed"
