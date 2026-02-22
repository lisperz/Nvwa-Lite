#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

# Ensure required directories exist
mkdir -p data/uploads logs

# Ensure virtual environment exists
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv sync
fi

echo "Starting Nvwa-Lite..."
uv run streamlit run src/ui/app.py \
    --server.port=8501 \
    --server.address=0.0.0.0
