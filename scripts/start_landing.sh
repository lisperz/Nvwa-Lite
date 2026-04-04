#!/bin/bash
# Start all services locally including the landing page
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

mkdir -p logs data/uploads

echo "Starting Nvwa local stack..."
docker-compose up -d landing redis nvwa-lite

echo ""
echo "All services running:"
echo "   Landing page : http://localhost"
echo "   App (token)  : http://localhost:8501"
echo "   Dashboard    : http://localhost:8502"
echo ""
echo "Flow: open http://localhost → click Request Access → http://localhost:8501"
echo ""
echo "Stop with: docker-compose down"
