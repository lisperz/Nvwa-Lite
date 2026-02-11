#!/bin/bash
set -e

echo "Stopping Nvwa-Lite..."
pkill -f "streamlit run src/ui/app.py" 2>/dev/null && echo "Stopped." || echo "Not running."
