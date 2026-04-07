#!/bin/bash
# Run all unit tests and save a timestamped report.
# Usage: bash scripts/run_unit_tests.sh

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT="local/reports/unit/unit_test_report_${TIMESTAMP}.md"

mkdir -p local/reports/unit

uv run python -m pytest tests/unit/ -v 2>&1 | tee "$REPORT"

echo ""
echo "Report saved: $REPORT"
