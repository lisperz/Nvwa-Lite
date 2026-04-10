#!/usr/bin/env bash
# Run Layer 1 unit tests against the Temple PI dataset.
#
# Usage:
#   bash scripts/run_layer1_tests.sh
#   bash scripts/run_layer1_tests.sh --tb short     # shorter tracebacks
#   bash scripts/run_layer1_tests.sh -k TestPlotting # run one section only
#
# Output:
#   Console: live pass/fail per test
#   File:    local/reports/unit/layer1_YYYYMMDD_HHMMSS.txt
#   Record:  copy to test-records/ manually when satisfied
#
# To also save plot PNGs to test-records/ for visual review:
#   SAVE_TEST_PLOTS=1 bash scripts/run_layer1_tests.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPORT_DIR="$REPO_ROOT/local/reports/unit"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
REPORT_FILE="$REPORT_DIR/layer1_${TIMESTAMP}.txt"
DATASET="${TEMPLE_DATASET_PATH:-/Users/zhuchen/Downloads/GSE223414_slim.h5ad}"

mkdir -p "$REPORT_DIR"

echo "========================================================"
echo "  Nvwa-Lite — Layer 1 Unit Tests"
echo "  Dataset : $DATASET"
echo "  Report  : $REPORT_FILE"
echo "========================================================"
echo ""

cd "$REPO_ROOT"

# Pass any extra args (e.g. -k, --tb) straight through to pytest
PYTHONPATH=. uv run pytest tests/unit/test_layer1_tools.py \
    -v \
    --tb=short \
    --durations=10 \
    "$@" \
    2>&1 | tee "$REPORT_FILE"

echo ""
echo "Report saved → $REPORT_FILE"
echo "To promote to test-records/:"
echo "  cp \"$REPORT_FILE\" test-records/$(date +%Y-%m-%d)_layer1-baseline.txt"
