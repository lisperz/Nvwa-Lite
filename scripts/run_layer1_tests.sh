#!/usr/bin/env bash
# Run Layer 1 unit tests against any h5ad dataset.
#
# Usage:
#   # Portable tests only (any dataset):
#   DATASET_PATH=/path/to/any.h5ad bash scripts/run_layer1_tests.sh
#
#   # Portable + Temple regression tests:
#   DATASET_PATH=/path/to/GSE223414_slim.h5ad DATASET_NAME=temple bash scripts/run_layer1_tests.sh
#
#   # Legacy default (Temple, backwards compatible):
#   bash scripts/run_layer1_tests.sh
#
#   # Shorter tracebacks:
#   DATASET_PATH=... bash scripts/run_layer1_tests.sh --tb short
#
#   # Run one section only:
#   DATASET_PATH=... bash scripts/run_layer1_tests.sh -k TestQCMetrics
#
#   # Also save plot PNGs to test-records/ for visual review:
#   SAVE_TEST_PLOTS=1 DATASET_PATH=... bash scripts/run_layer1_tests.sh
#
# Output:
#   Console: live pass/fail per test
#   File:    local/reports/unit/layer1_YYYYMMDD_HHMMSS.txt
#   Record:  copy to test-records/ manually when satisfied

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPORT_DIR="$REPO_ROOT/local/reports/unit"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

# Resolve dataset path — DATASET_PATH takes priority, then TEMPLE_DATASET_PATH (legacy), then default
DATASET="${DATASET_PATH:-${TEMPLE_DATASET_PATH:-/Users/zhuchen/Downloads/GSE223414_slim.h5ad}}"
DATASET_NAME_VAL="${DATASET_NAME:-}"

# If DATASET_NAME not set but using the default Temple path, auto-enable regression tests
if [ -z "$DATASET_NAME_VAL" ] && [ "$DATASET" = "/Users/zhuchen/Downloads/GSE223414_slim.h5ad" ]; then
    DATASET_NAME_VAL="temple"
fi

REPORT_FILE="$REPORT_DIR/layer1_${TIMESTAMP}.txt"

mkdir -p "$REPORT_DIR"

echo "========================================================"
echo "  Nvwa-Lite — Layer 1 Unit Tests"
echo "  Dataset     : $DATASET"
echo "  Dataset name: ${DATASET_NAME_VAL:-unknown (portable tests only)}"
if [ "$DATASET_NAME_VAL" = "temple" ]; then
    echo "  Regression  : Temple-specific regression tests ENABLED"
else
    echo "  Regression  : Skipped (set DATASET_NAME=temple to enable)"
fi
echo "  Report      : $REPORT_FILE"
echo "========================================================"
echo ""

cd "$REPO_ROOT"

# Choose which test files to run
if [ "$DATASET_NAME_VAL" = "temple" ]; then
    TEST_TARGETS="tests/unit/test_layer1_portable.py tests/unit/test_temple_regression.py"
else
    TEST_TARGETS="tests/unit/test_layer1_portable.py"
fi

# Pass any extra args (e.g. -k, --tb) straight through to pytest
DATASET_PATH="$DATASET" DATASET_NAME="$DATASET_NAME_VAL" \
    PYTHONPATH=. uv run pytest $TEST_TARGETS \
    -v \
    --tb=short \
    --durations=10 \
    "$@" \
    2>&1 | tee "$REPORT_FILE"

echo ""
echo "Report saved → $REPORT_FILE"
echo "To promote to test-records/:"
echo "  cp \"$REPORT_FILE\" test-records/$(date +%Y-%m-%d)_layer1-baseline.txt"
