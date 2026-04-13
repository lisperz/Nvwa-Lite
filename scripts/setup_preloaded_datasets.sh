#!/usr/bin/env bash
# scripts/setup_preloaded_datasets.sh
#
# Copies pre-registered dataset files into data/preloaded/ with size validation.
# Run this from the nvwa-lite/ directory before starting the app.
#
# Usage:
#   bash scripts/setup_preloaded_datasets.sh
#
# To copy from a specific source instead of data/uploads/:
#   bash scripts/setup_preloaded_datasets.sh /path/to/source/dir

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEST_DIR="$PROJECT_ROOT/data/preloaded"
SOURCE_DIR="${1:-$PROJECT_ROOT/data/uploads}"

echo "=== Nvwa-Lite Preloaded Dataset Setup ==="
echo "Source dir : $SOURCE_DIR"
echo "Dest dir   : $DEST_DIR"
echo ""

mkdir -p "$DEST_DIR"

# ---------------------------------------------------------------------------
# Dataset: GSE223414_slim.h5ad (expected ~1 GB)
# ---------------------------------------------------------------------------
FILENAME="GSE223414_slim.h5ad"
MIN_BYTES=500000000  # 500 MB minimum
SOURCE="$SOURCE_DIR/$FILENAME"
DEST="$DEST_DIR/$FILENAME"

echo "--- Processing: $FILENAME ---"

if [ ! -f "$SOURCE" ]; then
  echo "ERROR: Source file not found: $SOURCE"
  echo "  Options:"
  echo "    1. Download from S3: aws s3 cp s3://nvwa-mvp-pilot/preloadDataSet/beta_user_002/$FILENAME $SOURCE"
  echo "    2. Provide a different source dir: bash $0 /path/to/source"
  exit 1
fi

# Platform-compatible file size check
if command -v stat &> /dev/null; then
  SIZE_BYTES=$(stat -f%z "$SOURCE" 2>/dev/null || stat -c%s "$SOURCE" 2>/dev/null)
else
  SIZE_BYTES=$(wc -c < "$SOURCE")
fi

SIZE_MB=$(( SIZE_BYTES / 1048576 ))
echo "  Found: $SOURCE ($SIZE_MB MB)"

if [ "$SIZE_BYTES" -lt "$MIN_BYTES" ]; then
  echo ""
  echo "ERROR: File is only $SIZE_MB MB — expected >= $(( MIN_BYTES / 1048576 )) MB."
  echo "  This looks like a test/slim copy, NOT the full dataset."
  echo ""
  echo "  Files in source dir matching GSE223414_slim*:"
  ls -lh "$SOURCE_DIR"/GSE223414_slim* 2>/dev/null || echo "  (none found)"
  echo ""
  echo "  Refusing to copy. Please verify you have the correct file."
  exit 1
fi

if [ -f "$DEST" ]; then
  DEST_BYTES=$(stat -f%z "$DEST" 2>/dev/null || stat -c%s "$DEST" 2>/dev/null)
  if [ "$DEST_BYTES" -eq "$SIZE_BYTES" ]; then
    echo "  Already in place: $DEST ($SIZE_MB MB) — skipping copy."
  else
    echo "  Replacing existing (different size): $DEST"
    cp "$SOURCE" "$DEST"
    echo "  Done: $DEST ($SIZE_MB MB)"
  fi
else
  cp "$SOURCE" "$DEST"
  echo "  Copied: $DEST ($SIZE_MB MB)"
fi

echo ""
echo "=== Setup complete ==="
echo "Contents of $DEST_DIR:"
ls -lh "$DEST_DIR"
