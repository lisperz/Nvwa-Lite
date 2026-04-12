#!/bin/bash
# One-time CORS configuration for browser-direct S3 uploads
# Run this once before using the direct upload feature

set -e

source .env

if [ -z "$S3_BUCKET_NAME" ]; then
    echo "Error: S3_BUCKET_NAME not set in .env"
    exit 1
fi

export S3_BUCKET_NAME AWS_REGION AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY

echo "Configuring CORS on bucket: $S3_BUCKET_NAME"

uv run python -c "
import os
from src.storage.service import S3StorageService

s3 = S3StorageService(
    bucket_name=os.getenv('S3_BUCKET_NAME'),
    region=os.getenv('AWS_REGION', 'us-east-2'),
)
s3.configure_cors([
    'https://nvwa.bio',
    'http://localhost:8501',
    'http://localhost:8502',
])
print('✓ CORS configured successfully')
"
