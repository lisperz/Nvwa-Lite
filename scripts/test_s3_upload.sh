#!/bin/bash
# Test script for S3 direct upload functionality

set -e

echo "Testing S3 Direct Upload Implementation"
echo "========================================"

# Check if .env exists
if [ ! -f .env ]; then
    echo "Error: .env file not found"
    exit 1
fi

# Load environment variables
source .env

# Check required environment variables
if [ -z "$S3_BUCKET_NAME" ]; then
    echo "Error: S3_BUCKET_NAME not set in .env"
    exit 1
fi

if [ -z "$AWS_ACCESS_KEY_ID" ]; then
    echo "⚠ AWS_ACCESS_KEY_ID not set (will use IAM role if on EC2)"
else
    echo "✓ AWS credentials configured"
fi

echo "✓ S3_BUCKET_NAME configured: $S3_BUCKET_NAME"

# Test Python imports
echo ""
echo "Testing Python imports..."
export S3_BUCKET_NAME AWS_REGION AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY
uv run python -c "
from src.storage.service import S3StorageService
from src.ui.s3_uploader import s3_direct_upload_widget
import requests
print('✓ All imports successful')
"

# Test S3 service initialization
echo ""
echo "Testing S3 service initialization..."
uv run python -c "
import os
from src.storage.service import S3StorageService

s3_service = S3StorageService(
    bucket_name=os.getenv('S3_BUCKET_NAME'),
    region=os.getenv('AWS_REGION', 'us-east-2')
)
print('✓ S3 service initialized successfully')
"

# Test presigned URL generation
echo ""
echo "Testing presigned URL generation..."
uv run python -c "
import os
from src.storage.service import S3StorageService

s3_service = S3StorageService(
    bucket_name=os.getenv('S3_BUCKET_NAME'),
    region=os.getenv('AWS_REGION', 'us-east-2')
)

url, key = s3_service.generate_upload_url(
    user_id='test_user',
    session_id='test_session',
    filename='test.h5ad',
    expires_in=3600
)

print(f'✓ Presigned URL generated')
print(f'  S3 Key: {key}')
print(f'  URL length: {len(url)} chars')
"

echo ""
echo "========================================"
echo "✓ All tests passed!"
echo ""
echo "To enable direct S3 upload in production:"
echo "1. Ensure S3_BUCKET_NAME is set in .env"
echo "2. Ensure AWS credentials are configured"
echo "3. The app will automatically use direct upload"
echo ""
echo "To disable direct S3 upload (use proxy):"
echo "  Pass use_direct_s3=False to file_upload_widget()"
