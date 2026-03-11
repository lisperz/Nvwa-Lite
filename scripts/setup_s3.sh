#!/bin/bash
# Setup S3 bucket for Nvwa MVP pilot

set -e

BUCKET_NAME="${S3_BUCKET_NAME:-nvwa-mvp-pilot}"
REGION="${AWS_REGION:-us-west-2}"

echo "🪣 Setting up S3 bucket: $BUCKET_NAME"

# Check if bucket exists
if aws s3 ls "s3://$BUCKET_NAME" 2>&1 | grep -q 'NoSuchBucket'; then
    echo "Creating bucket..."
    aws s3 mb "s3://$BUCKET_NAME" --region "$REGION"
else
    echo "Bucket already exists"
fi

# Create lifecycle policy file
cat > /tmp/s3_lifecycle.json << 'EOF'
{
  "Rules": [
    {
      "Id": "DeleteOldPilotData",
      "Status": "Enabled",
      "Prefix": "users/",
      "Expiration": {
        "Days": 14
      }
    }
  ]
}
EOF

# Set lifecycle policy (14-day TTL)
echo "Setting lifecycle policy (14-day TTL)..."
aws s3api put-bucket-lifecycle-configuration \
  --bucket "$BUCKET_NAME" \
  --lifecycle-configuration file:///tmp/s3_lifecycle.json

# Enable versioning (for safety)
echo "Enabling versioning..."
aws s3api put-bucket-versioning \
  --bucket "$BUCKET_NAME" \
  --versioning-configuration Status=Enabled

# Set CORS policy for direct browser uploads
cat > /tmp/s3_cors.json << 'EOF'
{
  "CORSRules": [
    {
      "AllowedOrigins": ["*"],
      "AllowedMethods": ["GET", "PUT", "POST"],
      "AllowedHeaders": ["*"],
      "MaxAgeSeconds": 3000
    }
  ]
}
EOF

echo "Setting CORS policy..."
aws s3api put-bucket-cors \
  --bucket "$BUCKET_NAME" \
  --cors-configuration file:///tmp/s3_cors.json

# Clean up temp files
rm /tmp/s3_lifecycle.json /tmp/s3_cors.json

echo ""
echo "✅ S3 bucket configured successfully!"
echo ""
echo "Bucket details:"
echo "  Name: $BUCKET_NAME"
echo "  Region: $REGION"
echo "  TTL: 14 days"
echo "  Versioning: Enabled"
echo ""
echo "Add to your .env file:"
echo "  S3_BUCKET_NAME=$BUCKET_NAME"
echo "  AWS_REGION=$REGION"
