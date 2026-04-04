#!/bin/bash
# Deploy pilot_tokens.json to EC2 instance
# Usage: ./scripts/deploy-tokens-to-ec2.sh [path-to-ssh-key]

set -e

EC2_IP="3.150.203.87"
EC2_USER="ubuntu"
EC2_PATH="~/Nvwa-Lite"
LOCAL_TOKENS="pilot_tokens.json"

# Check if SSH key is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <path-to-ssh-key>"
    echo ""
    echo "Example:"
    echo "  $0 ~/.ssh/my-ec2-key.pem"
    echo ""
    echo "If you don't have the key path, try:"
    echo "  ls ~/.ssh/*.pem"
    echo "  ls ~/Downloads/*.pem"
    exit 1
fi

SSH_KEY="$1"

# Verify SSH key exists
if [ ! -f "$SSH_KEY" ]; then
    echo "Error: SSH key not found at: $SSH_KEY"
    exit 1
fi

# Verify tokens file exists
if [ ! -f "$LOCAL_TOKENS" ]; then
    echo "Error: pilot_tokens.json not found in current directory"
    echo "Please run this script from the nvwa-lite directory"
    exit 1
fi

echo "=========================================="
echo "Deploying pilot_tokens.json to EC2"
echo "=========================================="
echo "EC2 IP: $EC2_IP"
echo "SSH Key: $SSH_KEY"
echo ""

# Copy tokens file to EC2
echo "Step 1: Copying pilot_tokens.json to EC2..."
scp -i "$SSH_KEY" "$LOCAL_TOKENS" "${EC2_USER}@${EC2_IP}:${EC2_PATH}/"

if [ $? -eq 0 ]; then
    echo "✓ File copied successfully"
else
    echo "✗ Failed to copy file"
    exit 1
fi

# Restart the nvwa-lite container
echo ""
echo "Step 2: Restarting nvwa-lite container..."
ssh -i "$SSH_KEY" "${EC2_USER}@${EC2_IP}" "cd ${EC2_PATH} && docker-compose restart nvwa-lite"

if [ $? -eq 0 ]; then
    echo "✓ Container restarted successfully"
else
    echo "✗ Failed to restart container"
    exit 1
fi

echo ""
echo "=========================================="
echo "Deployment completed successfully!"
echo "=========================================="
echo ""
echo "You can now test authentication at:"
echo "http://${EC2_IP}:8501"
echo ""
echo "Try logging in with one of these tokens:"
grep -o '"token": "[^"]*"' "$LOCAL_TOKENS" | head -3 | sed 's/"token": /  - /'
