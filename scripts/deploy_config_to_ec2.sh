#!/bin/bash
# Deploy updated .env configuration to EC2

set -e

SSH_KEY="${1:-/Users/zhuchen/Downloads/nvwa-key.pem}"
EC2_HOST="ubuntu@3.150.203.87"
PROJECT_DIR="/home/ubuntu/Nvwa-Lite"

echo "Deploying configuration updates to EC2..."

# Add session config to EC2 .env if not present
ssh -i "$SSH_KEY" "$EC2_HOST" << 'EOF'
cd /home/ubuntu/Nvwa-Lite

# Backup current .env
cp .env .env.backup.$(date +%Y%m%d_%H%M%S)

# Add session config if not present
if ! grep -q "MAX_CONCURRENT_SESSIONS" .env; then
    echo "" >> .env
    echo "# Session Management Configuration" >> .env
    echo "MAX_CONCURRENT_SESSIONS=20" >> .env
    echo "MAX_SESSIONS_PER_USER=2" >> .env
    echo "SESSION_TIMEOUT_MINUTES=30" >> .env
    echo "Session config added to .env"
else
    echo "Session config already exists in .env"
fi

# Restart services to apply changes
docker-compose down
docker-compose up -d

echo "Services restarted with new configuration"
EOF

echo "✅ Configuration deployed successfully"
