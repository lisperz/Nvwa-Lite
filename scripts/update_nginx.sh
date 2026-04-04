#!/bin/bash
# Push updated nginx config to EC2 and reload.
# Fixes: client_max_body_size for .h5ad uploads, WebSocket proxy, streaming timeouts.
# Usage: bash scripts/update_nginx.sh [/path/to/ssh-key]

set -e

SSH_KEY="${1:-/Users/zhuchen/Downloads/nvwa-key.pem}"
EC2_HOST="ubuntu@3.150.203.87"

echo "=== Uploading nginx config to EC2 ==="
scp -i "$SSH_KEY" nginx/landing_ssl.conf "$EC2_HOST:/tmp/nvwa_nginx.conf"

echo "=== Installing and reloading nginx on EC2 ==="
ssh -i "$SSH_KEY" "$EC2_HOST" << 'EOF'
sudo cp /tmp/nvwa_nginx.conf /etc/nginx/sites-available/nvwa.bio
sudo nginx -t && sudo systemctl reload nginx
echo "✅ nginx reloaded successfully"
EOF

echo "=== Done ==="
echo "File uploads should now work up to 2 GB."
