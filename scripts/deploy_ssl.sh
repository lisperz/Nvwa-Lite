#!/bin/bash
# Deploy SSL configuration to EC2

set -e

SSH_KEY="/Users/zhuchen/Downloads/nvwa-key.pem"
EC2_HOST="ubuntu@3.150.203.87"

echo "Step 1: Uploading SSL setup script..."
scp -i "$SSH_KEY" scripts/setup_ssl.sh "$EC2_HOST:/tmp/"

echo "Step 2: Installing SSL certificate..."
ssh -i "$SSH_KEY" "$EC2_HOST" "bash /tmp/setup_ssl.sh"

echo "Step 3: Updating nginx configuration..."
scp -i "$SSH_KEY" infra/nginx/nvwa.bio.conf "$EC2_HOST:/tmp/"

ssh -i "$SSH_KEY" "$EC2_HOST" << 'EOF'
# Stop Docker landing service (if running)
cd /home/ubuntu/Nvwa-Lite
docker compose stop landing 2>/dev/null || true

# Install nginx on host
sudo apt-get update
sudo apt-get install -y nginx

# Copy SSL config
sudo cp /tmp/nvwa.bio.conf /etc/nginx/sites-available/nvwa.bio
sudo ln -sf /etc/nginx/sites-available/nvwa.bio /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default

# Copy landing page files
sudo mkdir -p /var/www/nvwa.bio
sudo cp -r /home/ubuntu/Nvwa-Lite/landing/* /var/www/nvwa.bio/

# Test and reload nginx
sudo nginx -t
sudo systemctl enable nginx
sudo systemctl restart nginx
EOF

echo "✅ SSL setup complete!"
