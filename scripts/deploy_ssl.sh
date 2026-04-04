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
scp -i "$SSH_KEY" nginx/landing_ssl.conf "$EC2_HOST:/tmp/"

ssh -i "$SSH_KEY" "$EC2_HOST" << 'EOF'
# Stop Docker nginx
cd /home/ubuntu/Nvwa-Lite
docker-compose stop landing

# Install nginx on host
sudo apt-get update
sudo apt-get install -y nginx

# Copy SSL config
sudo cp /tmp/landing_ssl.conf /etc/nginx/sites-available/nvwa.bio
sudo ln -sf /etc/nginx/sites-available/nvwa.bio /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default

# Copy landing page files
sudo mkdir -p /var/www/nvwa.bio
sudo cp -r /home/ubuntu/Nvwa-Lite/landing/* /var/www/nvwa.bio/

# Update nginx config to use correct path
sudo sed -i 's|/usr/share/nginx/html|/var/www/nvwa.bio|g' /etc/nginx/sites-available/nvwa.bio

# Test and reload nginx
sudo nginx -t
sudo systemctl enable nginx
sudo systemctl restart nginx
EOF

echo "✅ SSL setup complete!"
