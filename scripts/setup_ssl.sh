#!/bin/bash
# Install SSL certificate using Certbot and configure nginx

set -e

echo "=== Installing Certbot ==="
sudo apt-get update
sudo apt-get install -y certbot python3-certbot-nginx

echo "=== Obtaining SSL Certificate ==="
# Stop nginx temporarily to allow certbot standalone mode
sudo systemctl stop nginx

# Get certificate for nvwa.bio
sudo certbot certonly --standalone -d nvwa.bio -d www.nvwa.bio \
  --non-interactive --agree-tos --email admin@nvwa.bio

echo "=== Certificate obtained successfully ==="
sudo ls -la /etc/letsencrypt/live/nvwa.bio/

echo "=== Starting nginx ==="
sudo systemctl start nginx
