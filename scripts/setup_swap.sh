#!/bin/bash
# Add 4GB swap space to EC2 instance as safety buffer

set -e

echo "Setting up swap space on EC2..."

# Check if swap already exists
if sudo swapon --show | grep -q '/swapfile'; then
    echo "Swap already configured"
    sudo swapon --show
    exit 0
fi

# Create 4GB swap file
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# Make permanent
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

# Verify
echo "Swap configured successfully:"
sudo swapon --show
free -h
