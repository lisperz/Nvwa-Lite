#!/bin/bash
# Deployment script for EC2 instance
# Run this script on your EC2 server to update the application

set -e  # Exit on error

echo "🚀 Starting deployment..."

# Navigate to project directory
cd /home/ubuntu/Nvwa-Lite

# Pull latest code
echo "📥 Pulling latest code from GitHub..."
git pull origin main

# Verify the change
echo "✅ Verifying MAX_UPLOAD_MB setting..."
grep "MAX_UPLOAD_MB" src/ui/components.py

# Check .env file
if [ ! -f .env ]; then
    echo "⚠️  .env file not found!"
    echo "Please create .env file with your OPENAI_API_KEY"
    echo "Example: echo 'OPENAI_API_KEY=your-key-here' > .env"
    exit 1
else
    echo "✅ .env file exists"
    # Verify it has the API key
    if grep -q "OPENAI_API_KEY" .env; then
        echo "✅ OPENAI_API_KEY found in .env"
    else
        echo "⚠️  OPENAI_API_KEY not found in .env file!"
        exit 1
    fi
fi

# Rebuild Docker container
echo "🔨 Rebuilding Docker container..."
docker-compose down
docker-compose build --no-cache
docker-compose up -d

# Wait for container to start
echo "⏳ Waiting for container to start..."
sleep 5

# Check container status
echo "📊 Container status:"
docker-compose ps

# Show logs
echo "📝 Recent logs:"
docker-compose logs --tail=50

echo ""
echo "✅ Deployment complete!"
echo "🌐 Access your app at: http://3.150.203.87"
echo ""
echo "To view live logs, run: docker-compose logs -f"
