#!/bin/bash
# Quick start script for local testing

set -e

echo "🚀 Nvwa-Lite Local Testing Setup"
echo "================================"
echo ""

# Check if .env exists
if [ ! -f .env ]; then
    echo "📝 Creating .env file from .env.local.example..."
    cp .env.local.example .env
    echo "⚠️  Please edit .env and add your OPENAI_API_KEY"
    echo ""
    read -p "Press Enter after you've added your OpenAI API key to .env..."
fi

# Check if OpenAI key is set
if ! grep -q "OPENAI_API_KEY=sk-" .env 2>/dev/null; then
    echo "⚠️  Warning: OPENAI_API_KEY not found in .env"
    echo "Please add your OpenAI API key to .env file"
    exit 1
fi

# Generate test tokens if not exists
if [ ! -f pilot_tokens.json ]; then
    echo "🔑 Generating test pilot tokens..."
    ./scripts/generate_pilot_tokens.sh
    echo ""
fi

# Create logs directory
mkdir -p logs

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "❌ Docker is not running. Please start Docker and try again."
    exit 1
fi

echo "🐳 Starting services with Docker Compose..."
echo ""
docker-compose up --build -d

echo ""
echo "✅ Services started successfully!"
echo ""
echo "📊 Access points:"
echo "   Landing:   http://localhost"
echo "   Main App:  http://localhost:8501"
echo "   Dashboard: http://localhost:8502"
echo ""
echo "🔑 Test tokens available in: pilot_tokens.json"
echo ""
echo "📝 View logs:"
echo "   docker-compose logs -f"
echo ""
echo "🛑 Stop services:"
echo "   docker-compose down"
echo ""
echo "📖 Testing guide: docs/TESTING_GUIDE.md"
echo ""
