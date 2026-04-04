#!/bin/bash
# Check system resources and Docker container stats

set -e

echo "=== System Resources ==="
echo ""

# Memory
echo "Memory Usage:"
free -h | grep -E "Mem|Swap"
echo ""

# CPU
echo "CPU Info:"
nproc
echo ""

# Disk
echo "Disk Usage:"
df -h / | grep -v Filesystem
echo ""

echo "=== Docker Container Stats ==="
echo ""
docker stats --no-stream --format "table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}"
echo ""

echo "=== Active Sessions (Redis) ==="
echo ""
ACTIVE_COUNT=$(docker-compose exec -T redis redis-cli GET active_sessions_count 2>/dev/null || echo "0")
echo "Active sessions: ${ACTIVE_COUNT:-0}"
echo ""

echo "=== Container Health ==="
echo ""
docker-compose ps
