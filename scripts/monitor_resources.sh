#!/bin/bash
# Monitor resources and alert if thresholds exceeded

set -e

# Thresholds
MEM_THRESHOLD=85
CPU_THRESHOLD=90

# Get current usage
MEM_USAGE=$(free | grep Mem | awk '{printf "%.0f", $3/$2 * 100}')
CPU_USAGE=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1 | cut -d'.' -f1)

echo "=== Resource Monitor ==="
echo "Memory: ${MEM_USAGE}% (threshold: ${MEM_THRESHOLD}%)"
echo "CPU: ${CPU_USAGE}% (threshold: ${CPU_THRESHOLD}%)"
echo ""

# Check memory
if [ "$MEM_USAGE" -gt "$MEM_THRESHOLD" ]; then
    echo "⚠️  WARNING: Memory usage above ${MEM_THRESHOLD}%"
    echo "Consider reducing MAX_CONCURRENT_SESSIONS or upgrading instance"
fi

# Check CPU
if [ "$CPU_USAGE" -gt "$CPU_THRESHOLD" ]; then
    echo "⚠️  WARNING: CPU usage above ${CPU_THRESHOLD}%"
fi

# Show active sessions
ACTIVE_SESSIONS=$(docker-compose exec -T redis redis-cli GET active_sessions_count 2>/dev/null || echo "0")
echo "Active sessions: ${ACTIVE_SESSIONS:-0}/20"
