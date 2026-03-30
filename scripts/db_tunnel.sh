#!/bin/bash
# Open SSH tunnel to RDS PostgreSQL via EC2.
# Forwards localhost:5433 → RDS:5432 through the EC2 instance.
#
# TablePlus / psql connection settings:
#   Host:     localhost
#   Port:     5433
#   User:     nvwa_admin
#   Database: nvwa
#   SSL mode: PREFERRED
#
# To close the tunnel later:
#   pkill -f "ssh.*5433.*nvwa-pilot-db"

set -e

SSH_KEY="/Users/zhuchen/Downloads/nvwa-key.pem"
EC2_HOST="ubuntu@3.150.203.87"
RDS_HOST="nvwa-pilot-db.cj42ock2ykh4.us-east-2.rds.amazonaws.com"
LOCAL_PORT=5433
RDS_PORT=5432

# Kill any existing tunnel on this port
if lsof -ti tcp:$LOCAL_PORT &>/dev/null; then
    echo "Closing existing tunnel on port $LOCAL_PORT..."
    lsof -ti tcp:$LOCAL_PORT | xargs kill -9
    sleep 1
fi

echo "Opening SSH tunnel: localhost:$LOCAL_PORT → $RDS_HOST:$RDS_PORT via EC2..."
ssh -i "$SSH_KEY" \
    -L ${LOCAL_PORT}:${RDS_HOST}:${RDS_PORT} \
    "$EC2_HOST" \
    -N -f -o StrictHostKeyChecking=no -o ExitOnForwardFailure=yes

echo "✅ Tunnel is open on localhost:$LOCAL_PORT"
echo "   Connect with: psql -h localhost -p $LOCAL_PORT -U nvwa_admin -d nvwa"
echo "   To close:     pkill -f 'ssh.*5433.*nvwa-pilot-db'"
