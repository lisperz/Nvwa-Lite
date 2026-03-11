# Nvwa MVP Phase 1 Implementation Guide

## Overview

This guide walks you through implementing Phase 1 of the Nvwa MVP transformation, adding authentication, structured logging, S3 storage, session management, and monitoring to the existing Streamlit application.

## Prerequisites

- Existing Nvwa-Lite running on AWS EC2
- AWS account with S3 access
- Redis installed (or will use in-memory fallback)
- Python 3.11+
- uv package manager

## Phase 1 Architecture

```
┌─────────────────────────────────────────────────────────┐
│                  Streamlit App (Port 8501)               │
│                                                          │
│  ┌──────────┐  ┌──────────┐  ┌───────────────┐        │
│  │ Auth     │  │ Session  │  │ Analysis       │        │
│  │ Service  │  │ Manager  │  │ Engine         │        │
│  └────┬─────┘  └────┬─────┘  └────────┬───────┘        │
│       │             │                  │                 │
│  ┌────▼─────────────▼──────────────────▼──────────┐    │
│  │           Event Logger (structured JSON)        │    │
│  └────────────────────┬────────────────────────────┘    │
└───────────────────────┼─────────────────────────────────┘
                        │
           ┌────────────▼────────────┐
           │  Storage (S3 + local)   │
           └─────────────────────────┘
```

## Implementation Steps

### Step 1: Install Dependencies

```bash
cd /Users/zhuchen/Downloads/Nvwa\ Bio\ Technical\ Challange/nvwa-lite

# Install new dependencies
uv sync

# Verify installation
uv run python -c "import redis, boto3; print('✅ Dependencies installed')"
```

### Step 2: Setup AWS S3

```bash
# Set environment variables
export AWS_ACCESS_KEY_ID=your_access_key
export AWS_SECRET_ACCESS_KEY=your_secret_key
export AWS_REGION=us-west-2
export S3_BUCKET_NAME=nvwa-mvp-pilot

# Run S3 setup script
./scripts/setup_s3.sh
```

This creates an S3 bucket with:
- 14-day TTL for automatic cleanup
- Versioning enabled
- CORS configured for browser uploads

### Step 3: Generate Pilot Tokens

```bash
# Generate tokens for 10 pilot users
./scripts/generate_pilot_tokens.sh
```

This creates:
- `pilot_tokens.json` - Full token details with access URLs
- `pilot_tokens.env` - Environment variable format

### Step 4: Update .env File

```bash
# Copy example and edit
cp .env.example .env
nano .env
```

Add:
1. Your OpenAI API key
2. AWS credentials
3. S3 bucket name
4. Pilot user tokens from `pilot_tokens.env`

### Step 5: Setup Redis (Optional)

```bash
# Install Redis on EC2
sudo apt-get update
sudo apt-get install redis-server -y

# Configure Redis
sudo nano /etc/redis/redis.conf
# Set: maxmemory 256mb
# Set: maxmemory-policy allkeys-lru

# Start Redis
sudo systemctl start redis
sudo systemctl enable redis

# Test
redis-cli ping  # Should return PONG
```

**Note:** If Redis is not available, the system will use in-memory fallback automatically.

### Step 6: Update Docker Configuration

The new `docker-compose.yml` includes:
- Redis service
- Dashboard service (port 8502)
- Updated environment variables

```bash
# Rebuild containers
docker-compose down
docker-compose build
docker-compose up -d
```

### Step 7: Verify Installation

```bash
# Check services are running
docker-compose ps

# Check logs
docker-compose logs -f nvwa-lite

# Test authentication
curl "http://3.150.203.87/?token=YOUR_TOKEN_HERE"
```

## Next Steps

After completing the basic setup, you need to:

1. **Integrate authentication into Streamlit UI** (Task #8)
2. **Integrate logging into agent core** (Task #1)
3. **Create monitoring dashboard** (Task #3)
4. **Test end-to-end** (Task #10)

See `IMPLEMENTATION_TASKS.md` for detailed task breakdown.

## Monitoring

### Log Files

Structured JSON logs are written to `logs/`:
- `tool_execution.log` - Tool calls and results
- `user_interaction.log` - User messages and responses
- `system_metrics.log` - Token usage and performance

### CEO Dashboard

Access at: `http://3.150.203.87:8502`

Shows:
- Active users (30min window)
- Token usage (24h)
- Estimated cost
- Tool success rates
- Session depth metrics

## Troubleshooting

### Redis Connection Failed

If Redis is not available, the system automatically falls back to in-memory storage. This is fine for the pilot but won't persist across restarts.

### S3 Upload Failed

Check:
1. AWS credentials are correct in `.env`
2. S3 bucket exists: `aws s3 ls s3://nvwa-mvp-pilot`
3. CORS is configured: `aws s3api get-bucket-cors --bucket nvwa-mvp-pilot`

### Authentication Not Working

Check:
1. Tokens are in `.env` file
2. Token format is correct: `user_id:email:token`
3. Streamlit app restarted after adding tokens

## Cost Estimates

### AWS Costs (Monthly)

- **EC2 t3.xlarge**: ~$120 (on-demand) or ~$75 (1-year reserved)
- **S3 Storage**: ~$5-10 (100GB with 14-day TTL)
- **Data Transfer**: ~$5-10 (minimal for pilot)
- **Total**: ~$130-145/month (on-demand) or ~$90-100/month (reserved)

### OpenAI Costs

- **gpt-4o-mini**: $0.15/1M input tokens, $0.60/1M output tokens
- **Estimated**: ~$20-50/month for 5-10 active users

## Security Notes

1. **Tokens are pre-shared secrets** - Treat them like passwords
2. **S3 bucket is private** - Only accessible via pre-signed URLs
3. **No public endpoints** - All access requires valid token
4. **Logs contain user data** - Ensure proper access controls
5. **14-day TTL** - Data automatically deleted after 2 weeks

## Support

For issues or questions:
1. Check logs: `docker-compose logs nvwa-lite`
2. Review task list: `IMPLEMENTATION_TASKS.md`
3. Consult architecture doc: `ARCHITECTURE.md`

## Next Phase

Phase 2 (when needed):
- Extract FastAPI backend
- Add job queue with Docker isolation
- Migrate to PostgreSQL
- Add WebSocket for real-time updates
- Scale to 50+ users

Current Phase 1 design supports clean migration to Phase 2 without rewrites.
