# Quick Start - Local Testing with Real AWS

## Prerequisites
- Docker Desktop installed and running
- OpenAI API key
- AWS credentials (for S3 testing)

## Step 1: Create Test S3 Bucket (2 minutes)

```bash
cd nvwa-lite

# Create test bucket (separate from production)
./scripts/setup_s3.sh nvwa-test-data
```

## Step 2: Setup Environment (2 minutes)

```bash
# Copy environment template
cp .env.local.example .env

# Edit .env and add your credentials
nano .env  # or use your preferred editor

# Required changes:
# 1. OPENAI_API_KEY=sk-your-real-key
# 2. AWS_ACCESS_KEY_ID=your-key
# 3. AWS_SECRET_ACCESS_KEY=your-secret
# 4. S3_BUCKET_NAME=nvwa-test-data
```

## Step 3: Start Services (3 minutes)

```bash
# Quick start script handles everything
./scripts/start_local_test.sh
```

This will:
- Generate test tokens
- Start Redis, main app, and dashboard
- Create logs directory

**Access:**
- Main App: http://localhost:8501
- Dashboard: http://localhost:8502

## Step 3: Quick Smoke Test (5 minutes)

### Test Authentication
1. Open http://localhost:8501
2. Copy a token from `pilot_tokens.json`
3. Paste token and authenticate
4. Should see "Welcome, pilot_user_X@example.com"

### Test Basic Functionality
1. Upload a test .h5ad file
2. Ask: "What's in this dataset?"
3. Verify response appears
4. Check dashboard at http://localhost:8502
5. Should see 1 active user, 1 session

### Test S3 Integration (IMPORTANT!)
```bash
# Check that file was uploaded to S3
aws s3 ls s3://nvwa-test-data/users/ --recursive

# Should see: users/<user_id>/sessions/<session_id>/uploads/<filename>
```

This confirms real S3 integration works!

### Test Logging
```bash
# Check logs are being created
ls -lh logs/

# Should see:
# - user_interaction.log
# - tool_execution.log
# - system_metrics.log

# View a log file
tail -f logs/user_interaction.log
```

## Step 4: Test Concurrency (Optional, 5 minutes)

1. Keep first session open
2. Open incognito window
3. Try to authenticate as same user
4. Upload dataset → Should fail with "Unable to create session"
5. This confirms max 1 session per user works

## Step 5: Stop Services

```bash
docker-compose down
```

## If Everything Works

✅ You're ready to deploy to production EC2
✅ Follow `docs/SIMPLIFIED_DEPLOYMENT.md` for production deployment
✅ No need for separate test EC2 - local testing with real S3 is sufficient

## If Issues Occur

### Check Docker logs
```bash
docker-compose logs -f
```

### Check specific service
```bash
docker-compose logs nvwa-lite
docker-compose logs dashboard
docker-compose logs redis
```

### Restart services
```bash
docker-compose restart
```

### Full reset
```bash
docker-compose down -v  # Remove volumes
rm -rf logs/*           # Clear logs
./scripts/start_local_test.sh
```

## Common Issues

**"Redis connection failed"**
- Check: `docker-compose ps` - Redis should be running
- System will use in-memory fallback (still works)

**"OpenAI API key not configured"**
- Check .env file has: `OPENAI_API_KEY=sk-...`
- Restart: `docker-compose restart`

**"Port already in use"**
- Stop other services using 8501/8502
- Or change ports in docker-compose.yml

## Next Steps After Local Testing

1. **If tests pass:** Set up test EC2 environment
2. **If issues found:** Fix and retest locally
3. **Before production:** Test on separate EC2 first

See `docs/DEPLOYMENT_STRATEGY.md` for full deployment workflow.
