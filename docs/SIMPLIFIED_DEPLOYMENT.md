# Simplified Deployment Strategy - Local to Production

## Overview

For a pilot with 5-10 users, we can skip the test EC2 and go directly from local testing to production. This saves cost and complexity while still ensuring safety.

## Two-Tier Approach

### 1. Local Testing (With Real AWS Integration)
- Test on your machine with Docker Compose
- Use REAL S3 bucket for testing (separate from production)
- Test all features including S3 uploads
- Zero risk to production users

### 2. Production Deployment
- Deploy to existing EC2 after local validation
- Use production S3 bucket
- Monitor closely for first 24 hours

## Setup for Local Testing with Real S3

### Step 1: Create Test S3 Bucket

```bash
# Create a test bucket (different from production)
./scripts/setup_s3.sh nvwa-test-data

# This bucket is for LOCAL testing only
# Production will use: nvwa-mvp-pilot-data
```

### Step 2: Configure Local Environment

Create `.env` with REAL AWS credentials:

```bash
# OpenAI
OPENAI_API_KEY=sk-your-real-key

# AWS S3 - TEST BUCKET (for local testing)
AWS_ACCESS_KEY_ID=your-aws-key
AWS_SECRET_ACCESS_KEY=your-aws-secret
AWS_DEFAULT_REGION=us-east-1
S3_BUCKET_NAME=nvwa-test-data  # Test bucket, not production!

# Redis (Docker Compose will provide)
REDIS_HOST=redis
REDIS_PORT=6379
REDIS_DB=0

# Test tokens (will be replaced for production)
PILOT_TOKEN_1=test_token_1
PILOT_TOKEN_2=test_token_2
PILOT_TOKEN_3=test_token_3

# Session config
MAX_CONCURRENT_SESSIONS=2
SESSION_TIMEOUT_MINUTES=30
```

### Step 3: Comprehensive Local Testing

```bash
# Start services
./scripts/start_local_test.sh

# Access
# Main: http://localhost:8501
# Dashboard: http://localhost:8502
```

**Test Checklist:**

1. **Authentication**
   - ✅ Valid token works
   - ✅ Invalid token rejected
   - ✅ Logout works

2. **S3 Integration** (REAL AWS)
   - ✅ Upload dataset
   - ✅ Check S3 bucket for file
   - ✅ Verify path: `users/<user_id>/sessions/<session_id>/uploads/`
   - ✅ File accessible

3. **Session Management**
   - ✅ Max 1 session per user enforced
   - ✅ Max 2 concurrent sessions enforced
   - ✅ Session timeout works

4. **Logging**
   - ✅ Check `logs/user_interaction.log`
   - ✅ Check `logs/tool_execution.log`
   - ✅ Check `logs/system_metrics.log`
   - ✅ All have user_id and session_id

5. **Dashboard**
   - ✅ Metrics display correctly
   - ✅ User activity tracked
   - ✅ Tool usage shown
   - ✅ Errors logged

6. **Agent Functionality**
   - ✅ "What's in this dataset?" works
   - ✅ "Show me UMAP plot" works
   - ✅ Plots display correctly
   - ✅ No crashes

### Step 4: Verify S3 Test Data

```bash
# Check what was uploaded to test bucket
aws s3 ls s3://nvwa-test-data/users/ --recursive

# Should see files in structure:
# users/<user_id>/sessions/<session_id>/uploads/<filename>
```

### Step 5: Clean Up Test Data

```bash
# After testing, clean up test S3 bucket
aws s3 rm s3://nvwa-test-data/users/ --recursive

# Or keep it for future testing
```

## Production Deployment

### Prerequisites

- ✅ All local tests passed
- ✅ S3 integration verified locally
- ✅ Dashboard showing correct metrics
- ✅ No errors in logs

### Step 1: Prepare Production Environment

**On your local machine:**

```bash
# 1. Generate production tokens
./scripts/generate_pilot_tokens.sh

# Save pilot_tokens.json securely
# You'll distribute these to pilot users

# 2. Create production S3 bucket (if not exists)
./scripts/setup_s3.sh nvwa-mvp-pilot-data
```

### Step 2: Prepare Production .env

Create `.env.production` with production values:

```bash
# OpenAI
OPENAI_API_KEY=sk-your-real-key

# AWS S3 - PRODUCTION BUCKET
AWS_ACCESS_KEY_ID=your-aws-key
AWS_SECRET_ACCESS_KEY=your-aws-secret
AWS_DEFAULT_REGION=us-east-1
S3_BUCKET_NAME=nvwa-mvp-pilot-data  # PRODUCTION bucket!

# Redis
REDIS_HOST=localhost
REDIS_PORT=6379
REDIS_DB=0

# Production tokens (from pilot_tokens.json)
PILOT_TOKEN_1=<secure-token-1>
PILOT_TOKEN_2=<secure-token-2>
PILOT_TOKEN_3=<secure-token-3>
# ... etc for all pilot users

# Session config
MAX_CONCURRENT_SESSIONS=2
SESSION_TIMEOUT_MINUTES=30
```

### Step 3: Deploy to Production EC2

**Recommended: Deploy during off-hours (evening/weekend)**

```bash
# SSH to production EC2
ssh -i your-key.pem ubuntu@your-ec2-ip

# Navigate to project
cd nvwa-lite

# Backup current state
git branch backup-$(date +%Y%m%d-%H%M%S)
cp .env .env.backup

# Pull latest code
git pull origin main

# Update .env with production values
nano .env
# Copy contents from .env.production

# Restart services
docker-compose down
docker-compose up -d

# Verify services started
docker-compose ps

# Check logs for errors
docker-compose logs -f
```

### Step 4: Verify Production Deployment

**Immediate checks (first 5 minutes):**

```bash
# On EC2, check logs
docker-compose logs nvwa-lite | tail -50
docker-compose logs dashboard | tail -50

# Check services are running
docker-compose ps
# All should show "Up"

# Check Redis
docker-compose exec redis redis-cli ping
# Should return: PONG
```

**Access and test:**

1. Open http://your-ec2-ip:8501
2. Authenticate with a production token
3. Upload a test dataset
4. Ask a simple question
5. Check dashboard at http://your-ec2-ip:8502
6. Verify metrics appear

**Check S3:**

```bash
# Verify file uploaded to production bucket
aws s3 ls s3://nvwa-mvp-pilot-data/users/ --recursive
```

### Step 5: Monitor Closely (First 24 Hours)

**Every 2 hours for first day:**
- Check dashboard for errors
- Review logs for issues
- Verify user sessions working
- Check S3 uploads succeeding

**Dashboard URL:** http://your-ec2-ip:8502

**Log files on EC2:**
```bash
ssh your-ec2
cd nvwa-lite
tail -f logs/user_interaction.log
tail -f logs/tool_execution.log
```

## Rollback Plan

If issues occur in production:

### Quick Rollback

```bash
# SSH to EC2
ssh your-ec2

cd nvwa-lite

# Stop services
docker-compose down

# Rollback to previous commit
git log --oneline -5  # Find previous commit
git checkout <previous-commit-hash>

# Restore previous .env
cp .env.backup .env

# Restart with old version
docker-compose up -d

# Verify
docker-compose ps
docker-compose logs -f
```

### If Rollback Needed

1. Notify pilot users (if any are active)
2. Execute rollback (< 2 minutes)
3. Investigate issue locally
4. Fix and retest locally
5. Redeploy when ready

## Safety Checklist

Before deploying to production:

- [ ] All local tests passed
- [ ] S3 integration tested with real AWS
- [ ] Dashboard showing correct metrics
- [ ] No errors in local logs
- [ ] Production .env prepared
- [ ] Production tokens generated
- [ ] Rollback plan understood
- [ ] Deploying during off-hours
- [ ] Ready to monitor for 24 hours

## Cost Comparison

### With Test EC2
- Test EC2: $15/month
- Production EC2: $30/month
- S3: $3/month
- **Total: $48/month**

### Without Test EC2 (This Approach)
- Production EC2: $30/month
- S3: $3/month (test + production buckets)
- **Total: $33/month**
- **Savings: $15/month ($180/year)**

## Risk Assessment

### Risks of Skipping Test EC2

**Low Risk:**
- Local testing covers 95% of scenarios
- Graceful fallbacks built in
- Small pilot (5-10 users)
- Can rollback in < 2 minutes

**Mitigation:**
- Test with real S3 locally
- Deploy during off-hours
- Monitor closely first 24 hours
- Have rollback plan ready

### When You WOULD Need Test EC2

- Large user base (50+ users)
- Mission-critical application
- Complex infrastructure
- Multiple services/dependencies
- Compliance requirements

**For 5-10 user pilot:** Test EC2 is overkill.

## Recommended Timeline

### Day 1 (Today) - Local Testing
- Set up local environment with real S3
- Run all test cases
- Verify S3 integration
- Check dashboard metrics
- Review logs

### Day 2 (Tomorrow) - Production Deployment
- Generate production tokens
- Prepare production .env
- Deploy during evening (low usage)
- Verify deployment
- Monitor for 2-3 hours

### Day 3-7 - Monitoring
- Check dashboard daily
- Review logs for issues
- Collect user feedback
- Address any bugs

## Quick Reference

### Local Testing with Real S3
```bash
# Setup
./scripts/setup_s3.sh nvwa-test-data
cp .env.local.example .env
# Edit .env with real AWS credentials
./scripts/start_local_test.sh

# Test
# Open http://localhost:8501
# Upload dataset, verify S3

# Cleanup
docker-compose down
```

### Production Deployment
```bash
# On EC2
git pull origin main
# Update .env
docker-compose down && docker-compose up -d
docker-compose logs -f
```

### Rollback
```bash
# On EC2
docker-compose down
git checkout <previous-commit>
cp .env.backup .env
docker-compose up -d
```

## Conclusion

For a 5-10 user pilot, going directly from local testing to production is:
- ✅ Cost-effective (saves $180/year)
- ✅ Simpler (fewer environments to manage)
- ✅ Safe (with proper local testing and rollback plan)
- ✅ Practical (test EC2 provides limited additional value)

The key is **thorough local testing with real S3 integration** before production deployment.
