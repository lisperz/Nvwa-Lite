# Deployment Strategy - Test vs Production

## Overview

To safely test Phase 1 implementation without risking production, we'll use a three-tier approach:

1. **Local Development** - Docker Compose on your machine
2. **Test/Staging Environment** - Separate AWS EC2 instance
3. **Production Environment** - Current AWS EC2 instance

## Environment Configuration

### 1. Local Development (Recommended for Initial Testing)

**Purpose:** Safe testing with no AWS costs, no production impact

**Setup:**
```bash
# Use docker-compose for complete isolation
docker-compose up --build

# Access:
# - Main app: http://localhost:8501
# - Dashboard: http://localhost:8502
# - Redis: localhost:6379
```

**Configuration (.env.local):**
```bash
# OpenAI
OPENAI_API_KEY=sk-...

# Local Redis
REDIS_HOST=localhost
REDIS_PORT=6379

# Local storage (no S3)
# S3 variables can be omitted - system will use local fallback

# Test tokens
PILOT_TOKEN_1=test_token_1
PILOT_TOKEN_2=test_token_2
```

**Advantages:**
- ✅ Zero AWS costs
- ✅ No production impact
- ✅ Fast iteration
- ✅ Easy to reset/restart

**Limitations:**
- ❌ Doesn't test S3 integration
- ❌ Doesn't test real network conditions

### 2. Test/Staging Environment (AWS EC2)

**Purpose:** Full AWS integration testing before production deployment

**Infrastructure:**
- Separate EC2 instance: `nvwa-test` or `nvwa-staging`
- Separate S3 bucket: `nvwa-test-data`
- Can share Redis or use separate instance
- Different security group (optional)

**Setup:**

1. **Create Test S3 Bucket:**
```bash
# Use different bucket name
./scripts/setup_s3.sh nvwa-test-data
```

2. **Launch Test EC2:**
```bash
# Same instance type as production
# t3.medium or t3.small
# Tag: Environment=test
```

3. **Configuration (.env.test):**
```bash
# OpenAI
OPENAI_API_KEY=sk-...

# AWS S3 - TEST BUCKET
AWS_ACCESS_KEY_ID=your-key
AWS_SECRET_ACCESS_KEY=your-secret
AWS_DEFAULT_REGION=us-east-1
S3_BUCKET_NAME=nvwa-test-data  # Different from production!

# Redis - can use different database
REDIS_HOST=localhost
REDIS_PORT=6379
REDIS_DB=1  # Production uses 0, test uses 1

# Test tokens (different from production)
PILOT_TOKEN_1=test_token_abc123
PILOT_TOKEN_2=test_token_def456
```

4. **Deploy to Test:**
```bash
# SSH to test EC2
ssh -i your-key.pem ubuntu@test-ec2-ip

# Clone repo
git clone <repo-url>
cd nvwa-lite

# Copy test environment file
cp .env.test .env

# Start services
docker-compose up -d

# Access:
# http://test-ec2-ip:8501
# http://test-ec2-ip:8502
```

**Advantages:**
- ✅ Tests full AWS integration
- ✅ Tests real network conditions
- ✅ No production impact
- ✅ Can share with team for testing

**Costs:**
- EC2 t3.small: ~$15/month (can stop when not testing)
- S3: ~$1/month
- Total: ~$16/month (stop EC2 to save costs)

### 3. Production Environment (AWS EC2)

**Purpose:** Live pilot deployment for real users

**Infrastructure:**
- Current EC2 instance
- Production S3 bucket: `nvwa-prod-data` or `nvwa-mvp-pilot-data`
- Production Redis
- Production security group

**Configuration (.env.production):**
```bash
# OpenAI
OPENAI_API_KEY=sk-...

# AWS S3 - PRODUCTION BUCKET
AWS_ACCESS_KEY_ID=your-key
AWS_SECRET_ACCESS_KEY=your-secret
AWS_DEFAULT_REGION=us-east-1
S3_BUCKET_NAME=nvwa-mvp-pilot-data  # Production bucket!

# Redis
REDIS_HOST=localhost
REDIS_PORT=6379
REDIS_DB=0  # Production uses database 0

# Production tokens (generated for real pilot users)
PILOT_TOKEN_1=prod_secure_token_1
PILOT_TOKEN_2=prod_secure_token_2
# ... etc
```

**Deployment Process:**
```bash
# SSH to production EC2
ssh -i your-key.pem ubuntu@prod-ec2-ip

# Pull latest code
cd nvwa-lite
git pull origin main

# Backup current .env
cp .env .env.backup

# Update .env if needed
# nano .env

# Restart services
docker-compose down
docker-compose up -d

# Verify
docker-compose ps
docker-compose logs -f
```

## Recommended Testing Workflow

### Phase 1: Local Testing (Today)

1. **Setup local environment:**
```bash
cd nvwa-lite
cp .env.example .env.local
# Edit .env.local with test values
ln -sf .env.local .env
```

2. **Generate test tokens:**
```bash
./scripts/generate_pilot_tokens.sh
# Use these tokens for local testing
```

3. **Start services:**
```bash
docker-compose up --build
```

4. **Run test cases:**
- Follow `docs/TESTING_GUIDE.md`
- Test all 8 scenarios
- Verify dashboard metrics
- Check logs for errors

5. **Iterate and fix:**
- If issues found, fix and restart
- No impact on production
- Fast feedback loop

### Phase 2: Test Environment (Tomorrow)

1. **Create test infrastructure:**
- Launch test EC2 instance
- Create test S3 bucket
- Configure test .env

2. **Deploy to test:**
```bash
# On test EC2
git clone <repo>
docker-compose up -d
```

3. **Full integration testing:**
- Test with real S3 uploads
- Test with real network latency
- Share with 1-2 internal testers
- Monitor dashboard for 24 hours

4. **Validate:**
- All features work
- No errors in logs
- Performance acceptable
- Dashboard shows correct metrics

### Phase 3: Production Deployment (After Validation)

1. **Generate production tokens:**
```bash
./scripts/generate_pilot_tokens.sh
# Save securely for pilot users
```

2. **Deploy to production:**
```bash
# On production EC2
git pull origin main
docker-compose down
docker-compose up -d
```

3. **Monitor closely:**
- Watch dashboard for first 24 hours
- Check logs frequently
- Be ready to rollback if needed

4. **Rollback plan:**
```bash
# If issues occur
docker-compose down
git checkout <previous-commit>
docker-compose up -d
```

## Environment Variables Checklist

### Critical Differences Between Environments

| Variable | Local | Test | Production |
|----------|-------|------|------------|
| S3_BUCKET_NAME | (omit) | nvwa-test-data | nvwa-mvp-pilot-data |
| REDIS_DB | 0 | 1 | 0 |
| PILOT_TOKEN_* | test_token_* | test_token_* | prod_secure_* |
| Log Level | DEBUG | INFO | INFO |

### Security Considerations

**Never:**
- ❌ Use production tokens in test environment
- ❌ Use test tokens in production
- ❌ Commit .env files to git
- ❌ Share production credentials

**Always:**
- ✅ Keep .env files separate
- ✅ Use different S3 buckets
- ✅ Generate new tokens for each environment
- ✅ Backup production .env before changes

## Cost Optimization

### Development Phase
- Use local Docker (free)
- Stop test EC2 when not testing
- Delete test S3 data regularly

### Production Phase
- Keep production EC2 running
- Monitor S3 costs (14-day TTL helps)
- Use t3.small if t3.medium is overkill

## Monitoring

### Local Development
- Check logs: `docker-compose logs -f`
- Dashboard: http://localhost:8502

### Test Environment
- SSH and check logs
- Dashboard: http://test-ec2-ip:8502
- CloudWatch (optional)

### Production
- Dashboard: http://prod-ec2-ip:8502
- Set up CloudWatch alarms
- Daily log review
- Weekly metrics review

## Quick Reference

### Start Local Testing
```bash
docker-compose up --build
# Main: http://localhost:8501
# Dashboard: http://localhost:8502
```

### Deploy to Test
```bash
ssh test-ec2
cd nvwa-lite && git pull
docker-compose restart
```

### Deploy to Production
```bash
ssh prod-ec2
cd nvwa-lite && git pull
docker-compose down && docker-compose up -d
```

### Rollback Production
```bash
ssh prod-ec2
cd nvwa-lite
git log --oneline  # Find previous commit
git checkout <commit-hash>
docker-compose restart
```

## Next Steps

1. **Today:** Test locally with docker-compose
2. **Tomorrow:** Set up test EC2 if local tests pass
3. **Day 3:** Deploy to production after test validation
4. **Week 1:** Monitor production closely

This approach ensures zero production risk during testing phase.
