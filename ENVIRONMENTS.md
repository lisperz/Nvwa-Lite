# Environment Strategy Summary

## Simplified Two-Tier Approach (Recommended for Pilot)

For a 5-10 user pilot, we use a simplified approach that saves cost and complexity:

### 1. Local (Your Machine) - **START HERE**
- **Purpose:** Safe testing with real AWS integration
- **Cost:** Free (only S3 storage ~$1/month)
- **Setup:** `./scripts/start_local_test.sh`
- **Access:** http://localhost:8501
- **S3 Bucket:** nvwa-test-data (for testing)
- **When:** All development and testing

### 2. Production (AWS EC2) - **AFTER LOCAL VALIDATION**
- **Purpose:** Real pilot users
- **Cost:** ~$33/month
- **Setup:** Current EC2 + production S3 bucket
- **Access:** http://prod-ec2-ip:8501
- **S3 Bucket:** nvwa-mvp-pilot-data (production)
- **When:** After thorough local testing

## Why Skip Test EC2?

**For small pilots (5-10 users):**
- ✅ Saves $180/year
- ✅ Simpler infrastructure
- ✅ Local testing with real S3 covers 95% of scenarios
- ✅ Fast rollback if issues occur
- ✅ Graceful fallbacks built in

**Test EC2 is overkill when:**
- User base is small
- Can test with real AWS locally
- Have good rollback plan
- Can monitor closely after deployment

## Key Differences

| Aspect | Local | Production |
|--------|-------|------------|
| S3 Bucket | nvwa-test-data | nvwa-mvp-pilot-data |
| Tokens | test_token_* | prod_secure_* |
| Risk | Zero | Low (with monitoring) |
| Cost | ~$1/mo (S3) | ~$33/mo |
| AWS Integration | Real S3 | Real S3 |

## Workflow

```
Local Testing with Real S3 (Today)
    ↓ [If passes]
Production EC2 (Tomorrow)
    ↓
Monitor & Iterate
```

## Safety Rules

❌ **NEVER:**
- Deploy untested code to production
- Use production tokens in local testing
- Skip S3 integration testing locally

✅ **ALWAYS:**
- Test with real S3 locally first
- Use separate S3 buckets (test vs production)
- Keep .env files separate
- Have rollback plan ready
- Monitor closely after deployment

## Quick Commands

### Local Testing with Real S3
```bash
./scripts/setup_s3.sh nvwa-test-data
# Edit .env with AWS credentials
./scripts/start_local_test.sh
# Test at http://localhost:8501
# Verify S3: aws s3 ls s3://nvwa-test-data/users/ --recursive
docker-compose down
```

### Deploy to Production
```bash
ssh prod-ec2
cd nvwa-lite && git pull
# Update .env with production values
docker-compose down && docker-compose up -d
```

### Rollback Production
```bash
ssh prod-ec2
cd nvwa-lite
docker-compose down
git checkout <previous-commit>
cp .env.backup .env
docker-compose up -d
```

## Documentation

- **Quick Start:** `QUICKSTART.md` - 5-minute setup with real S3
- **Full Guide:** `docs/SIMPLIFIED_DEPLOYMENT.md` - Complete deployment workflow
- **Testing:** `docs/TESTING_GUIDE.md` - Comprehensive test cases
- **Architecture:** `docs/ARCHITECTURE.md` - Technical details

## Current Status

✅ Code complete and committed
⏳ Ready for local testing with real S3
⏳ Production deployment after local validation

## Cost Savings

**Two-Tier (Recommended):** $33/month
**Three-Tier (With Test EC2):** $48/month
**Savings:** $180/year

For a 5-10 user pilot, the simplified approach is more practical.
