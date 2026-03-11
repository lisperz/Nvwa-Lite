# Environment Strategy Summary

## Three-Tier Approach

### 1. Local (Your Machine) - **START HERE**
- **Purpose:** Safe testing, zero production risk
- **Cost:** Free
- **Setup:** `./scripts/start_local_test.sh`
- **Access:** http://localhost:8501
- **When:** Initial testing, development, debugging

### 2. Test/Staging (AWS EC2) - **BEFORE PRODUCTION**
- **Purpose:** Full AWS integration testing
- **Cost:** ~$16/month (can stop when not testing)
- **Setup:** Separate EC2 + S3 bucket (nvwa-test-data)
- **Access:** http://test-ec2-ip:8501
- **When:** After local tests pass, before production

### 3. Production (AWS EC2) - **LIVE PILOT**
- **Purpose:** Real pilot users
- **Cost:** ~$37/month
- **Setup:** Current EC2 + S3 bucket (nvwa-mvp-pilot-data)
- **Access:** http://prod-ec2-ip:8501
- **When:** After test environment validation

## Key Differences

| Aspect | Local | Test | Production |
|--------|-------|------|------------|
| S3 Bucket | None (local) | nvwa-test-data | nvwa-mvp-pilot-data |
| Tokens | test_token_* | test_token_* | prod_secure_* |
| Risk | Zero | Low | High |
| Cost | Free | $16/mo | $37/mo |

## Workflow

```
Local Testing (Today)
    ↓ [If passes]
Test EC2 (Tomorrow)
    ↓ [If passes]
Production EC2 (Day 3)
    ↓
Monitor & Iterate
```

## Safety Rules

❌ **NEVER:**
- Deploy untested code to production
- Use production tokens in test
- Test on production EC2

✅ **ALWAYS:**
- Test locally first
- Use separate S3 buckets
- Keep .env files separate
- Have rollback plan

## Quick Commands

### Local Testing
```bash
./scripts/start_local_test.sh
# Test at http://localhost:8501
docker-compose down
```

### Deploy to Test
```bash
ssh test-ec2
cd nvwa-lite && git pull
docker-compose up -d
```

### Deploy to Production (After Test Validation)
```bash
ssh prod-ec2
cd nvwa-lite && git pull
docker-compose down && docker-compose up -d
```

## Documentation

- **Quick Start:** `QUICKSTART.md`
- **Full Strategy:** `docs/DEPLOYMENT_STRATEGY.md`
- **Testing Guide:** `docs/TESTING_GUIDE.md`
- **Architecture:** `docs/ARCHITECTURE.md`

## Current Status

✅ Code complete and committed
⏳ Ready for local testing
⏳ Test environment not yet set up
⏳ Production deployment pending validation
