# Phase 1 Testing Guide

This guide walks through end-to-end testing of the Nvwa-Lite MVP Phase 1 implementation.

## Prerequisites

Before testing, ensure you have:

1. AWS credentials configured (for S3 testing)
2. Redis installed locally or Docker available
3. OpenAI API key
4. Python 3.11+ with uv installed

## Setup Steps

### 1. Generate Pilot Tokens

```bash
cd /Users/zhuchen/Downloads/Nvwa\ Bio\ Technical\ Challange/nvwa-lite
./scripts/generate_pilot_tokens.sh
```

This creates:
- `pilot_tokens.json` - Token database
- `pilot_tokens.env` - Environment variables
- Prints access URLs for distribution

### 2. Configure S3 Bucket

```bash
# Set your AWS credentials
export AWS_ACCESS_KEY_ID="your-key"
export AWS_SECRET_ACCESS_KEY="your-secret"
export AWS_DEFAULT_REGION="us-east-1"

# Run S3 setup
./scripts/setup_s3.sh nvwa-mvp-pilot-data
```

This configures:
- 14-day lifecycle policy
- Versioning enabled
- CORS for pre-signed URLs

### 3. Create .env File

```bash
cp .env.example .env
```

Edit `.env` with your values:

```bash
# OpenAI
OPENAI_API_KEY=sk-...

# AWS S3
AWS_ACCESS_KEY_ID=your-key
AWS_SECRET_ACCESS_KEY=your-secret
AWS_DEFAULT_REGION=us-east-1
S3_BUCKET_NAME=nvwa-mvp-pilot-data

# Redis (local)
REDIS_HOST=localhost
REDIS_PORT=6379

# Pilot Tokens (copy from pilot_tokens.env)
PILOT_TOKEN_1=...
PILOT_TOKEN_2=...
# ... etc
```

### 4. Start Services

#### Option A: Docker Compose (Recommended)

```bash
docker-compose up --build
```

This starts:
- Main app on http://localhost:8501
- Dashboard on http://localhost:8502
- Redis on localhost:6379

#### Option B: Local Development

Terminal 1 - Redis:
```bash
redis-server
```

Terminal 2 - Main App:
```bash
./scripts/run_dev.sh
```

Terminal 3 - Dashboard:
```bash
uv run streamlit run src/monitoring/dashboard.py --server.port=8502
```

## Test Cases

### Test 1: Authentication Flow

**Objective:** Verify token-based authentication works correctly.

**Steps:**
1. Open http://localhost:8501
2. You should see authentication form
3. Try invalid token → Should show error
4. Enter valid token from `pilot_tokens.json`
5. Should see "Welcome, [email]!" and redirect to main UI
6. Verify user email shown in sidebar
7. Click "Logout" button
8. Should return to authentication form

**Expected Results:**
- ✅ Invalid tokens rejected
- ✅ Valid tokens accepted
- ✅ User info displayed
- ✅ Logout clears session

### Test 2: Session Management

**Objective:** Verify concurrency limits and session tracking.

**Steps:**
1. Authenticate as User 1
2. Upload a dataset (e.g., pbmc3k.h5ad)
3. Verify session created in logs
4. Open incognito window
5. Authenticate as User 1 again
6. Try to upload dataset
7. Should see error: "Unable to create session"
8. Authenticate as User 2 in incognito
9. Upload dataset → Should work (User 2's first session)
10. Open another incognito window
11. Authenticate as User 3
12. Try to upload → Should fail (system at capacity: 2 concurrent)

**Expected Results:**
- ✅ Max 1 session per user enforced
- ✅ Max 2 concurrent sessions enforced
- ✅ Clear error messages

### Test 3: Structured Logging

**Objective:** Verify all events are logged with user_id and session_id.

**Steps:**
1. Authenticate and upload dataset
2. Ask: "What's in this dataset?"
3. Ask: "Show me the UMAP plot"
4. Check log files in `logs/` directory

**Expected Files:**
- `user_interaction.log` - User messages
- `tool_execution.log` - Tool calls with timing
- `system_metrics.log` - Token usage

**Verify Each Log Entry Has:**
```json
{
  "timestamp": "2026-03-11T...",
  "user_id": "pilot_user_1",
  "session_id": "uuid-here",
  "task_type": "analysis",
  "event": "...",
  "payload": {...}
}
```

**Expected Results:**
- ✅ All logs in JSON format
- ✅ user_id and session_id present
- ✅ Tool execution times recorded
- ✅ Token usage tracked

### Test 4: Tool Execution and Logging

**Objective:** Verify tools execute correctly and are logged.

**Steps:**
1. Upload pbmc3k dataset
2. Ask: "Preprocess the data"
3. Wait for completion
4. Ask: "Show me the UMAP plot"
5. Verify plot displays
6. Check `logs/tool_execution.log`

**Expected Log Entries:**
- `preprocess_data` tool call with duration
- `plot_umap` tool call with duration
- Status: "success" for both
- No errors

**Expected Results:**
- ✅ Tools execute successfully
- ✅ Timing data recorded
- ✅ No errors in logs

### Test 5: Error Handling

**Objective:** Verify errors are logged and handled gracefully.

**Steps:**
1. Upload dataset
2. Ask: "Show violin plot for NONEXISTENT_GENE"
3. Should see error message in UI
4. Check `logs/tool_execution.log`

**Expected Log Entry:**
```json
{
  "event": "tool_execution",
  "payload": {
    "tool_name": "plot_violin",
    "status": "error",
    "error_message": "Gene not found...",
    "duration_ms": 123
  }
}
```

**Expected Results:**
- ✅ Error displayed to user
- ✅ Error logged with context
- ✅ System remains stable

### Test 6: Monitoring Dashboard

**Objective:** Verify dashboard displays real-time metrics.

**Steps:**
1. Perform Tests 1-5 to generate activity
2. Open http://localhost:8502
3. Verify dashboard shows:
   - Active Users count
   - Total Sessions count
   - Tool usage statistics
   - Token usage
   - Recent errors (if any)
   - User activity table

**Expected Results:**
- ✅ Metrics update in real-time
- ✅ All KPIs display correctly
- ✅ User activity tracked
- ✅ Errors shown if present

### Test 7: S3 Upload (Optional)

**Objective:** Verify S3 integration works (if configured).

**Note:** This requires AWS credentials and S3 bucket setup.

**Steps:**
1. Authenticate and upload dataset
2. Check S3 bucket for uploaded file
3. Verify path: `users/<user_id>/sessions/<session_id>/uploads/<filename>`
4. Verify file accessible via pre-signed URL

**Expected Results:**
- ✅ File uploaded to correct S3 path
- ✅ Per-user namespace isolation
- ✅ Pre-signed URLs work

### Test 8: Session Timeout

**Objective:** Verify sessions expire after 30 minutes.

**Steps:**
1. Authenticate and create session
2. Wait 30 minutes (or modify timeout in code for testing)
3. Try to send message
4. Session should be expired

**Expected Results:**
- ✅ Session expires after timeout
- ✅ User notified to re-authenticate

## Validation Checklist

After completing all tests, verify:

- [ ] Authentication works with valid/invalid tokens
- [ ] Concurrency limits enforced (1 per user, 2 system-wide)
- [ ] All events logged with user_id and session_id
- [ ] Tool executions logged with timing
- [ ] Errors logged with full context
- [ ] Dashboard displays real-time metrics
- [ ] S3 uploads work (if configured)
- [ ] Session timeout works
- [ ] Logout clears session properly
- [ ] No crashes or unhandled exceptions

## Troubleshooting

### Redis Connection Failed

If you see "Redis connection failed" in logs:
- Check Redis is running: `redis-cli ping`
- System will fall back to in-memory session storage
- Dashboard may not show accurate metrics

### S3 Upload Failed

If S3 uploads fail:
- Verify AWS credentials in `.env`
- Check bucket exists and is accessible
- System will fall back to local storage

### Authentication Not Working

If authentication fails:
- Verify tokens in `.env` match `pilot_tokens.json`
- Check token format: `PILOT_TOKEN_N=token_value`
- Restart application after changing `.env`

### Dashboard Not Updating

If dashboard shows no data:
- Verify log files exist in `logs/` directory
- Check log files are valid JSON (one entry per line)
- Ensure you've performed some actions to generate logs

## Success Criteria

Phase 1 is validated when:

1. ✅ All 8 test cases pass
2. ✅ No unhandled exceptions in logs
3. ✅ Dashboard displays accurate metrics
4. ✅ Concurrency limits work correctly
5. ✅ All events logged with proper structure
6. ✅ System handles errors gracefully
7. ✅ Performance is acceptable (< 5s response time)

## Next Steps

After successful validation:

1. Generate production tokens for pilot users
2. Deploy to AWS EC2
3. Share access URLs with pilot users
4. Monitor dashboard for issues
5. Collect feedback for Phase 2

## Notes

- Keep `pilot_tokens.json` secure - it contains authentication credentials
- Monitor `logs/` directory size - implement log rotation if needed
- Dashboard auto-refreshes every 30 seconds
- Redis data persists in Docker volume `redis-data`
