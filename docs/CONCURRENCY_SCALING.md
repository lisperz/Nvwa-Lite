# Concurrency Scaling Guide

## Overview

Nvwa-Lite now supports **5-10 concurrent users**, each with up to **2 active sessions**.

## Configuration

Set these environment variables in `.env`:

```bash
# System-wide limit (total sessions across all users)
MAX_CONCURRENT_SESSIONS=20

# Per-user limit (sessions per individual user)
MAX_SESSIONS_PER_USER=2

# Session timeout (minutes of inactivity before auto-cleanup)
SESSION_TIMEOUT_MINUTES=30
```

## Capacity Planning

### Current Settings
- **20 total sessions** system-wide
- **2 sessions per user**
- Supports **10 users** at full capacity (10 users × 2 sessions = 20 sessions)

### Resource Requirements

**Per session:**
- Memory: ~500MB - 2GB (depends on .h5ad file size)
- CPU: Moderate during analysis, low when idle

**For 20 concurrent sessions:**
- Recommended: 16GB+ RAM, 4+ vCPUs
- Minimum: 8GB RAM, 2 vCPUs (may experience slowdowns)

## EC2 Instance Recommendations

| Users | Sessions | Recommended Instance | RAM | vCPUs |
|-------|----------|---------------------|-----|-------|
| 5     | 10       | t3.large            | 8GB | 2     |
| 10    | 20       | t3.xlarge           | 16GB| 4     |
| 15    | 30       | t3.2xlarge          | 32GB| 8     |

## Monitoring

Check current resource usage:
```bash
bash scripts/check_resources.sh
```

View active sessions in admin dashboard:
- Navigate to `https://nvwa.bio/admin`
- Check "Sessions" tab for active session count

## Scaling Up

To support more users:

1. **Increase limits** in `.env`:
   ```bash
   MAX_CONCURRENT_SESSIONS=30
   MAX_SESSIONS_PER_USER=2
   ```

2. **Upgrade EC2 instance** if needed

3. **Restart services**:
   ```bash
   docker-compose down
   docker-compose up -d
   ```

## Session Lifecycle

- **Active**: User is analyzing data
- **Idle**: No activity for < 30 minutes (still counts toward limit)
- **Expired**: Auto-cleaned after 30 minutes of inactivity
- **Completed**: User manually ended session

## Troubleshooting

**"System at capacity" error:**
- Wait for idle sessions to expire (30 min)
- Admin can manually flush Redis: `docker-compose exec redis redis-cli FLUSHALL`

**High memory usage:**
- Reduce `MAX_CONCURRENT_SESSIONS`
- Upgrade EC2 instance
- Reduce `SESSION_TIMEOUT_MINUTES` for faster cleanup
