# Session Summary — 2026-03-29 (Latest Update)

## Recent Changes This Session (2026-03-29)

### 1. Feedback Widget - Auto-Popup & Local Testing ✅ COMPLETED

**Context:** Completed feedback widget implementation with auto-popup after 10 seconds and full local testing against production RDS database.

**Final Implementation:**
- **Auto-popup mechanism:** Uses `st.fragment(run_every=2)` to check timer every 2 seconds
- **Timer logic:** Tracks plot generation timestamp, shows feedback after 10 seconds
- **Manual trigger:** Feedback icon button in bottom-left corner (80px from bottom)
- **Skip behavior:** Closes dialog and resets timer for another 10-second cycle
- **Database:** Connected to production RDS via SSH tunnel for local testing

**Files modified:**
- `src/ui/feedback_trigger.py` — Rewrote with fragment-based timer (Option 1 approach)
- `src/ui/feedback_dialog.py` — Fixed Skip button to reset timer properly
- `src/ui/app.py` — Integrated fragment timer, fixed icon positioning, added DATABASE_URL
- `src/session/manager.py` — Added orphaned session cleanup in `_get_user_session_count()`
- `.env` — Added DATABASE_URL with URL-encoded password, increased MAX_CONCURRENT_SESSIONS to 20, added ADMIN_PASSWORD

**Database Migration:**
- ✅ `feedback_responses` table created in production RDS via TablePlus
- ✅ Indexes created: `idx_feedback_session`, `idx_feedback_user`, `idx_feedback_created`

**Status:**
- ✅ Auto-popup works (10 seconds after plot generation)
- ✅ Skip button works (resets timer for re-showing)
- ✅ Manual feedback icon works
- ✅ Database logging works (tested locally against production RDS)
- ⏳ Production deployment pending
- ⏳ Change timer to 10 minutes for production (currently 10 seconds for testing)

---

### 2. Admin Dashboard - Feedback Tab ✅ COMPLETED

**Context:** Added feedback monitoring to admin dashboard for tracking user responses.

**Implementation:**
- New "💬 Feedback" tab in admin dashboard
- Statistics: Total responses, average score, positive percentage (4-5 stars)
- Feedback table: Shows all responses with user, session, score, time saved, and comments
- Time window filtering (uses existing sidebar filter)

**Files modified:**
- `src/monitoring/analytics.py` — Added `get_feedback_responses()` and `get_feedback_stats()`
- `src/monitoring/dashboard.py` — Added Feedback tab with stats and table
- `.env` — Added ADMIN_PASSWORD for dashboard access

**Status:**
- ✅ Feedback tab working locally
- ✅ Statistics display correctly
- ✅ Feedback table shows all responses
- ⏳ Production deployment pending

---

## Previous Session Changes (2026-03-28)

### 1. Feedback Widget Implementation ✅ LOCAL TESTING (IN PROGRESS)

**Context:** Need to collect user feedback after analysis sessions for Demo Day (May 15). Every response = one data point toward funding story.

**Implementation Approach:**
- **Initial attempt:** HTML/JS widget with iframe - had overlay issues (iframe can't break out to cover parent content)
- **Final approach:** Native Streamlit dialog modal - proper overlay behavior

**Current Implementation:**
- Feedback button appears after plot generation (💬 Give Feedback)
- Streamlit `@st.dialog` modal opens on click
- 3 questions: Q1 (star rating 1-5), Q2 (time saved 4 options), Q3 (optional text)
- Skip or Submit options
- Timer changed to 10 seconds for testing (production: 10 minutes)

**Files created:**
- `src/ui/feedback_dialog.py` — Streamlit dialog-based feedback widget
- `src/ui/feedback_widget.py` — HTML/JS widget (deprecated, kept for reference)
- `migrations/002_add_feedback_responses.sql` — database schema
- `docs/FEEDBACK_WIDGET.md` — implementation documentation

**Files modified:**
- `src/ui/app.py` — dialog integration, feedback button, session state management
- `src/agent/tools.py` — plot generation tracking callback
- `src/db/logger.py` — `log_feedback()` method
- `scripts/run_migration.py` — run all migrations in order

**Database:**
- New table: `feedback_responses` (response_id, session_id, user_id, q1_score, q2_time_saved, q3_open_text, created_at)

**Status:**
- ✅ Database schema created
- ✅ Backend logging implemented
- ✅ Dialog-based UI implemented
- ⏳ UI refinements in progress (positioning, timing)
- ⏳ Migration pending (needs DATABASE_URL on EC2)
- ⏳ Production deployment pending

**Known Issues:**
- Need to adjust feedback button positioning/timing
- Timer currently set to 10s for testing (change to 10min for production)

---

## Previous Session Changes (2026-03-27)

### 1. Concurrency Scaling ✅ DEPLOYED

**Context:** User needed to support 5-10 concurrent users with 2 sessions each. Previous limit was 2 total sessions system-wide.

**Changes:**
- Updated `SessionManager` to support `max_concurrent_sessions=20` and `max_sessions_per_user=2`
- Added environment variables: `MAX_CONCURRENT_SESSIONS`, `MAX_SESSIONS_PER_USER`, `SESSION_TIMEOUT_MINUTES`
- Updated `src/ui/app.py` to read config from environment
- Added 4GB swap space on EC2 as safety buffer

**Capacity:**
- 20 total sessions system-wide
- 2 sessions per user
- Supports 10 users at full capacity

**Files modified:**
- `src/session/manager.py` — added per-user limits, environment config
- `src/ui/app.py` — reads session config from env vars
- `.env.example` — added session management variables

**New files:**
- `docs/CONCURRENCY_SCALING.md` — scaling guide with EC2 recommendations
- `scripts/check_resources.sh` — monitor system resources
- `scripts/monitor_resources.sh` — alert on high usage

**Deployment:**
- ✅ Configuration deployed to EC2
- ✅ 4GB swap space added
- ✅ All services running with new limits

---

### 2. SSL Certificate Setup (Let's Encrypt) ✅ DEPLOYED

**Context:** Site was down with Cloudflare Error 521. Cloudflare was in "Full (strict)" mode but EC2 had no SSL certificate.

**Solution:** Installed Let's Encrypt SSL certificate on EC2 with nginx as reverse proxy.

**Implementation:**
- Obtained SSL certificate via Certbot for nvwa.bio and www.nvwa.bio
- Configured nginx on host (not Docker) to handle HTTPS on port 443
- Set up reverse proxy: nginx → Streamlit (8501) and dashboard (8502)
- Configured auto-renewal using webroot method

**SSL Setup:**
- Certificate: `/etc/letsencrypt/live/nvwa.bio/`
- Expires: June 25, 2026 (auto-renews ~May 26)
- Auto-renewal: Certbot timer runs twice daily
- Renewal test: PASSED ✅

**Nginx Configuration:**
- Port 80: Redirects to HTTPS
- Port 443: SSL termination + reverse proxy
- Locations: `/` (landing), `/app` (Streamlit), `/admin` (dashboard)

**Files created:**
- `nginx/landing_ssl.conf` — nginx SSL configuration
- `scripts/setup_ssl.sh` — SSL installation script
- `scripts/deploy_ssl.sh` — deployment script
- `docs/FIX_SSL_MODE.md` — SSL troubleshooting guide
- `docs/CLOUDFLARE_521_FIX.md` — Cloudflare error guide

**Deployment:**
- ✅ SSL certificate installed and verified
- ✅ Nginx configured with HTTPS
- ✅ Site accessible at https://nvwa.bio
- ✅ Cloudflare "Full (strict)" mode working

---

### 3. Local Database Access via SSH Tunnel ✅

**Context:** User needed to view RDS PostgreSQL tables from local machine. RDS security group only allows connections from EC2 instance, not public internet.

**Solution:** Set up SSH tunnel through EC2 instance to forward local port to RDS database.

**SSH Tunnel Command:**
```bash
ssh -i "/Users/zhuchen/Downloads/nvwa-key.pem" \
  -L 5433:nvwa-pilot-db.cj42ock2ykh4.us-east-2.rds.amazonaws.com:5432 \
  ubuntu@3.150.203.87 -N -f
```

**TablePlus Connection Settings:**
- Host: `localhost` (not the RDS hostname)
- Port: `5433` (local tunnel port, not 5432)
- User: `nvwa_admin`
- Password: `NvwaBio#2026!X9qL`
- Database: `nvwa` (not `nvwa_pilot`)
- SSL mode: `PREFERRED`

**To close tunnel later:**
```bash
pkill -f "ssh.*5433.*nvwa-pilot-db"
```

**Database credentials location:** EC2 `.env` file contains `DATABASE_URL` with full connection string

---

## Previous Session Changes (2026-03-25)

### 1. S3 Integration for File Uploads and Artifacts ✅ DEPLOYED

**Context:** Wired S3 storage for user uploads and generated artifacts (plots/tables). Successfully deployed to EC2 and migration completed.

**Implementation:**
- **File uploads**: Modified `file_upload_widget()` to upload .h5ad files to S3 while keeping local copy for immediate loading
- **Artifact storage**: Updated `DatabaseLogger.log_artifacts()` to upload plots (PNG) and tables (CSV) to S3
- **Database schema**: Added `s3_key` column to `message_artifacts` table via migration
- **Fallback behavior**: System gracefully falls back to local storage if S3 is unavailable

**S3 key structure:**
```
users/{user_id}/sessions/{session_id}/
  ├── uploads/{filename}.h5ad
  └── results/
      ├── plot/plot_{message_id}_{idx}.png
      └── table/table_{message_id}_{idx}.csv
```

**Files modified:**
- `src/ui/components.py` — `file_upload_widget()` now returns `(path, s3_key)` and uploads to S3
- `src/ui/app.py` — updated to handle new upload signature and generate session_id early
- `src/db/logger.py` — `log_artifacts()` uploads plots/tables to S3 and stores s3_key

**New files:**
- `migrations/001_add_s3_key_to_artifacts.sql` — migration to add s3_key column
- `scripts/run_migration.py` — Python migration runner
- `scripts/migrate_db.sh` — shell script to run migrations
- `docs/S3_INTEGRATION.md` — S3 integration documentation

**Deployment:**
- ✅ Code pushed to GitHub (commit 215059f)
- ✅ Deployed to EC2 via `docker-compose build --no-cache && docker-compose up -d`
- ✅ Database migration executed: `s3_key` column added to `message_artifacts`
- ✅ Services running: nvwa-lite (8501), dashboard (8502), redis (healthy)
- ✅ S3 bucket: `nvwa-mvp-pilot` (us-east-2) with IAM role attached

**Note:** Local `.env` file is gitignored (correct). EC2 instance already has `AWS_REGION=us-east-2` and `S3_BUCKET_NAME=nvwa-mvp-pilot` configured.

---

### 2. Verified DB Logging Live on EC2

**Context:** Confirmed all 5 RDS tables have live data after real user interactions.

**Verified:**
- `users`: 10 pilot users seeded
- `analysis_sessions`: session rows created correctly
- `tool_executions`: all tool calls logged with duration and status
- `chat_messages`: user messages logged
- `token_usage`: token counts logged per turn

---

### 2. Full Chat Log in Admin Dashboard (User + Assistant)

**Context:** Admin could see user messages but not agent responses. Fixed two bugs and added assistant logging.

**Bugs fixed:**
1. `get_recent_sessions()` in `analytics.py` queried `started_at` — column doesn't exist, it's `created_at` → fixed throughout; this was causing "No sessions found for this user" in the Users tab
2. `dashboard.py` referenced `started_at` in the Sessions tab dataframe → fixed to `created_at`

**New feature:**
- `DatabaseLogger.log_assistant_message()` — new method, returns the inserted row `id` (via `RETURNING id`)
- `core.py` now calls `log_assistant_message()` after every agent turn (unconditionally, even if `final_text` is empty — logs `"(no text response)"` as fallback)
- `chat_messages` table now stores both `role='user'` and `role='assistant'` rows in chronological order
- Dashboard session drill-down renders the full conversation as a chat thread

**Files modified:**
- `src/db/logger.py` — added `log_assistant_message()` returning `int | None`
- `src/agent/core.py` — call `log_assistant_message()` unconditionally after each turn
- `src/monitoring/analytics.py` — fixed `started_at` → `created_at` in `get_recent_sessions()`
- `src/monitoring/dashboard.py` — fixed `started_at` → `created_at` in Sessions tab display

---

### 3. Plots, Tables, and Generated Code in Admin Dashboard

**Context:** Admin could see text chat but not the plots/tables/code the agent generated. Added a new `message_artifacts` table and full artifact logging pipeline.

**New RDS table:**
```sql
CREATE TABLE message_artifacts (
    id            SERIAL PRIMARY KEY,
    message_id    INTEGER NOT NULL REFERENCES chat_messages(id) ON DELETE CASCADE,
    session_id    TEXT NOT NULL,
    user_id       TEXT NOT NULL,
    artifact_type TEXT NOT NULL,  -- 'plot' | 'table'
    title         TEXT,
    image_b64     TEXT,           -- base64-encoded PNG for plots
    csv_data      TEXT,           -- raw CSV for tables
    display_df    TEXT,           -- markdown table preview
    code          TEXT,           -- generated Python code
    created_at    TIMESTAMPTZ DEFAULT NOW()
)
```

**New method:** `DatabaseLogger.log_artifacts(message_id, session_id, user_id, plot_results, table_results)`
- Encodes `PlotResult.image` bytes as base64 and stores in `image_b64`
- Stores `TableResult.csv_data`, `display_df`, `code` for tables

**Pipeline in `core.py`:**
1. `log_assistant_message()` → returns `msg_id`
2. `log_artifacts(msg_id, ..., get_plot_results(), get_table_results())` — called immediately after

**Dashboard rendering (foldable expanders under each assistant message):**
- 📊 Plot artifacts: shows generated code + image (`io.BytesIO` decode)
- 📋 Table artifacts: shows generated code + markdown preview + CSV download button

**Bugs fixed:**
1. Stale `.pyc` bytecode — container was running old compiled code that didn't include `log_artifacts`; fixed by clearing all `__pycache__` inside container and restarting
2. `st.image(raw_bytes)` raises `PIL.UnidentifiedImageError` — Streamlit requires `io.BytesIO(bytes)`, not raw bytes

**Files modified:**
- `src/db/logger.py` — added `log_artifacts()`
- `src/agent/core.py` — call `log_artifacts()` after `log_assistant_message()`
- `src/monitoring/analytics.py` — `get_session_messages()` now joins `message_artifacts` and attaches `artifacts` list to each message
- `src/monitoring/dashboard.py` — renders artifacts as foldable expanders with plot/table/code/download

---

## Current System State

### Working Features
✓ Three distinct DE analysis modes with correct intent disambiguation
✓ `n_genes=0` support for computing ALL markers (not just top 20)
✓ Cluster resolution supporting both numeric IDs and cell type names
✓ Scope propagation from analysis to export
✓ System prompt with original examples and behavior
✓ Docker live code mounting for development
✓ Sidebar: Upload + Pipeline + Example Queries only
✓ Dotplot/heatmap shows correct per-cluster gene count
✓ `https://nvwa.bio` → landing page with SSL
✓ `https://nvwa.bio/app` → Streamlit app (token auth)
✓ `https://nvwa.bio/admin` → Admin dashboard (password-protected)
✓ RDS PostgreSQL DB logging (dual-write: file + DB)
✓ 10 pilot users seeded in RDS `users` table
✓ Admin can view full chat logs (user + assistant) per session
✓ Admin can view plots, tables, and generated code per assistant turn (foldable)
✓ S3 integration for file uploads and artifacts
✓ Concurrency: 20 sessions, 2 per user (supports 10 concurrent users)
✓ SSL certificate with auto-renewal
✓ Nginx reverse proxy with HTTPS

### Known Issues
- None currently blocking

### Tool Count
25 tools total

### RDS Schema (6 tables)
- `users` — pilot user tokens (SHA-256 hashed)
- `analysis_sessions` — one row per dataset upload session
- `tool_executions` — every tool call with args, result preview, duration, status
- `token_usage` — LLM token counts per turn
- `chat_messages` — user + assistant messages with `role` column
- `message_artifacts` — plots (base64 PNG) and tables (CSV + markdown) linked to assistant messages

### Infrastructure
- **EC2:** `ubuntu@3.150.203.87`, project at `/home/ubuntu/Nvwa-Lite`
- **Instance:** t3.xlarge (16GB RAM, 4 vCPUs) + 4GB swap
- **SSH key:** `/Users/zhuchen/Downloads/nvwa-key.pem`
- **RDS:** `nvwa-pilot-db.cj42ock2ykh4.us-east-2.rds.amazonaws.com`
- **S3 bucket:** `nvwa-mvp-pilot` (us-east-2) — IAM role attached, wired for uploads/artifacts
- **GitHub:** `https://github.com/lisperz/Nvwa-Lite`
- **Domain:** `https://nvwa.bio` via Cloudflare (Full strict mode)
- **SSL:** Let's Encrypt certificate (auto-renews), expires June 25, 2026
- **Nginx:** Host-level reverse proxy on ports 80/443

### EC2 .env keys (summary)
- `DATABASE_URL` — full PostgreSQL DSN
- `ADMIN_PASSWORD` — dashboard login (`NvwaAdmin2026Secure`)
- `AWS_REGION=us-east-2`
- `S3_BUCKET_NAME=nvwa-mvp-pilot`
- `OPENAI_API_KEY` — agent LLM
- `MAX_CONCURRENT_SESSIONS=20`
- `MAX_SESSIONS_PER_USER=2`
- `SESSION_TIMEOUT_MINUTES=30`

## Next Session Priorities

### High Priority
1. **Change admin password** — Replace default `NvwaAdmin2026Secure` with something personal

### Medium Priority
1. **Narrow IAM S3 policy** — Replace `AmazonS3FullAccess` with scoped JSON policy on `nvwa-ec2-role`
2. **Monitor resource usage** — Track memory/CPU under load with multiple users

### Low Priority
1. **Streamlit deprecation** — Update `use_container_width` to `width` parameter
2. **Session expiry** — Add TTL to in-memory session fallback

## Key Learnings (Cumulative)

1. **Prompt Stability:** Even conservative prompt changes significantly impact agent behavior; prefer tool-level fixes
2. **Scope Propagation:** Critical to pass target identifiers through entire pipeline (resolve → analyze → export)
3. **Intent Disambiguation:** Explicit mode definitions and examples are necessary for correct tool selection
4. **Brace Escaping:** Python `.format()` requires `{{{{` to produce literal `{{` in output
5. **Agent Output Clarity:** Include explicit directives like "IMPORTANT: Use ALL genes below" to prevent LLM from truncating
6. **Docker State:** Redis `FLUSHALL` + container restart resets all login/session state
7. **n_genes=0 Pattern:** Use 0 as sentinel for "all genes", convert to `adata.n_vars` at the analysis layer
8. **`docker-compose restart` ≠ recreate:** Changing volumes requires stop + rm + up to take effect
9. **nginx home dir access:** `www-data` needs `o+x` on all parent dirs of the webroot
10. **`#` in DATABASE_URL:** Safe in `.env` file; breaks shell interpolation — never pass as CLI arg
11. **`token_hash` not `token`:** RDS `users` table stores SHA-256 hash, not raw token
12. **PostgreSQL Decimal:** `SUM()` returns `decimal.Decimal`, not `float` — cast at the DB layer in `_query()`, not at call sites
13. **`__init__.py` import drift:** When renaming/removing public classes, always update the package `__init__.py` immediately or leave it empty
14. **RDS column name drift:** `analysis_sessions` uses `created_at` not `started_at` — always verify schema before writing queries
15. **`log_assistant_message` must return row id:** Use `RETURNING id` so the caller can link artifacts to the message
16. **`st.image()` needs BytesIO:** `st.image(raw_bytes)` raises `PIL.UnidentifiedImageError`; always wrap with `io.BytesIO(bytes)`
17. **Stale `.pyc` bytecode:** Volume-mounted source changes are live immediately, but if `.pyc` files exist from a previous build they can shadow the new `.py` files — clear with `find /app/src -name '*.pyc' -delete` inside the container and restart
18. **Cloudflare SSL modes:** "Flexible" = HTTP to origin, "Full" = HTTPS to origin (any cert), "Full (strict)" = HTTPS to origin (valid cert required)
19. **Let's Encrypt renewal:** Certbot with `--standalone` requires stopping nginx; use `webroot` method for zero-downtime renewals
20. **Nginx as reverse proxy:** Host-level nginx can proxy to Docker containers on localhost ports without exposing them directly
21. **Redis purpose:** Session management for multi-user apps — tracks active sessions, enforces concurrency limits, stores session state
22. **Docker Compose version bug:** Old docker-compose (1.29.2) has `ContainerConfig` KeyError when recreating containers with volume changes — use `down` then `up` instead of restart
