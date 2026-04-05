# Session Summary — 2026-04-04 (Latest Update)

## Session Summary — 2026-04-04

### 1. Pushed Production Snapshot to Shared Repo (yzhou-nvwa/nvwa-mvp) ✅

**Context:** The other developer (Yuxin Chen / chenyux3) has a shared GitHub repo at `https://github.com/yzhou-nvwa/nvwa-mvp` that is intended to become the canonical deployed repo. Chen's local production code needed to be pushed there for comparison and reconciliation.

**What was done:**
- Remote `boss` already pointed to `https://github.com/yzhou-nvwa/nvwa-mvp.git`
- Created branch `reconcile-from-production` from local `main`
- Committed previously uncommitted local production files (nginx config, scripts, docs, discuss)
- Excluded: `src/monitoring/dashboard.py.tmp` (temp), `src/ui/feedback_widget.py` (deprecated)
- Pushed to `boss/reconcile-from-production`
- Opened PR #3: https://github.com/yzhou-nvwa/nvwa-mvp/pull/3 (review only, not for direct merge)

**Files committed in that branch beyond origin/main:**
- `nginx/landing_ssl.conf` — WebSocket keepalive, upload size, SSL proxy settings
- `.env.example` — session management env vars
- `Dockerfile` — production tweaks
- `scripts/` — SSL setup, swap space, resource monitoring, EC2 deploy scripts
- `docs/` — CLOUDFLARE_521_FIX, CONCURRENCY_SCALING, FEEDBACK_WIDGET, FIX_SSL_MODE
- `discuss/` — architecture proposals, bug fix notes, roadmap documents

---

### 2. Reviewed Conflict Resolution in PR #4 (yzhou-nvwa/nvwa-mvp) ✅ CONFIRMED

**Context:** Yuxin Chen (chenyux3) created PR #4 (`merge/chen-reconcile` → `main`) to reconcile Chen's production snapshot with their main branch (which had a regression test harness + task router). He asked Chen to confirm `src/agent/` conflict resolutions.

**PR #4:** https://github.com/yzhou-nvwa/nvwa-mvp/pull/4

**Analysis performed:** Diffed `reconcile-from-production` vs `boss/merge/chen-reconcile` for all `src/agent/` files.

**Results:**

| File | Status | Detail |
|---|---|---|
| `core.py` | ✅ Identical to production | No changes |
| `tools.py` | ✅ Identical to production | No changes |
| `analysis_tools.py` | ✅ Confirmed OK | Added `seurat_clusters`, `cluster`, `clusters` to `_get_cluster_key()` fallback list |
| `prompts.py` | ✅ Confirmed OK | 3 additions: gene-missing behavior now requires STOP + user confirmation before retrying |
| `router.py` | 🆕 New file (their addition) | Rule-based intent classifier, not yet wired into core.py |

**Key prompts.py changes (their additions on top of Chen's version):**
- Gene Search intent mapping: "STOP — report error, list 1–3 alternatives, do NOT auto-select"
- Resilience Protocol: "STOP immediately. Do NOT retry with a different gene."
- Self-Healing split: gene-missing → stop + ask user; other failures → self-heal as before

**Confirmed and approved.** Yuxin will merge PR #4.

**Status:** ⏳ Waiting for Yuxin to merge PR #4 and make `yzhou-nvwa/nvwa-mvp` the canonical deployed repo.

**Note on GitHub diff UI:** The `prompts.py` changes appear buried in a large green block in the GitHub PR view because the entire file was new relative to `boss/main`. The precise diff was extracted by comparing the two branches directly locally.

---

## Previous Session Changes (2026-04-03)

### 4. Cell Composition Cross-Tabulation Bug Fix ✅ DEPLOYED (PR #4)

**Problem:** When users asked "How many cells per cell type in each condition?", the agent:
1. Incorrectly called `differential_expression()` first (irrelevant DEG table)
2. Fabricated a composition table with round numbers (4,500, 3,200...) instead of real counts (10,784, 2,891, 1,455...)
3. Said "cannot access data" when asked for exact counts after showing a plot
4. Showed raw serialized dataframe text in the UI ("Generated Code" section)

**Root Causes:**
1. No cross-tabulation tool existed — `inspect_metadata()` only does single-column `value_counts()`
2. `composition_plot()` and `cell_composition()` were separate tools with no data sharing — plot computed data then discarded it
3. `TableResult` constructor used wrong field name (`data=` instead of `csv_data=`, missing `display_df=`)
4. `display_df` used `to_string()` instead of `to_markdown()` → raw text rendered in Streamlit
5. Tool return value contained raw dataframe dumps, causing agent to echo them in response
6. No anti-fabrication guardrail — when composition tool failed, agent fell back to `inspect_metadata()` marginals and fabricated the full contingency table (mathematically invalid)

**Solution (Unified Pipeline):**
- `src/analysis/composition.py` — NEW: `cross_tabulate_metadata()` — Seurat-equivalent `groupby().size().unstack()`
- `src/plotting/composition.py` — NEW: `plot_composition()` — stacked bar chart, accepts pre-computed crosstab
- Replaced two separate tools (`cell_composition` + `composition_plot`) with ONE unified `composition_analysis()` tool:
  - Computes crosstab ONCE (single source of truth)
  - Returns both count table + stacked bar chart from same computation
  - Fixed `TableResult` constructor (`csv_data=`, `display_df=`)
  - `display_df` now uses `to_markdown()` for proper Streamlit rendering
  - Tool return value is short confirmation — no raw data dumps
- System prompt updates:
  - CRITICAL RULE #5: ANTI-FABRICATION GUARDRAIL — exact counts only from executed aggregation
  - Marginal distributions do NOT determine cross-tabulation — no fallback from `inspect_metadata()`
  - "Are my samples balanced?" now distinguishes total cell count balance vs cell-type composition balance
  - "Can you show me the exact number?" (follow-up) → refer to existing table, do NOT re-run
  - Composition balance interpretation must cite row-normalized percentages and largest shifts

**Files modified:**
- `src/analysis/composition.py` — NEW (83 lines)
- `src/plotting/composition.py` — NEW (103 lines), added `crosstab` parameter
- `src/agent/tools.py` — replaced `cell_composition` + `composition_plot` with `composition_analysis`; added `import io`, `import pandas as pd`
- `src/agent/prompts.py` — anti-fabrication rule + intent mapping + balance guidance

**PR:** https://github.com/lisperz/Nvwa-Lite/pull/4 (merged to main, squash)

**Deployment:** ✅ Deployed to EC2 2026-04-03

**Tested:**
- "How many cells per cell type in each condition?" → `composition_analysis()` called, real numbers shown ✓
- "Can you show me the exact number?" → refers to existing table, no re-run ✓
- "Are my samples balanced?" → asks user to clarify total counts vs composition ✓
- Table renders as proper Markdown table, no raw dumps ✓
- Anti-fabrication: if tool fails, agent reports error — does NOT synthesize contingency table ✓

---

### 3. Multi-Turn Visualization State Tracking + Split UMAP Legend Fix ✅ DEPLOYED

**Problem 1:** When users iteratively refine a UMAP plot across chat turns, the agent drops previously established constraints (e.g., `split_by` lost when user says "color by cell type instead").

**Problem 2:** Split UMAP plots showed incomplete legends (only 2 cell types instead of all 12) because the legend only displayed cell types present in the last panel.

**Root causes:**
1. Chat history is text-only — tool call parameters not preserved across turns
2. No visualization state tracking (unlike `DatasetState` for preprocessing)
3. No prompt guidance for refinement protocol
4. Split panel legend hardcoded to show only last panel's cell types

**Solution:**
- Created `src/agent/viz_state.py` — `VisualizationState` dataclass + singleton API (`bind_viz_state`, `get_viz_state`, `update_viz_state`)
- Updated 6 plot tools to call `update_viz_state()` after successful plot generation
- Added `LAST VISUALIZATION STATE` and `VISUALIZATION REFINEMENT PROTOCOL` sections to system prompt
- Threaded `viz_state` through `create_agent()` and persisted in `st.session_state`
- Fixed split UMAP legend to show ALL cell types across all panels with unified legend on right side

**Files modified:**
- `src/agent/viz_state.py` — NEW file (81 lines)
- `src/agent/tools.py` — added `update_viz_state()` calls in 6 plot tools
- `src/agent/prompts.py` — added viz state block + refinement protocol to template
- `src/agent/core.py` — threaded `viz_state` parameter through agent creation
- `src/ui/app.py` — persist `viz_state` in session state, reset on dataset load
- `src/plotting/executor.py` — unified legend for split UMAP plots showing all cell types

**PR:** https://github.com/lisperz/Nvwa-Lite/pull/3 (merged to main, squash)

**Deployment:** ✅ Deployed to EC2 via SSH at ~10:20 UTC 2026-04-03

**Tested:**
- Multi-turn refinement: "Show UMAP split by condition" → "Color by cell type instead" → preserves split ✓
- Legend display: Shows all 12 cell types in unified legend on right side ✓
- State reset: Upload new dataset → no stale viz state ✓

---

### 2. Fixed Inconsistent Response for Cell Distribution Queries ✅ DEPLOYED

**Root cause:** The word "distribution" is ambiguous — the LLM sometimes interpreted "distribution of cells across conditions" as QC metric violin plots (nCount_RNA/nFeature_RNA) instead of cell count breakdowns per condition.

**Fix:** Added explicit intent mapping in `src/agent/prompts.py` under `USER INTENT MAPPING`:
- Maps "cell distribution / how many cells per condition / are samples balanced / cell composition" queries to `inspect_metadata()`
- Includes negative instruction: "Do NOT use violin_plot for this — violin plots show expression/QC metric distributions, not cell composition"
- 5 concrete examples to anchor the LLM's behavior

**Files modified:**
- `src/agent/prompts.py` — added intent mapping entry

**Also fixed:** `scripts/start.sh` — added `source .env` so local dev starts with env vars loaded

**PR:** https://github.com/lisperz/Nvwa-Lite/pull/2 (merged to main, squash)

**Deployment:** ✅ Deployed to EC2 via SSH. `nvwa-lite` and `redis` containers running. Landing container exits (expected — port 80 held by host nginx).

**Tested:** Confirmed fix works locally before deploying.

---

### 1. Fixed Token Usage Counting and Cost Calculation ✅ DEPLOYED

**Root cause:** Two bugs were causing ~30x inflation in reported token usage and costs:

1. **Fan-out JOIN bug in analytics queries** — `get_overview()` and `get_user_breakdown()` in `analytics.py` joined `chat_messages` (many rows per session) with aggregated `token_usage`/`tool_executions` subqueries (one row per session). This multiplied token/tool counts by the number of messages per session. Example: 10 messages × 50K tokens = 500K counted instead of 50K.

2. **Token logging bug in core.py** — Extracted `total_tokens` from OpenAI response but stored it entirely as `output_tokens` with `input_tokens=0`. This caused all tokens to be charged at the output rate ($0.60/M) instead of splitting between input ($0.15/M) and output ($0.60/M) rates.

**Fix:**
- `src/monitoring/analytics.py` — rewrote `get_overview()` and `get_user_breakdown()` to use separate queries without JOINs, eliminating row multiplication
- `src/agent/core.py` — extract `prompt_tokens` and `completion_tokens` separately from OpenAI response metadata (`usage.get("prompt_tokens")` and `usage.get("completion_tokens")`) and log them correctly

**Files modified:**
- `src/monitoring/analytics.py` — fixed fan-out JOIN in overview and user breakdown queries
- `src/agent/core.py` — fixed token extraction to properly split input/output tokens

**Result:** Dashboard now matches OpenAI billing ($0.45 actual vs $14.59 previously reported for 7-day window)

**Deployment:** ✅ Deployed to EC2 via SSH, services rebuilt and restarted at 01:57 UTC 2026-04-03

**Note:** Historical data in database still shows `input_tokens=0` because it was logged before the fix. Only new sessions created after deployment will have correct input/output token split.

**PR:** https://github.com/lisperz/Nvwa-Lite/pull/1 (merged to main)

---

## Previous Session Changes (2026-03-31)

### 1. Fixed nginx Upload Limit in /app Location Block ✅ DEPLOYED

### 1. Fixed nginx Upload Limit in /app Location Block ✅ DEPLOYED

**Root cause:** `client_max_body_size 2000m;` was set at the server level but **missing inside the `/app` location block**. Nginx location blocks do not reliably inherit this setting, so uploads > 1 MB through `/app` were being rejected with HTTP 413.

**Fix:**
- Added `client_max_body_size 2000m;` explicitly inside the `location /app` block
- Also added `client_body_timeout 300s;` and `client_body_buffer_size 128k;` for slow uploads
- Increased `proxy_connect_timeout` from 60s to 300s

**Files modified:**
- `nginx/landing_ssl.conf` — added upload directives inside `/app` location block

**Deployment:** ✅ Deployed to EC2 via `bash scripts/update_nginx.sh`

---

### 2. Diagnosed Slow/Stuck Uploads — Root Cause: Cloudflare 100 MB Limit ✅ RESOLVED

**Root cause:** Cloudflare (Free plan) has a **100 MB upload size limit**. Since `nvwa.bio` was proxied through Cloudflare (orange cloud), any file > 100 MB was silently dropped or stalled, even though nginx and Streamlit were both configured for 2 GB.

**Investigation findings:**
- `curl -I https://nvwa.bio/app` confirmed `server: cloudflare` — traffic was Cloudflare-proxied
- Current upload flow is **proxy upload**: Browser → Cloudflare → nginx → Streamlit → S3 (double transfer, hits CF limit)
- `storage/service.py` already has `generate_upload_url()` for presigned S3 URLs, but it is **not wired to the UI** — this is the gap to fix in future

**Resolution (immediate):** User turned off Cloudflare proxy for `nvwa.bio` — set DNS record to "gray cloud" (DNS only). Traffic now goes directly to EC2 without Cloudflare proxying.
- ✅ Upload size limit removed (nginx limit is 2 GB, Streamlit limit is 2 GB)
- ⚠️ Cloudflare DDoS protection disabled (acceptable for pilot phase)

---

### 3. Upload Architecture Analysis & Future Plan 📋

**Current upload flow (proxy upload — inefficient):**
```
Browser → nginx → Streamlit → S3
```
File is received by Streamlit first (consumes EC2 bandwidth), then re-uploaded to S3.

**Key finding in codebase:**
- `src/storage/service.py` — `S3StorageService.generate_upload_url()` (lines 48-80) already generates presigned PUT URLs for direct browser→S3 uploads
- `src/ui/components.py` — `file_upload_widget()` uses `st.file_uploader()` (proxy upload) and calls `s3_service.upload_file()` (backend-to-S3 push)
- The presigned URL method exists but is **not connected to the UI**

**Future plan (Option 2 — Direct S3 Upload):**
- Replace `st.file_uploader()` with a custom HTML/JS component
- Browser requests presigned URL from backend
- Browser uploads directly to S3 (bypasses nginx/Cloudflare/Streamlit entirely)
- Backend loads file from S3 for processing
- Benefits: no size limits (5 TB), zero server bandwidth, scales to any number of users
- Estimated cost: ~$6/month vs $200/month for Cloudflare Business plan

**Reference documents:**
- `discuss/UPLOAD_ARCHITECTURE_PROPOSAL.md` — original architecture proposal
- `discuss/DIRECT_S3_UPLOAD_IMPLEMENTATION.md` — detailed implementation plan with code examples

---

## Previous Session Changes (2026-03-30)

### 1. Fixed nginx Upload Limit for .h5ad Files ✅ DEPLOYED

**Root cause:** nginx `client_max_body_size` defaulted to 1 MB, blocking any `.h5ad` file > 1 MB at the proxy layer (413 error), causing browser uploads to hang indefinitely (red dot in Streamlit sidebar).

**Fix:**
- Rewrote `nginx/landing_ssl.conf` with complete production configuration:
  - `client_max_body_size 2000m;` — allows 2 GB file uploads
  - Added WebSocket upgrade map: `map $http_upgrade $connection_upgrade`
  - `/app` proxy block with WebSocket headers, 3600s timeouts, `proxy_buffering off`, `proxy_request_buffering off`, `proxy_socket_keepalive on`
  - `/admin` proxy block with WebSocket headers, 300s timeouts
- Created `scripts/update_nginx.sh` — pushes config to EC2 and reloads nginx (no downtime)
- Created `scripts/db_tunnel.sh` — SSH tunnel script for local RDS access on port 5433

**Files modified:**
- `nginx/landing_ssl.conf` — complete rewrite with WebSocket support
- `scripts/update_nginx.sh` — nginx deployment script
- `scripts/db_tunnel.sh` — RDS SSH tunnel script

**Deployment:** ✅ Deployed to EC2, file uploads now work up to 2 GB

---

### 2. Fixed Feedback Widget for Production ✅ DEPLOYED & WORKING

**Root cause:** Feedback timer worked locally but not in production. Multiple issues:
1. **WebSocket connections dropping** — nginx was timing out idle WebSocket connections after 600s, preventing `st.fragment(run_every=2)` from auto-rerunning
2. **Missing callback invocation** — `_on_plot_generated()` callback wasn't being called when plots were generated
3. **Missing `log_feedback()` method** — DatabaseLogger was missing the feedback logging method, causing AttributeError on submission

**Fix:**
- Updated nginx WebSocket configuration:
  - Increased `proxy_read_timeout` and `proxy_send_timeout` to 3600s (1 hour)
  - Added `proxy_socket_keepalive on` to prevent idle connection drops
  - Changed `Connection` header from hardcoded `"upgrade"` to dynamic `$connection_upgrade` variable
- Added debug logging to trace callback flow:
  - `src/ui/feedback_trigger.py` — logs fragment execution, timer checks, and dialog triggers
  - `src/agent/tools.py` — logs callback invocation in `_store_and_return()`
- Added missing `log_feedback()` method to `DatabaseLogger`:
  - Inserts feedback responses into `feedback_responses` table
  - Handles q1_score (1-5 stars), q2_time_saved (4 options), q3_open_text (optional)

**Files modified:**
- `nginx/landing_ssl.conf` — WebSocket keepalive configuration
- `src/ui/feedback_trigger.py` — added debug logging
- `src/agent/tools.py` — added callback logging
- `src/db/logger.py` — added `log_feedback()` method

**Status:**
- ✅ Feedback dialog appears after 5 minutes in production
- ✅ Feedback submission works and saves to RDS
- ✅ Fragment auto-rerun works with proper WebSocket configuration
- ✅ Timer set to 5 minutes (300 seconds) as requested for production use

---

## Previous Session Changes (2026-03-29)

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
✓ Cell composition analysis: `composition_analysis()` unified tool — computes crosstab once, returns count table + stacked bar chart + optional % table
✓ Anti-fabrication guardrail: exact counts only from executed aggregation; fallback from marginals is blocked
✓ "Are my samples balanced?" distinguishes total cell count balance vs cell-type composition balance

### Known Issues
- None currently blocking

### Tool Count
27 tools total (added `composition_analysis`; removed `cell_composition` and `composition_plot`)

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

---

## Recent Changes (2026-03-30)

### Fix: .h5ad Upload Stuck on Red Dot ✅ DEPLOYED

**Root cause:** nginx `client_max_body_size` defaulted to 1 MB. Any `.h5ad` file > 1 MB was rejected at the nginx proxy layer (413 error), causing the browser upload to hang indefinitely (red dot in Streamlit sidebar). The Streamlit limit was already 2000 MB but nginx blocked it first.

**Secondary issue:** `landing_ssl.conf` in repo was missing the `/app` and `/admin` proxy blocks (they had been added manually on EC2 and drifted out of sync).

**Fix:**
- Rewrote `nginx/landing_ssl.conf` with complete config:
  - `client_max_body_size 2000m;` — allows 2 GB file uploads
  - `/app` proxy block with WebSocket headers, 600s timeouts, `proxy_buffering off`, `proxy_request_buffering off`
  - `/admin` proxy block with WebSocket headers, 300s timeouts
- Created `scripts/update_nginx.sh` — pushes config to EC2 and reloads nginx (no downtime)

**Files modified:**
- `nginx/landing_ssl.conf` — complete rewrite
- `scripts/update_nginx.sh` — new deploy script

**Deployment:** ✅ Deployed to EC2 via `bash scripts/update_nginx.sh`

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
