# Session Summary — 2026-03-24 (Latest Update)

## Recent Changes This Session (2026-03-24)

### 1. Landing Page Cleanup

**Changes:**
- Removed Testimonials section (3 review cards) and Pricing section (3 tier cards)
- Removed unused CSS blocks and "Pricing" navbar link
- Replaced 5 institution logos with only "Northwestern University"
- Updated all "Request Access" button hrefs from `http://localhost:8501` → `https://nvwa.bio/app`

**Files Modified:**
- `landing/index.html`

---

### 2. Production Deployment: nvwa.bio Live

**Context:** Made `https://nvwa.bio` serve the landing page, with `/app` routing to the Streamlit app.

**Changes:**
- `nvwa-lite/.streamlit/config.toml` — Added `baseUrlPath = "/app"` to `[server]` section
- `nvwa-lite/docker-compose.yml` — Added `.streamlit` volume mount so config reaches container
- EC2 system nginx (`/etc/nginx/sites-enabled/nvwa`) — Serves landing at `/`, proxies `/app` → port 8501, proxies `/admin` → port 8502

**Files Modified:**
- `nvwa-lite/.streamlit/config.toml`
- `nvwa-lite/docker-compose.yml`
- EC2 `/etc/nginx/sites-enabled/nvwa` (SSH only, not in git)

---

### 3. AWS RDS + DatabaseLogger Integration

**Context:** Set up RDS PostgreSQL and wired event logging into the agent.

**DB Details:**
- Endpoint: `nvwa-pilot-db.cj42ock2ykh4.us-east-2.rds.amazonaws.com`
- Database: `nvwa`, User: `nvwa_admin`
- Password in EC2 `.env` as full DSN in `DATABASE_URL`
- 5 tables: `users`, `analysis_sessions`, `tool_executions`, `token_usage`, `chat_messages`

**New Files:**
- `src/db/__init__.py` — Package marker
- `src/db/client.py` — `ThreadedConnectionPool` singleton, `get_conn()` context manager
- `src/db/logger.py` — `DatabaseLogger` mirroring `EventLogger`, fire-and-forget
- `scripts/seed_users.py` — Seeds 10 pilot users (SHA-256 hashed tokens) into `users` table

**Modified Files:**
- `src/agent/core.py` — Dual-writes to file logs AND RDS
- `src/ui/app.py` — `agent.invoke()` passes `filename=ds_state.filename`
- `pyproject.toml` — Added `psycopg2-binary>=2.9.0`

**Key learnings:**
- `users` table stores `token_hash` (SHA-256), not raw token
- `#` in DATABASE_URL value is fine in `.env` file but breaks shell interpolation
- `docker-compose restart` does NOT pick up new volumes — must stop + rm + up

---

### 4. Admin Monitoring Dashboard

**Context:** Rewrote the dashboard to query RDS directly instead of parsing log files.

**URL:** `https://nvwa.bio/admin` (password-protected)
**Admin password:** `NvwaAdmin2026Secure` (stored in EC2 `.env` as `ADMIN_PASSWORD`)

**4-tab layout:**
- **Overview** — 6 KPI tiles (users, sessions, messages, tool calls, errors, estimated cost) + hourly bar chart + token breakdown
- **Users** — Per-user table with sessions, messages, tool calls, tokens, cost, last seen; drill into sessions and replay chat history
- **Tools** — Bar chart by call count, error rate per tool, avg/max latency, recent error log
- **Sessions** — Last 50 sessions: dataset filename, message count, duration, status

**New/modified files:**
- `src/monitoring/analytics.py` — Rewritten to query RDS (replaces file-log parsing)
- `src/monitoring/dashboard.py` — Rewritten with 4-tab layout + admin password gate
- `src/monitoring/__init__.py` — Cleared stale `AnalyticsService` import
- `docker-compose.yml` — Dashboard command now includes `--server.baseUrlPath=/admin`

**Bugs fixed during this work:**
1. `ImportError: cannot import name 'AnalyticsService'` — `__init__.py` still imported the deleted class → cleared to bare package marker
2. `TypeError: unsupported operand type(s) for *: 'decimal.Decimal' and 'float'` — PostgreSQL `SUM()` returns `decimal.Decimal`; fixed in `_query()` by casting all `Decimal` values to `float` as rows come back

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
✓ `https://nvwa.bio` → landing page (Northwestern University logo only)
✓ `https://nvwa.bio/app` → Streamlit app (token auth)
✓ `https://nvwa.bio/admin` → Admin dashboard (password-protected)
✓ RDS PostgreSQL DB logging (dual-write: file + DB)
✓ 10 pilot users seeded in RDS `users` table

### Known Issues
- None currently blocking

### Tool Count
21 tools total

### Infrastructure
- **EC2:** `ubuntu@3.150.203.87`, project at `/home/ubuntu/Nvwa-Lite`
- **SSH key:** `/Users/zhuchen/Downloads/nvwa-key.pem`
- **RDS:** `nvwa-pilot-db.cj42ock2ykh4.us-east-2.rds.amazonaws.com`
- **S3 bucket:** `nvwa-mvp-pilot` (us-east-2) — IAM role attached, not yet wired into file uploads
- **GitHub:** `https://github.com/lisperz/Nvwa-Lite`
- **Domain:** `https://nvwa.bio` via Cloudflare

### EC2 .env keys (summary)
- `DATABASE_URL` — full PostgreSQL DSN
- `ADMIN_PASSWORD` — dashboard login (`NvwaAdmin2026Secure`)
- `AWS_REGION=us-east-2`
- `S3_BUCKET_NAME=nvwa-mvp-pilot`
- `OPENAI_API_KEY` — agent LLM

## Next Session Priorities

### High Priority
1. **Verify DB logging live** — After a real user interaction, confirm rows appear in `tool_executions`, `token_usage`, `chat_messages` via the admin dashboard
2. **Change admin password** — Replace default `NvwaAdmin2026Secure` with something personal

### Medium Priority
1. **Wire S3 for file uploads** — Replace local `data/uploads/` with `S3StorageService` (already written)
2. **Narrow IAM S3 policy** — Replace `AmazonS3FullAccess` with scoped JSON policy on `nvwa-ec2-role`

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
