# Session Summary — 2026-03-23 (Latest Update)

## Recent Changes This Session (2026-03-23)

### 1. Landing Page Created

**Context:** Built a public-facing marketing landing page for the Nvwa service.

**Solution:**
- Created `landing/index.html` — standalone static HTML landing page
- Tech: Plain HTML/CSS, no framework, no build step
- Font: Plus Jakarta Sans (Google Fonts)
- Design system: Blue `#2563EB` primary, orange `#F97316` CTA, clean minimalist SaaS style
- Sections: Navbar, Hero (with animated demo window), Social proof logos, Features (6 cards), How it works (3 steps), Testimonials (3), Pricing (3 tiers), Final CTA, Footer
- All "Request Access" / "Request Early Access" / bottom CTA buttons link to `http://localhost:8501` (local) — update to production URL when deploying

**Files Created:**
- `landing/index.html` — full landing page

---

### 2. Landing Page Served via Nginx (Local Docker)

**Context:** Wired the landing page into the local Docker stack so the full user flow works locally.

**Solution:**
- Added `nginx/landing.conf` — nginx server block serving `landing/` on port 80
- Added `landing` service to `docker-compose.yml` (nginx:alpine, port 80, mounts `../landing/` and `nginx/landing.conf`)
- Added `scripts/start_landing.sh` — convenience script to start landing + redis + nvwa-lite together
- Updated `scripts/start_local_test.sh` to show landing URL in output

**Local user flow:**
```
http://localhost  →  landing page
  click "Request Access"
    →  http://localhost:8501  (Streamlit token auth)
```

**Files Created/Modified:**
- `nginx/landing.conf` — new nginx config
- `docker-compose.yml` — added `landing` service
- `scripts/start_landing.sh` — new convenience start script
- `scripts/start_local_test.sh` — updated access point output

---

### 3. Docker Cleanup — Freed ~9.1 GB

**Context:** Removed outdated and unused Docker images and build cache to free disk space.

**What was removed:**
- `nvwa-lite-dashboard:latest` image — 3.5 GB (dashboard container was unhealthy and not needed)
- Full build cache (17 layers) — 5.6 GB

**Commands used:**
```bash
docker stop nvwa-lite-dashboard-1 && docker rm nvwa-lite-dashboard-1
docker image rm nvwa-lite-dashboard:latest
docker builder prune -af
```

**Current active images:**
- `nvwa-lite-nvwa-lite:latest` — 3.5 GB (main app)
- `nginx:alpine` — 92 MB (landing page)
- `redis:7-alpine` — 61 MB (session store)

---

## Previous Session (2026-03-21)

### 1. Support "All Markers" in Marker Table Export

**Context:** When the user asked to "show all markers" for each cluster, the system only returned top 20 per cluster because `sc.tl.rank_genes_groups()` was called with a hardcoded `n_genes=20`.

**Root Cause:**
- `differential_expression` tool in `src/agent/tools.py` called `run_differential_expression()` without an `n_genes` parameter
- `run_differential_expression()` in `src/analysis/differential.py` defaulted to `n_genes=20`
- Scanpy only stores the top-N genes in `adata.uns["rank_genes_groups"]`, so "all" still meant only 20

**Solution:**
- Added `n_genes: int = 20` parameter to `differential_expression()` and `get_cluster_degs()` tools
- Added `n_genes=0` → `adata.n_vars` (all genes) guard at the top of `run_differential_expression()`
- Updated system prompt (`src/agent/prompts.py`) to guide agent to use `n_genes=0` when user asks for "all markers", "complete marker table", etc.
- Added workflow examples: `differential_expression(n_genes=0)` → `get_de_results_table()`

**Files Modified:**
- `src/analysis/differential.py` — Added `n_genes <= 0` → `adata.n_vars` guard
- `src/agent/tools.py` — Added `n_genes` param to `differential_expression` and `get_cluster_degs`
- `src/agent/prompts.py` — Updated intent mapping, DE decision logic, mode examples, and clarification protocol

---

### 2. Sidebar Cleanup — Removed Dataset Info Sections

**Context:** User requested removal of "Dataset Info", "Observation Keys", and "Common marker genes" sidebar sections. Only "Upload", "Pipeline", and "Example Queries" should remain.

**Solution:**
- Replaced `dataset_info_panel(adata, state)` with a lean `pipeline_panel(state)` function that only shows the processing pipeline status (e.g., `Normalized → PCA → UMAP → Clustered (louvain) → DE`)
- Removed unused `AnnData` import from `components.py`
- Updated `app.py` import and call site

**Files Modified:**
- `src/ui/components.py` — Replaced `dataset_info_panel` with `pipeline_panel`
- `src/ui/app.py` — Updated import and call from `dataset_info_panel(adata, state=ds_state)` to `pipeline_panel(ds_state)`

---

### 3. Fix Dotplot Showing Too Few Genes ("Top N Per Cluster")

**Context:** When user asked "Generate a dot plot for the top 3 markers of every cluster", the dotplot only showed 3 genes total (all from B cells), instead of top 3 per cluster × number of clusters.

**Root Cause:**
- `get_top_markers` returned a flat comma-separated list with a generic message
- The LLM agent interpreted "top 3 marker genes" as "3 genes total" and passed only 3 to `dotplot()`
- The tool output did not clearly show the per-cluster breakdown

**Solution:**
- Updated `get_top_markers` in `src/agent/tools.py` to:
  1. Show a per-cluster gene breakdown (e.g., `B cells: CD79A, CD74, CD79B`) before the full list
  2. Add `IMPORTANT: Use ALL genes below for dotplot/heatmap (not just one cluster's genes)` directive
  3. Still return the full flat comma-separated gene list at the end (unchanged for downstream use)
- Updated system prompt rule 4 to clarify: `get_top_markers(n_genes_per_cluster=N)` returns N genes **per cluster**, so all returned genes must be passed to dotplot/heatmap

**Files Modified:**
- `src/agent/tools.py` — Updated `get_top_markers` output format
- `src/agent/prompts.py` — Updated MANDATORY HEATMAP/DOTPLOT SOP note

---

### 4. Deployment

All three changes were committed (commit `1de4bc1`) and deployed to EC2 (`ubuntu@3.150.203.87`) via `scripts/deploy-to-ec2.sh`. Old Docker images were pruned locally before rebuild.

---

## Previous Session (2026-03-20)

### 1. Production Session Cache Reset (2026-03-20)

**Context:** User reported login issues when using different tokens. Users were being blocked from creating new sessions.

**Root Cause:**
- `SessionManager._has_active_session()` returns `True` for users with stale sessions (e.g., browser closed without logout)
- In-memory fallback has no TTL, so stale sessions persist until server restart
- `@st.cache_resource` on `get_auth_service()` and `get_session_manager()` in `src/ui/app.py:64,70` means token/session state is frozen until cache is cleared
- `create_session()` silently returns `None` when a user already has an active session, blocking login with no clear error

**Solution:**
- SSH'd into EC2 (`ubuntu@3.150.203.87`) using key at `/Users/zhuchen/Downloads/nvwa-key.pem`
- Flushed all Redis session data: `docker-compose exec -T redis redis-cli FLUSHALL`
- Restarted app container: `docker-compose restart nvwa-lite`
- This cleared both Redis sessions and Streamlit's `@st.cache_resource` cache

**Commands Used:**
```bash
ssh -i /Users/zhuchen/Downloads/nvwa-key.pem -T ubuntu@3.150.203.87 \
  "cd ~/Nvwa-Lite && docker-compose exec -T redis redis-cli FLUSHALL && docker-compose restart nvwa-lite"
```

**Key Learnings:**
- Redis `FLUSHALL` + container restart is the standard procedure to reset all login state in production
- The in-memory session fallback has no TTL — if Redis is down, stale sessions persist forever until restart
- Consider adding a session expiry mechanism to the in-memory fallback, or surfacing a clearer error when `create_session()` returns `None`

**Files Modified:** None (operational fix only)

---

## Previous Session (2026-03-18)

### 1. Production Auth Fix: pilot_tokens.json Not Accessible in Docker Container (2026-03-18)

**Context:** After deploying to AWS EC2, users could not log in with valid tokens. The same tokens worked fine on local machine.

**Root Cause:**
- `pilot_tokens.json` is listed in `.gitignore` (intentionally, as it contains sensitive auth data)
- Therefore it was never pushed to GitHub and never pulled to EC2 during deployment
- Even after manually creating the file on EC2, authentication still failed
- The real issue: `docker-compose.yml` did not mount `pilot_tokens.json` as a volume, so the file existed on the EC2 host but was invisible inside the Docker container
- The app (`src/ui/app.py:67`) looks for `pilot_tokens.json` using a relative path `Path("pilot_tokens.json")`, which resolves to `/app/pilot_tokens.json` inside the container

**Solution:**
- Added volume mount to `docker-compose.yml` for the `nvwa-lite` service:
  ```yaml
  - ./pilot_tokens.json:/app/pilot_tokens.json
  ```
- Copied corrected `docker-compose.yml` to EC2 via SCP using SSH key at `/Users/zhuchen/Downloads/nvwa-key.pem`
- Restarted containers with `docker-compose down && docker-compose up -d`

**Files Modified:**
- `docker-compose.yml` — Added `./pilot_tokens.json:/app/pilot_tokens.json` volume mount to `nvwa-lite` service

**Key Learnings:**
- `pilot_tokens.json` must be manually created on any new deployment target (EC2, etc.) since it is gitignored
- The file must also be volume-mounted in docker-compose.yml to be accessible inside the container
- SSH key for EC2 access is located at `/Users/zhuchen/Downloads/nvwa-key.pem`
- EC2 instance IP: `3.150.203.87`, user: `ubuntu`, project path: `~/Nvwa-Lite`

**New Script Created:**
- `scripts/deploy-tokens-to-ec2.sh` — Helper script to copy `pilot_tokens.json` to EC2 and restart the container. Usage: `./scripts/deploy-tokens-to-ec2.sh <path-to-ssh-key>`

---

## Previous Session (2026-03-14)

### 1. Critical Bug Fixes: DE Intent Disambiguation & Cluster Resolution (2026-03-14)

**Context:** Fixed two critical logic bugs in the single-cell analysis agent that were causing incorrect tool selection and scope propagation.

**Bug 1: Report/Export Scope Defaulting to Wrong Cluster**
- **Problem:** When user requested "DEGs for Cluster X", the system would run analysis for Cluster X but then export results for a default cluster (often B cells)
- **Root Cause:** `get_de_results_table()` had a hardcoded default `target_cluster="B cells"` parameter
- **Solution:**
  - Added `target_group` parameter to `run_differential_expression()` in `src/analysis/differential.py`
  - Created centralized `resolve_analysis_scope()` function in `src/analysis/cluster_resolution.py`
  - Added new `get_cluster_degs(cluster="X")` tool for single-cluster one-vs-rest analysis
  - Updated system prompt to clarify three distinct DE modes:
    - Mode 1: One-vs-rest for ALL clusters → `differential_expression()`
    - Mode 2: One-vs-rest for ONE cluster → `get_cluster_degs(cluster="X")`
    - Mode 3: Pairwise comparison → `compare_groups_de(group1, group2)`
- **Files Modified:**
  - `src/analysis/differential.py` — Added target_group parameter
  - `src/analysis/cluster_resolution.py` — New file with scope resolution logic
  - `src/agent/tools.py` — Added get_cluster_degs tool
  - `src/agent/prompts.py` — Updated DE intent disambiguation section

**Bug 2: "DEGs for Cluster X" Misinterpreted as Pairwise**
- **Problem:** User request "DEGs for Cluster 3" was being interpreted as pairwise comparison instead of one-vs-rest
- **Root Cause:** Ambiguous intent mapping in system prompt
- **Solution:**
  - Clarified intent recognition rules in system prompt:
    - "DEGs for X" = Mode 2 (one-vs-rest for X)
    - "Compare X and Y" = Mode 3 (pairwise)
  - Added explicit workflow examples for each mode
  - Emphasized scope propagation: resolve → analyze → export with same target

**Validation:**
- Created test scripts: `scripts/test_intent_disambiguation.py`, `scripts/test_cluster_resolution.py`
- Documented in: `docs/DE_INTENT_FIX_REPORT.md`, `docs/DE_INTENT_QUICK_REFERENCE.md`

### 2. System Prompt Refactor Attempt & Rollback (2026-03-14)

**Context:** Attempted to reduce example anchoring in system prompt but reverted due to behavior changes.

**What Happened:**
1. **Initial Refactor:** Removed specific biological examples (CD4/CD8 T cells, B cells, MS4A1) and replaced with generic placeholders
2. **Formatting Errors:** Encountered Python `.format()` IndexError due to unescaped braces in LaTeX example
3. **Root Cause:** Line 254 had `{count}` and `{total}` as single braces instead of escaped `{{{{count}}}}{{{{total}}}}`
4. **Conservative Refactor:** Created minimal refactor removing only most repetitive examples
5. **Behavior Impact:** User reported conservative refactor still changed agent behavior too much
6. **Rollback Decision:** Restored original working prompt with only the minimal brace escaping fix

**Final State:**
- System prompt restored to original version with all concrete examples
- Only fix preserved: Line 254 brace escaping to prevent format() errors
- Backup file deleted after successful restoration
- All refactor artifacts cleaned up

**Files Modified:**
- `src/agent/prompts.py` — Restored to original with brace fix

**Artifacts Cleaned Up:**
- Deleted 6 documentation files: CONSERVATIVE_REFACTOR_SUMMARY.md, PROMPT_REFACTORING_REPORT.md, PROMPT_REFACTORING_SUMMARY.md, PROMPT_RESTORATION_REPORT.md, REFACTORED_PROMPT_TESTING_GUIDE.md, STRING_FORMATTING_FIX.md
- Deleted backup: `src/agent/prompts_original_backup.py`

### 3. Streamlit Configuration Update (2026-03-14)

**Context:** Updated Streamlit configuration to support large file uploads and messages.

**Changes:**
- Added `maxMessageSize = 2000` to `.streamlit/config.toml`
- Kept `maxUploadSize = 2000` (already set)
- Ensures large datasets and analysis results can be handled

**Files Modified:**
- `.streamlit/config.toml` — Added maxMessageSize parameter

### 4. Docker Configuration Update (2026-03-14)

**Context:** Updated Docker setup to support live code updates during local testing.

**Changes:**
- Added source code volume mount: `./src:/app/src` to both nvwa-lite and dashboard services
- Rebuilt Docker images from scratch with `--no-cache`
- Cleaned up unused Docker resources (freed 11.79GB)

**Benefits:**
- Code changes now reflected immediately without rebuilding
- Faster development iteration
- Better local testing experience

**Files Modified:**
- `docker-compose.yml` — Added src volume mounts

## Current System State

### Working Features
✓ Three distinct DE analysis modes with correct intent disambiguation
✓ `n_genes=0` support for computing ALL markers (not just top 20)
✓ Cluster resolution supporting both numeric IDs and cell type names
✓ Scope propagation from analysis to export
✓ QC metrics handling (separate from gene analysis)
✓ System prompt with original examples and behavior
✓ Docker live code mounting for development
✓ Sidebar: Upload + Pipeline + Example Queries only (Dataset Info removed)
✓ Dotplot/heatmap shows correct per-cluster gene count

### Known Issues
- None currently blocking

### Tool Count
21 tools total (including get_cluster_degs added 2026-03-14)

### Documentation
- `docs/BUG_FIX_SUMMARY.md` — Overview of both bugs and fixes
- `docs/DE_INTENT_FIX_REPORT.md` — Detailed DE intent disambiguation fix
- `docs/DE_INTENT_QUICK_REFERENCE.md` — Quick reference for three DE modes
- `docs/CLUSTER_RESOLUTION_VALIDATION.md` — Cluster resolution validation

## Next Session Priorities

### High Priority
1. **Monitor "all markers" feature** — Validate `n_genes=0` produces correct full-gene results in production
2. **Monitor dotplot fix** — Confirm top-N-per-cluster dotplots show correct number of genes

### Medium Priority
1. **S3 Migration** — Local file storage on EC2 is not production-ready; migrate to S3
2. **Monitoring Setup** — Configure CloudWatch logs and error alerts

### Low Priority
1. **Streamlit Deprecation** — Update `use_container_width` to `width` parameter (warning in logs)
2. **Session Expiry** — Add TTL to in-memory session fallback to prevent stale session issues

## Key Infrastructure

- **EC2 IP:** `3.150.203.87`
- **EC2 User:** `ubuntu`
- **EC2 Project Path:** `/home/ubuntu/Nvwa-Lite`
- **SSH Key:** `/Users/zhuchen/Downloads/nvwa-key.pem`
- **Deploy Script:** `scripts/deploy-to-ec2.sh` (git pull + docker rebuild on EC2)
- **Local Test:** `scripts/start_local_test.sh`
- **GitHub Remote:** `origin` → `https://github.com/lisperz/Nvwa-Lite`

## Key Learnings (Cumulative)

1. **Prompt Stability:** Even conservative prompt changes can significantly impact agent behavior; prefer tool-level fixes
2. **Scope Propagation:** Critical to pass target identifiers through entire pipeline (resolve → analyze → export)
3. **Intent Disambiguation:** Explicit mode definitions and examples are necessary for correct tool selection
4. **Brace Escaping:** Python `.format()` requires `{{{{` to produce literal `{{` in output
5. **Agent Output Clarity:** When a tool returns data the agent must pass verbatim to another tool, include explicit directives like "IMPORTANT: Use ALL genes below" to prevent the LLM from summarizing or truncating
6. **Docker State:** Redis `FLUSHALL` + container restart resets all login/session state in production
7. **n_genes=0 Pattern:** Using 0 as a sentinel for "all genes" is clean — convert to `adata.n_vars` at the analysis layer, not the tool layer


## Recent Changes This Session (2026-03-20)

### 1. Production Session Cache Reset (2026-03-20)

**Context:** User reported login issues when using different tokens. Users were being blocked from creating new sessions.

**Root Cause:**
- `SessionManager._has_active_session()` returns `True` for users with stale sessions (e.g., browser closed without logout)
- In-memory fallback has no TTL, so stale sessions persist until server restart
- `@st.cache_resource` on `get_auth_service()` and `get_session_manager()` in `src/ui/app.py:64,70` means token/session state is frozen until cache is cleared
- `create_session()` silently returns `None` when a user already has an active session, blocking login with no clear error

**Solution:**
- SSH'd into EC2 (`ubuntu@3.150.203.87`) using key at `/Users/zhuchen/Downloads/nvwa-key.pem`
- Flushed all Redis session data: `docker-compose exec -T redis redis-cli FLUSHALL`
- Restarted app container: `docker-compose restart nvwa-lite`
- This cleared both Redis sessions and Streamlit's `@st.cache_resource` cache

**Commands Used:**
```bash
ssh -i /Users/zhuchen/Downloads/nvwa-key.pem -T ubuntu@3.150.203.87 \
  "cd ~/Nvwa-Lite && docker-compose exec -T redis redis-cli FLUSHALL && docker-compose restart nvwa-lite"
```

**Key Learnings:**
- Redis `FLUSHALL` + container restart is the standard procedure to reset all login state in production
- The in-memory session fallback has no TTL — if Redis is down, stale sessions persist forever until restart
- Consider adding a session expiry mechanism to the in-memory fallback, or surfacing a clearer error when `create_session()` returns `None`

**Files Modified:** None (operational fix only)

---

## Previous Session (2026-03-18)

## Recent Changes This Session (2026-03-18)

### 1. Production Auth Fix: pilot_tokens.json Not Accessible in Docker Container (2026-03-18)

**Context:** After deploying to AWS EC2, users could not log in with valid tokens. The same tokens worked fine on local machine.

**Root Cause:**
- `pilot_tokens.json` is listed in `.gitignore` (intentionally, as it contains sensitive auth data)
- Therefore it was never pushed to GitHub and never pulled to EC2 during deployment
- Even after manually creating the file on EC2, authentication still failed
- The real issue: `docker-compose.yml` did not mount `pilot_tokens.json` as a volume, so the file existed on the EC2 host but was invisible inside the Docker container
- The app (`src/ui/app.py:67`) looks for `pilot_tokens.json` using a relative path `Path("pilot_tokens.json")`, which resolves to `/app/pilot_tokens.json` inside the container

**Solution:**
- Added volume mount to `docker-compose.yml` for the `nvwa-lite` service:
  ```yaml
  - ./pilot_tokens.json:/app/pilot_tokens.json
  ```
- Copied corrected `docker-compose.yml` to EC2 via SCP using SSH key at `/Users/zhuchen/Downloads/nvwa-key.pem`
- Restarted containers with `docker-compose down && docker-compose up -d`

**Files Modified:**
- `docker-compose.yml` — Added `./pilot_tokens.json:/app/pilot_tokens.json` volume mount to `nvwa-lite` service

**Key Learnings:**
- `pilot_tokens.json` must be manually created on any new deployment target (EC2, etc.) since it is gitignored
- The file must also be volume-mounted in docker-compose.yml to be accessible inside the container
- SSH key for EC2 access is located at `/Users/zhuchen/Downloads/nvwa-key.pem`
- EC2 instance IP: `3.150.203.87`, user: `ubuntu`, project path: `~/Nvwa-Lite`

**New Script Created:**
- `scripts/deploy-tokens-to-ec2.sh` — Helper script to copy `pilot_tokens.json` to EC2 and restart the container. Usage: `./scripts/deploy-tokens-to-ec2.sh <path-to-ssh-key>`

---

## Previous Session (2026-03-14)

### 1. Critical Bug Fixes: DE Intent Disambiguation & Cluster Resolution (2026-03-14)

**Context:** Fixed two critical logic bugs in the single-cell analysis agent that were causing incorrect tool selection and scope propagation.

**Bug 1: Report/Export Scope Defaulting to Wrong Cluster**
- **Problem:** When user requested "DEGs for Cluster X", the system would run analysis for Cluster X but then export results for a default cluster (often B cells)
- **Root Cause:** `get_de_results_table()` had a hardcoded default `target_cluster="B cells"` parameter
- **Solution:**
  - Added `target_group` parameter to `run_differential_expression()` in `src/analysis/differential.py`
  - Created centralized `resolve_analysis_scope()` function in `src/analysis/cluster_resolution.py`
  - Added new `get_cluster_degs(cluster="X")` tool for single-cluster one-vs-rest analysis
  - Updated system prompt to clarify three distinct DE modes:
    - Mode 1: One-vs-rest for ALL clusters → `differential_expression()`
    - Mode 2: One-vs-rest for ONE cluster → `get_cluster_degs(cluster="X")`
    - Mode 3: Pairwise comparison → `compare_groups_de(group1, group2)`
- **Files Modified:**
  - `src/analysis/differential.py` — Added target_group parameter
  - `src/analysis/cluster_resolution.py` — New file with scope resolution logic
  - `src/agent/tools.py` — Added get_cluster_degs tool
  - `src/agent/prompts.py` — Updated DE intent disambiguation section

**Bug 2: "DEGs for Cluster X" Misinterpreted as Pairwise**
- **Problem:** User request "DEGs for Cluster 3" was being interpreted as pairwise comparison instead of one-vs-rest
- **Root Cause:** Ambiguous intent mapping in system prompt
- **Solution:**
  - Clarified intent recognition rules in system prompt:
    - "DEGs for X" = Mode 2 (one-vs-rest for X)
    - "Compare X and Y" = Mode 3 (pairwise)
  - Added explicit workflow examples for each mode
  - Emphasized scope propagation: resolve → analyze → export with same target

**Validation:**
- Created test scripts: `scripts/test_intent_disambiguation.py`, `scripts/test_cluster_resolution.py`
- Documented in: `docs/DE_INTENT_FIX_REPORT.md`, `docs/DE_INTENT_QUICK_REFERENCE.md`

### 2. System Prompt Refactor Attempt & Rollback (2026-03-14)

**Context:** Attempted to reduce example anchoring in system prompt but reverted due to behavior changes.

**What Happened:**
1. **Initial Refactor:** Removed specific biological examples (CD4/CD8 T cells, B cells, MS4A1) and replaced with generic placeholders
2. **Formatting Errors:** Encountered Python `.format()` IndexError due to unescaped braces in LaTeX example
3. **Root Cause:** Line 254 had `{count}` and `{total}` as single braces instead of escaped `{{{{count}}}}{{{{total}}}}`
4. **Conservative Refactor:** Created minimal refactor removing only most repetitive examples
5. **Behavior Impact:** User reported conservative refactor still changed agent behavior too much
6. **Rollback Decision:** Restored original working prompt with only the minimal brace escaping fix

**Final State:**
- System prompt restored to original version with all concrete examples
- Only fix preserved: Line 254 brace escaping to prevent format() errors
- Backup file deleted after successful restoration
- All refactor artifacts cleaned up

**Files Modified:**
- `src/agent/prompts.py` — Restored to original with brace fix

**Artifacts Cleaned Up:**
- Deleted 6 documentation files: CONSERVATIVE_REFACTOR_SUMMARY.md, PROMPT_REFACTORING_REPORT.md, PROMPT_REFACTORING_SUMMARY.md, PROMPT_RESTORATION_REPORT.md, REFACTORED_PROMPT_TESTING_GUIDE.md, STRING_FORMATTING_FIX.md
- Deleted backup: `src/agent/prompts_original_backup.py`

### 3. Streamlit Configuration Update (2026-03-14)

**Context:** Updated Streamlit configuration to support large file uploads and messages.

**Changes:**
- Added `maxMessageSize = 2000` to `.streamlit/config.toml`
- Kept `maxUploadSize = 2000` (already set)
- Ensures large datasets and analysis results can be handled

**Files Modified:**
- `.streamlit/config.toml` — Added maxMessageSize parameter

### 4. Docker Configuration Update (2026-03-14)

**Context:** Updated Docker setup to support live code updates during local testing.

**Changes:**
- Added source code volume mount: `./src:/app/src` to both nvwa-lite and dashboard services
- Rebuilt Docker images from scratch with `--no-cache`
- Cleaned up unused Docker resources (freed 11.79GB)

**Benefits:**
- Code changes now reflected immediately without rebuilding
- Faster development iteration
- Better local testing experience

**Files Modified:**
- `docker-compose.yml` — Added src volume mounts

## Current System State

### Working Features
✓ Three distinct DE analysis modes with correct intent disambiguation
✓ Cluster resolution supporting both numeric IDs and cell type names
✓ Scope propagation from analysis to export
✓ QC metrics handling (separate from gene analysis)
✓ System prompt with original examples and behavior
✓ Docker live code mounting for development

### Known Issues
- None currently blocking

### Tool Count
21 tools total (including get_cluster_degs added in this session)

### Documentation
- `docs/BUG_FIX_SUMMARY.md` — Overview of both bugs and fixes
- `docs/DE_INTENT_FIX_REPORT.md` — Detailed DE intent disambiguation fix
- `docs/DE_INTENT_QUICK_REFERENCE.md` — Quick reference for three DE modes
- `docs/CLUSTER_RESOLUTION_VALIDATION.md` — Cluster resolution validation

## Next Session Priorities

### High Priority
1. **Test DE Intent Disambiguation:** Validate that all three DE modes work correctly with real user queries
2. **Monitor Scope Propagation:** Ensure target cluster is correctly passed through entire pipeline
3. **Verify Cluster Resolution:** Test with both numeric IDs and cell type names

### Medium Priority
1. **Example Anchoring:** If anchoring issues persist, consider tool-level improvements rather than prompt changes
2. **Performance Optimization:** Monitor agent response time with restored prompt

### Low Priority
1. **Streamlit Deprecation:** Update `use_container_width` to `width` parameter (warning in logs)

## Key Learnings

1. **Prompt Stability:** Even conservative prompt changes can significantly impact agent behavior; prefer tool-level fixes
2. **Scope Propagation:** Critical to pass target identifiers through entire pipeline (resolve → analyze → export)
3. **Intent Disambiguation:** Explicit mode definitions and examples are necessary for correct tool selection
4. **Brace Escaping:** Python `.format()` requires `{{{{` to produce literal `{{` in output
5. **Testing First:** Always validate behavior changes before committing to refactors

## Files Modified This Session

### Core Logic
- `src/analysis/differential.py` — Added target_group parameter
- `src/analysis/cluster_resolution.py` — New file for scope resolution
- `src/agent/tools.py` — Added get_cluster_degs tool
- `src/agent/prompts.py` — Fixed brace escaping, restored original content

### Configuration
- `docker-compose.yml` — Added src volume mounts

### Documentation (New)
- `docs/BUG_FIX_SUMMARY.md`
- `docs/DE_INTENT_FIX_REPORT.md`
- `docs/DE_INTENT_QUICK_REFERENCE.md`
- `docs/CLUSTER_RESOLUTION_VALIDATION.md`

### Test Scripts (New)
- `scripts/test_intent_disambiguation.py`
- `scripts/test_cluster_resolution.py`
- `scripts/validate_bug_fixes.py`

### Cleanup
- Deleted 6 prompt refactor documentation files
- Deleted `src/agent/prompts_original_backup.py`

## Session Statistics

- **Bugs Fixed:** 2 critical logic bugs
- **New Tools Added:** 1 (get_cluster_degs)
- **New Files Created:** 7 (1 core logic, 3 test scripts, 3 docs)
- **Files Deleted:** 7 (6 docs, 1 backup)
- **Docker Cleanup:** 11.79GB freed
- **Prompt Iterations:** 3 (original → conservative refactor → rollback)
- **Final State:** Stable with original prompt + minimal fix
- **Deployment:** Successfully deployed to AWS EC2 production (172.31.16.191)

## Production Deployment (2026-03-14)

### Deployment Details
**Status:** ✅ Successfully deployed to AWS EC2

**Instance Information:**
- **Environment:** Production
- **Instance:** EC2 (172.31.16.191)
- **Location:** `/home/ubuntu/Nvwa-Lite`
- **Deployment Method:** SSH + deployment script
- **Deployment Script:** `scripts/deploy-to-ec2.sh`

**Deployment Process:**
1. Resolved git conflict with `.streamlit/config.toml`
2. Updated local config to match EC2 (added maxMessageSize = 2000)
3. Committed and pushed config changes to GitHub
4. Executed deployment script on EC2:
   - Reset local changes: `git reset --hard HEAD`
   - Pulled latest code: `git pull origin main`
   - Rebuilt Docker containers: `docker-compose build --no-cache`
   - Restarted services: `docker-compose up -d`

**Deployed Version:**
- **Commit:** `d644be8` (Update Streamlit config: add maxMessageSize = 2000)
- **Previous Commit:** `8533895` (Fix critical DE intent bugs and restore system prompt)

**What's Live in Production:**
- ✅ Bug fixes for DE intent disambiguation
- ✅ Cluster resolution with scope propagation
- ✅ Restored system prompt with original examples
- ✅ 21 tools including new `get_cluster_degs`
- ✅ Docker volume mounts for live code updates
- ✅ Streamlit config: 2000 MB upload/message limits
- ⚠️ Local file storage (S3 not configured yet)

**Post-Deployment Verification:**
- Container status: Running
- Services: nvwa-lite, dashboard, redis
- Logs: No errors reported
- Configuration: maxUploadSize = 2000, maxMessageSize = 2000

### Known Limitations in Production
1. **File Storage:** Using local filesystem (`data/uploads/`), not S3
   - Files stored on EC2 instance disk
   - No automatic cleanup
   - Not suitable for multi-instance scaling
   - Migration to S3 recommended before multi-user pilot

2. **Monitoring:** Basic Docker logs only
   - No CloudWatch integration
   - No automated alerts
   - Manual log review required

3. **Backup:** No automated backup configured
   - Code backed up in GitHub
   - Uploaded data not backed up
   - Manual backup recommended

### Production Access
- **Main App:** http://172.31.16.191:8501 (internal IP)
- **Dashboard:** http://172.31.16.191:8502 (internal IP)
- **Note:** External access depends on security group configuration

### Rollback Plan
If issues occur in production:
```bash
ssh ubuntu@172.31.16.191
cd ~/Nvwa-Lite
git log --oneline  # Find previous working commit
git checkout 8533895  # Or earlier commit
docker-compose down
docker-compose up -d
```

### Next Production Steps
1. **Monitor for 24-48 hours:**
   - Check logs regularly: `docker-compose logs -f`
   - Monitor dashboard metrics
   - Watch for any errors or performance issues

2. **S3 Migration (High Priority):**
   - Create S3 bucket
   - Configure AWS credentials
   - Update code to use S3StorageService
   - Test upload/download flow
   - Migrate existing files

3. **Monitoring Setup (Medium Priority):**
   - Configure CloudWatch logs
   - Set up error alerts
   - Create performance dashboards

4. **Backup Strategy (Medium Priority):**
   - Automated S3 backup for uploaded data
   - Database backup if using persistent storage
   - Configuration backup

5. **Security Hardening (Before Multi-User Pilot):**
   - Review security group rules
   - Enable HTTPS/SSL
   - Configure proper authentication
   - Set up rate limiting

