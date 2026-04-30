# Infrastructure Architecture

How the prod stack is wired today (snapshot: 2026-04-26). Two layers: an in-Docker compose stack defined in this repo, and a host-level fronting layer (nginx + SSL) that lives only on the EC2 box.

```
  Internet (nvwa.bio)
       │
       ▼ :443 (HTTPS)
  ┌─────────────────────────────────────────┐
  │ EC2 host (3.150.203.87)                 │
  │                                         │
  │   nginx  ← lives ON the host, NOT in    │
  │     │      Docker. Outside the repo.    │
  │     │                                   │
  │     ├── /app         → :8501  ┐         │
  │     ├── /admin       → :8502  │ Docker  │
  │     └── /api/cleanup → :8503  │ compose │
  │                               ▼         │
  │   ┌───────────────────────────────────┐ │
  │   │ docker-compose                    │ │
  │   │   nvwa-lite  (Streamlit + cleanup)│ │
  │   │   dashboard  (Streamlit admin)    │ │
  │   │   redis      (session store)      │ │
  │   └───────────────────────────────────┘ │
  └─────────────────────────────────────────┘
       │
       ├──▶ S3 bucket    (file storage)
       ├──▶ RDS Postgres (analytics, chat history)
       └──▶ OpenAI API   (LLM calls)
```

## nginx — what it does, why it's there

nginx is a reverse proxy on the EC2 host. It does four things Streamlit can't do well on its own:

1. **HTTPS termination** — owns the Let's Encrypt cert for `nvwa.bio` + `www.nvwa.bio`. Browser talks HTTPS to nginx; nginx talks plain HTTP inside the box.
2. **Path routing under one domain** — `/app` → main app (8501), `/admin` → dashboard (8502), `/api/cleanup` → cleanup endpoint (8503), `/` → static landing page on disk.
3. **WebSocket proxying** — Streamlit needs a persistent WebSocket back to the browser; nginx is configured with `Upgrade` headers to proxy it.
4. **Upload size + timeouts** — `client_max_body_size 2000m` (2 GB, matches .h5ad uploads), 1-hour read/send timeouts so long jobs aren't cut off.

The host nginx config currently lives at `/etc/nginx/sites-enabled/nvwa.bio` on the EC2 box. The source config is version-controlled in this repo at `infra/nginx/nvwa.bio.conf`.

## Streamlit plumbing — the dance

Streamlit re-runs your Python script every time the user clicks anything. To make that work behind nginx, three things have to line up:

1. **`baseUrlPath = "/app"`** ([.streamlit/config.toml](../.streamlit/config.toml)) — Streamlit needs to know it's served at `/app/`, not `/`, otherwise asset URLs come out wrong.
2. **WebSocket** — every interaction goes over a WebSocket back to the Python process. nginx must proxy upgrades. Misconfigured = app loads but nothing reacts.
3. **CORS / XSRF** — `enableXsrfProtection=true` + `enableCORS=false` is the safe combo when behind a single domain. Don't change them without understanding both.

The `dashboard` service is a separate Streamlit instance — same image, launched with `--server.baseUrlPath=/admin` and `--server.port=8502`.

**Operational implication:** any deploy change touching nginx routes, `baseUrlPath`, or WebSocket headers can break the app silently. Test the full path (browser → nginx → Streamlit) after any such edit, not just `curl localhost:8501`.

## Redis — what it does, why it's there

Streamlit's `st.session_state` is per-process, in-memory. That breaks if (a) multiple containers split state, (b) the container restarts (state evaporates), (c) another service needs the session (e.g., cleanup endpoint).

Redis solves this. Per-user/session state lives in keys like `session:<session_id>`; `SessionManager` is the writer.

**What's in Redis today:**
- Session metadata (user_id, session_id, dataset name, timestamps)
- NOT chat history — that's `st.session_state` (live) + RDS Postgres (persistence)
- NOT analytics — RDS

**Caveat:** the code has a Redis branch AND an in-memory fallback. Redis being healthy ≠ Redis actually being used.

## End-to-end request flow (upload + analyze)

1. Browser → `https://nvwa.bio/app`
2. nginx :443 terminates SSL, sees `/app`, proxies to `localhost:8501`
3. nvwa-lite (Streamlit) returns HTML; browser opens WebSocket back to `/app/_stcore/stream`
4. nginx upgrades the WebSocket
5. User picks a 1.5 GB .h5ad → Streamlit's `file_uploader` reads bytes server-side
6. Streamlit writes to `data/uploads` volume, then `S3StorageService.upload_file` → S3
7. SessionManager writes session state to Redis (`session:<id>`)
8. RDS Postgres logs `user_message` + `session_event`
9. Streamlit calls OpenAI API for LLM response
10. Tool runs (Scanpy/AnnData), result rendered as a plot
11. RDS logs `tool_calls`, `chat_messages` (assistant), `artifacts`

Step 5-6 (server-side double write to disk + S3) is currently a known bottleneck that has caused EC2 disk-fill incidents.

## Port 8503 — Cleanup Endpoint

The cleanup server runs inside the `nvwa-lite` container as a background HTTP server thread (started from `src/ui/app.py:109`). It listens on port 8503 internally.

In `docker-compose.yml`, port 8503 is bound to `127.0.0.1:8503:8503` (localhost-only). This means:
- The cleanup endpoint is NOT publicly reachable from the internet
- Only processes on the EC2 host can reach `localhost:8503`
- Host nginx proxies `https://nvwa.bio/api/cleanup` → `http://localhost:8503/cleanup`

Browser `sendBeacon` calls hit `https://nvwa.bio/api/cleanup`, which nginx proxies to the cleanup server inside the Docker network.

**Security note:** Port 8503 must remain bound to `127.0.0.1` only. Binding to `0.0.0.0:8503` would expose the raw cleanup endpoint publicly, bypassing nginx.

## Pieces NOT in the repo

The following live on the EC2 host or in external services, not version control:

- Host nginx config — source is in `infra/nginx/nvwa.bio.conf`; deployed to `/etc/nginx/sites-enabled/nvwa.bio` on EC2
- SSL certs (Let's Encrypt at `/etc/letsencrypt/live/nvwa.bio/`)
- Static landing page (`/var/www/nvwa.bio/index.html`)
- Populated `.env` with secrets (`OPENAI_API_KEY`, `DATABASE_URL`, `S3_BUCKET_NAME`, `ADMIN_PASSWORD`)
- AWS-side state — S3 bucket policy, lifecycle rules, RDS config, security groups
- Deploy procedure (manual `git pull + docker compose restart`)
