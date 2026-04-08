# Nvwa-Bio MVP — Claude Code Context

## Project
scRNA-seq visualization agent. Users upload `.h5ad` files and query the agent in natural language to produce plots, DE analysis, QC summaries, and biological interpretation.

**Stack:** Python, Streamlit, LangChain, OpenAI gpt-4o-mini, Scanpy, AnnData, Redis, S3.

---

## Team
- **CTO (you)** — architecture, quality systems, task router, regression harness
- **Chen (GitHub: lisperz)** — repo owner, full-stack, Layer 1 owner (tools, analysis, data pipeline), graduating May 2026. Reviews all PRs touching `src/agent/` or `src/analysis/`
- **CEO (Yalu)** — product, customer delivery, biological interpretation review (Layer 4)

---

## Architecture — Layer System

```
Incoming prompt
      ↓
  Task Router  ← entry point (src/agent/router.py)
   /        \
 2a          2b
  ↓            ↓
Tool call    Linguistic response
  ↓            ↓
Artifact     Human review only
  ↓          (QA Checklist — local/product/Nvwa_Layer2b_QA_Checklist.docx.md)
Gatekeeper
(automated checks — not yet built)
```

**Layer 1** — Tool execution (Scanpy, DE analysis, plotting). Owned by Chen.
**Layer 2a** — Tool selection and parameter passing. Automated via regression suite + gatekeeper (pending).
**Layer 2b** — Intent understanding, clarification, output honesty. Human QA checklist only.
**Layer 3** — UI/session layer.
**Layer 4** — Scientific correctness. Human-only, requires domain expertise. CEO reviews before each customer delivery.

**Key principle:** Layer 2a and 2b are two dimensions of a failure, not two buckets. A 2b root cause (misunderstood intent) can produce a 2a symptom (wrong tool called). Fix root cause, not symptom.

---

## Task Router (src/agent/router.py)
- Rule-based, two-tier: detect 2b linguistic patterns first, then match tool keywords
- Returns `RouterResult(layer, task_type, confidence, matched_on)`
- v0: classification + logging only — no behavior change to the agent
- v1 (pending): wire into `core.py`, log via EventLogger
- v2 (pending): behavioral forking — 2b path requires clarification before acting

---

## Key Files
- `src/agent/router.py` — task router (new)
- `src/agent/core.py` — agent loop (LangChain, tool-calling)
- `src/agent/tools.py` — 30+ tools, `get_all_tools()` returns full list
- `src/agent/prompts.py` — system prompt and intent rules
- `src/agent/analysis_tools.py` — advanced reasoning tools
- `src/logging/service.py` — EventLogger (log_tool_execution, log_session_event, etc.)
- `tests/integration/tests.yaml` — regression test cases
- `tests/unit/test_router.py` — router unit tests (84 cases, all passing)
- `scripts/run_tests.py` — regression runner, saves to local/reports/regression/
- `scripts/run_unit_tests.sh` — pytest wrapper, saves to local/reports/unit/

---

## PR Rules
- All changes on feature branches — main is always stable
- Branch naming: `feat/`, `fix/`, `chore/`
- Tag reviewer on Slack when PR opens
- 24-hour review SLA — self-merge allowed after 24h for tests/docs/scripts
- Mandatory approval from Chen for any `src/agent/` or `src/analysis/` changes
- PR description must include: what, why, how to test
- **Rebase before opening a PR:** always run `git fetch origin && git rebase origin/main` immediately before pushing and opening the PR. If main moved while the PR is open and touches the same files, rebase again before merging.

**Ownership rule:** Chen owns Layer 1 (accountability + review rights, not write lock). CTO can fix clearly scoped Layer 1 bugs; Chen reviews the PR.

---

## Repo Conventions
- **Timezone:** Pacific time for all docs, meeting notes, and daily summary filenames; UTC for code, logs, and DB timestamps.
- `tests/integration/` — integration test cases (tests.yaml)
- `tests/unit/` — unit tests (test_router.py, future: test_gatekeeper.py)
- `local/` — gitignored, machine-local only (reports, strategy docs, product docs)
- `test-records/` — tracked, curated test snapshots for key milestones

---

## Constraints
- Never add Co-authored-by trailers; never commit; CTO commits manually
- Do not refactor architecture without discussion
- Keep PRs small and focused — one logical change per PR
- Rule-based router stays as the foundation; LLM layer added on top later
- Do not change agent behavior until router + gatekeeper are wired in

---

## Daily Summary Habit
When asked to write a work summary for today:
1. First read the previous day's summary from `local/daily/` to avoid overlap or gaps
2. Write today's summary to `local/daily/YYYY-MM-DD.md` (Pacific time date)
3. Cover: what was done, decisions made (with reasoning), and what's pending

## Current Status (as of 2026-04-06)
- All work merged to main: regression test harness v0, task router v0, Chen's production reconcile
- **Pending next:** wire router into core.py (logging only), build gatekeeper
- **Known bugs:** HP-03/OOS-03 — agent adds unsolicited feature plot when gene mentioned (Layer 2b, discuss with Chen before fixing)
- **Fixed:** HP-06 — QC tool now returns TableResult instead of plain text
