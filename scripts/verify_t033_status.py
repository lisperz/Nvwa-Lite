"""T-033 verification — real LLM, real agent, real erroring tool call.

Step 5 of [t033_truthful_tool_status_plan.md](../local/product/t033_truthful_tool_status_plan.md).

Invokes a real agent with a prompt that deterministically triggers a tool
returning an ``Error:`` string, then inspects the EventLogger JSONL for that
invocation's ``status`` field. Exit code 0 if any tool_execution log row has
``status='error'``; 1 otherwise.

Usage:
    conda activate nvwa-ai
    python scripts/verify_t033_status.py

The script writes tool_execution.log into a fresh tmp dir per run so runs
don't contaminate each other.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


# Use a prompt that deterministically resolves to a single tool call that errors.
# Asking for a UMAP colored by a nonexistent obs column is unambiguous — there's
# no preflight lookup (unlike genes, which route through lookup_gene first) — and
# plot_umap raises ValueError, which umap_plot catches and turns into "Error: ...".
ERRORING_PROMPT = (
    "Please make a UMAP colored by the column named "
    "'nvwa_t033_fake_col_xyz'. I know it exists, just make the plot."
)

DATASET_PATH = REPO_ROOT / "local" / "data" / "pbmc_test.h5ad"


def _reset_event_logger_handlers() -> None:
    """EventLogger uses module-level named loggers; clear handlers so a fresh
    log_dir takes effect. (Mirrors tests/unit/test_logging.py fixture.)"""
    for name in ["nvwa.tool_execution", "nvwa.user_interaction", "nvwa.system_metrics"]:
        lg = logging.getLogger(name)
        for h in lg.handlers[:]:
            h.close()
            lg.removeHandler(h)


def main() -> int:
    if "OPENAI_API_KEY" not in os.environ:
        print("ERROR: OPENAI_API_KEY not set. Activate nvwa-ai conda env.", file=sys.stderr)
        return 2
    if not DATASET_PATH.exists():
        print(f"ERROR: dataset not found: {DATASET_PATH}", file=sys.stderr)
        return 2

    from src.agent.core import create_agent
    from src.domain.analysis.h5ad_loader import load_h5ad
    from src.platform.observability import events as logging_service
    from src.core.types import detect_dataset_state

    _reset_event_logger_handlers()

    with tempfile.TemporaryDirectory(prefix="t033-verify-") as tmp:
        tmp_path = Path(tmp)
        # Redirect EventLogger's default log dir into the tmp dir.
        original_default = logging_service.EventLogger.__init__.__defaults__
        logging_service.EventLogger.__init__.__defaults__ = (tmp_path,)
        try:
            adata = load_h5ad(DATASET_PATH)
            state = detect_dataset_state(
                adata, source=str(DATASET_PATH), filename=DATASET_PATH.name
            )
            agent = create_agent(
                adata,
                api_key=os.environ["OPENAI_API_KEY"],
                model="gpt-4o-mini",
                dataset_state=state,
                user_id="t033-verify",
                session_id="t033-verify-session",
            )
            resp = agent.invoke(ERRORING_PROMPT, filename=DATASET_PATH.name)
        finally:
            logging_service.EventLogger.__init__.__defaults__ = original_default

        log_file = tmp_path / "tool_execution.log"
        if not log_file.exists():
            print("FAIL: no tool_execution.log was written — agent did not call any tool.")
            print(f"Agent response: {resp.text[:300]}")
            return 1

        entries = [
            json.loads(line)
            for line in log_file.read_text().splitlines()
            if line.strip()
        ]
        if not entries:
            print("FAIL: tool_execution.log is empty.")
            return 1

        print(f"Captured {len(entries)} tool_execution log row(s):")
        for i, e in enumerate(entries, 1):
            payload = e.get("payload", {})
            tool = e.get("task_type")
            status = payload.get("status")
            error = payload.get("error")
            result_head = (payload.get("result") or "")[:120]
            print(f"  [{i}] tool={tool} status={status!r} error={error!r}")
            print(f"        result[:120]={result_head!r}")

        error_rows = [e for e in entries if e["payload"].get("status") == "error"]
        if error_rows:
            print(
                f"\nPASS: {len(error_rows)} row(s) correctly logged status='error'."
            )
            return 0
        print(
            "\nFAIL: all rows logged status='success' despite Error: tool return. "
            "This is the pre-T-033 bug."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
