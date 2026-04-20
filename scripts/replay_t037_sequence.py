"""T-037 hypothesis check — replay Yalu's 10-turn sequence N times, record tool calls.

Asserts at turns 8 and 10 (the Fabp4 prompts that failed in Yalu's session 2026-04-17)
whether the agent selects a subset_* variant with a bogus subset_value.

Usage:
    conda activate nvwa-ai
    python scripts/replay_t037_sequence.py --iter 2 --out local/reports/t037_main.md
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

YALU_SEQUENCE = [
    "Hi can you tell me more about this dataset",
    "can you tell me the quality of this data",
    "Can you show me the violin plots of the QC?",
    "Can you show me the QC split by sample?",
    "Looks good. Can you show me the umap ",
    "I want umap with cell type name",
    "Show Fabp4 mRNA level",
    "Can you draw the Fabp4 in violin plot with cell type",
    "Can you draw Fabp4 violin plot with cell type",
    "can you show Fabp4 in feature plot, but split by condition",
]
BUG_TURNS = {8: "subset_violin_plot", 10: "subset_feature_plot"}


class ToolCallCapture(logging.Handler):
    def __init__(self) -> None:
        super().__init__()
        self.calls: list[dict] = []

    def emit(self, record: logging.LogRecord) -> None:
        if isinstance(record.msg, str) and record.msg.startswith("Tool call:"):
            args = record.args
            if args and len(args) >= 2:
                self.calls.append({"tool": args[0], "args": args[1]})


def classify_turn(
    turn_idx: int,
    calls: list[dict],
    celltypes: set[str],
    var_names: set[str],
) -> tuple[str, str]:
    """Return (verdict, detail) for a given turn.

    verdict: 'PASS' | 'BUG' | 'NO_TOOL' | 'OTHER'
    """
    if not calls:
        return "NO_TOOL", "no tool invoked"

    bug_tool = BUG_TURNS.get(turn_idx)
    if bug_tool is None:
        return "PASS", ", ".join(c["tool"] for c in calls)

    for c in calls:
        if c["tool"] == bug_tool:
            sv = c["args"].get("subset_value", "")
            if sv in var_names or sv not in celltypes:
                return "BUG", f"{bug_tool} subset_value={sv!r}"
    tools = ", ".join(c["tool"] for c in calls)
    return "PASS", tools


def run_once(dataset_path: Path, model: str) -> list[list[dict]]:
    """Run Yalu's sequence once. Return list of tool_calls per turn."""
    from src.agent.core import create_agent
    from src.agent.tools import clear_plot_results, clear_table_results
    from src.analysis.h5ad_loader import load_h5ad
    from src.types import detect_dataset_state

    clear_plot_results()
    clear_table_results()

    adata = load_h5ad(dataset_path)
    state = detect_dataset_state(adata, source=str(dataset_path), filename=dataset_path.name)
    agent = create_agent(
        adata,
        api_key=os.environ["OPENAI_API_KEY"],
        model=model,
        dataset_state=state,
    )

    capture = ToolCallCapture()
    core_logger = logging.getLogger("src.agent.core")
    core_logger.addHandler(capture)
    core_logger.setLevel(logging.INFO)

    per_turn: list[list[dict]] = []
    history: list[tuple[str, str]] = []
    try:
        for prompt in YALU_SEQUENCE:
            capture.calls.clear()
            resp = agent.invoke(prompt, chat_history=history, filename=dataset_path.name)
            per_turn.append(list(capture.calls))
            history.append(("user", prompt))
            history.append(("assistant", resp.text))
    finally:
        core_logger.removeHandler(capture)
    return per_turn


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iter", type=int, default=2)
    parser.add_argument(
        "--sleep-between-iter",
        type=int,
        default=75,
        help="Seconds to sleep between iterations to stay under OpenAI TPM cap (one 10-turn iter burns ~150-190k of the 200k/min window).",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        default=REPO_ROOT / "local" / "data" / "lung_wt_cko_tko_slim.h5ad",
    )
    parser.add_argument("--model", default="gpt-4o-mini")
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "local" / "reports" / f"t037_replay_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
    )
    args = parser.parse_args()

    if not args.dataset.exists():
        print(f"ERROR: dataset not found: {args.dataset}", file=sys.stderr)
        return 2
    if "OPENAI_API_KEY" not in os.environ:
        print("ERROR: OPENAI_API_KEY not set", file=sys.stderr)
        return 2

    from src.analysis.h5ad_loader import load_h5ad

    adata = load_h5ad(args.dataset)
    celltypes = set(adata.obs["celltype"].astype(str).unique()) if "celltype" in adata.obs.columns else set()
    var_names = set(adata.var_names)
    del adata

    commit = os.popen("git rev-parse --short HEAD").read().strip()

    rows: list[dict] = []
    started = datetime.now(timezone.utc)
    for i in range(args.iter):
        if i > 0 and args.sleep_between_iter > 0:
            print(f"[iter {i+1}/{args.iter}] sleeping {args.sleep_between_iter}s for TPM window…", flush=True)
            time.sleep(args.sleep_between_iter)
        print(f"[iter {i+1}/{args.iter}] running…", flush=True)
        try:
            per_turn = run_once(args.dataset, args.model)
        except Exception as e:
            traceback.print_exc()
            rows.append({"iter": i + 1, "crashed": str(e), "turns": []})
            continue

        classified = []
        for idx, calls in enumerate(per_turn, start=1):
            v, d = classify_turn(idx, calls, celltypes, var_names)
            classified.append({"turn": idx, "verdict": v, "detail": d})
        rows.append({"iter": i + 1, "crashed": None, "turns": classified})

        for r in classified:
            if r["turn"] in BUG_TURNS:
                print(f"  turn {r['turn']:>2}: {r['verdict']:<7} {r['detail']}", flush=True)

    turn_rates = {
        t: sum(
            1
            for r in rows
            if not r["crashed"]
            and any(x["turn"] == t and x["verdict"] == "BUG" for x in r["turns"])
        )
        for t in BUG_TURNS
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        f.write(f"# T-037 replay report\n\n")
        f.write(f"- Commit: `{commit}`\n")
        f.write(f"- Model: `{args.model}`\n")
        f.write(f"- Dataset: `{args.dataset.name}`\n")
        f.write(f"- Iterations: {args.iter}\n")
        f.write(f"- Started (UTC): {started.isoformat()}\n\n")
        f.write("## Aggregate bug rate\n\n")
        for t, count in turn_rates.items():
            f.write(f"- Turn {t} ({BUG_TURNS[t]}): **{count}/{args.iter}** bug\n")
        f.write("\n## Per-iteration turn detail\n\n")
        f.write("| iter | " + " | ".join(f"t{t}" for t in BUG_TURNS) + " | all turns |\n")
        f.write("|---|" + "---|" * (len(BUG_TURNS) + 1) + "\n")
        for r in rows:
            if r["crashed"]:
                f.write(f"| {r['iter']} | CRASH | CRASH | {r['crashed']} |\n")
                continue
            bug_cells = []
            for t in BUG_TURNS:
                match = next((x for x in r["turns"] if x["turn"] == t), None)
                bug_cells.append(f"{match['verdict']}: {match['detail']}" if match else "-")
            all_tools = "; ".join(
                f"t{x['turn']}={x['detail']}" for x in r["turns"]
            )
            f.write(f"| {r['iter']} | " + " | ".join(bug_cells) + f" | {all_tools} |\n")

    print(f"\nReport: {args.out}")
    for t, count in turn_rates.items():
        print(f"Turn {t}: {count}/{args.iter} bug")
    return 0


if __name__ == "__main__":
    sys.exit(main())
