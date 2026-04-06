#!/usr/bin/env python3
"""
In-process regression test runner for nvwa-mvp scRNA visualization.

Loads each test case from tests/tests.yaml, invokes the agent directly
(no HTTP, no server), checks must_contain / must_not_contain rules against
the final response text, validates actual artifacts (PNG bytes, CSV rows),
and writes a markdown report.

Usage
-----
  python scripts/run_tests.py --dataset /abs/path/to/pbmc3k.h5ad
  python scripts/run_tests.py --dataset /abs/path/to/pbmc3k.h5ad --limit 3
  python scripts/run_tests.py --dataset /abs/path --model gpt-4o --report out.md

Environment
-----------
  OPENAI_API_KEY   Required.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

# Suppress noisy-but-harmless warnings from scanpy internals and anndata.
# PerformanceWarning fires during rank_genes_groups on fragmented DataFrames;
# it is a scanpy internal issue and does not affect correctness.
warnings.filterwarnings("ignore", message="DataFrame is highly fragmented")
warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")

# ── Paths ──────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent

# Add repo root so `src.*` imports resolve regardless of cwd.
sys.path.insert(0, str(REPO_ROOT))


def _check_deps() -> None:
    """Fail fast with a helpful message if project deps are not installed."""
    missing = []
    for pkg in ("langchain_core", "langchain_openai", "anndata", "scanpy"):
        try:
            __import__(pkg)
        except ModuleNotFoundError:
            missing.append(pkg)
    if missing:
        print(
            "ERROR: Missing packages: " + ", ".join(missing) + "\n"
            "       Run:  uv sync\n"
            "       Then: uv run python scripts/run_tests.py --dataset ...\n"
            "       (or activate the project venv before calling python directly)",
            file=sys.stderr,
        )
        sys.exit(2)


DEFAULT_CONFIG = REPO_ROOT / "tests" / "integration" / "tests.yaml"
DEFAULT_REPORT_DIR = REPO_ROOT / "local" / "reports" / "regression"
DEFAULT_MODEL = "gpt-4o-mini"

# PNG file signature — first 4 bytes of any valid PNG
PNG_MAGIC = b"\x89PNG"
# A real matplotlib plot at 150 dpi is always well above 5 KB;
# anything smaller is a blank/empty figure.
MIN_PLOT_BYTES = 5_000


# ── Data model ─────────────────────────────────────────────────────────────────

@dataclass
class CaseResult:
    case_id: str
    category: str
    prompt: str
    expected_status: str          # PASS | FAIL | WARN
    status: str = "pending"       # pass | fail | warn | error | skip
    duration_s: float = 0.0
    final_text: str = ""
    tool_called: bool = False
    plot_count: int = 0           # number of valid PlotResults produced
    table_count: int = 0          # number of valid TableResults produced
    failures: list[str] = field(default_factory=list)
    notes: str = ""               # propagated from YAML + runner notes
    failure_category: str = ""    # gene_not_found | missing_key | out_of_scope | …


# ── Artifact validators ────────────────────────────────────────────────────────

def _valid_plot(p: Any) -> bool:
    """Return True if PlotResult contains a real PNG image (not blank/empty)."""
    return (
        isinstance(p.image, (bytes, bytearray))
        and p.image[:4] == PNG_MAGIC
        and len(p.image) >= MIN_PLOT_BYTES
    )


def _valid_table(t: Any) -> bool:
    """Return True if TableResult contains at least one data row beyond the header."""
    lines = t.csv_data.strip().splitlines()
    return len(lines) >= 2


# ── Dataset loader ─────────────────────────────────────────────────────────────

def load_adata(dataset_path: Path) -> Any:
    """Load and return an AnnData, suppressing legacy format warnings."""
    import anndata as ad
    return ad.read_h5ad(dataset_path)


# ── Single-test runner ─────────────────────────────────────────────────────────

def run_one(
    tc: dict[str, Any],
    dataset_path: Path,
    api_key: str,
    model: str,
) -> tuple[str, bool, list, list]:
    """
    Load a fresh AnnData, create the agent, invoke the prompt, and collect artifacts.

    Returns (final_text, tool_called, plot_results, table_results).
    Raises on hard failures (caller catches and marks as error).
    """
    from src.agent.core import create_agent
    from src.agent.tools import clear_plot_results, clear_table_results, get_plot_results, get_table_results
    from src.types import detect_dataset_state

    # Clear artifact buffers from any previous test before loading the agent.
    clear_plot_results()
    clear_table_results()

    adata = load_adata(dataset_path)
    state = detect_dataset_state(
        adata,
        source=str(dataset_path),
        filename=dataset_path.name,
    )
    agent = create_agent(adata, api_key=api_key, model=model, dataset_state=state)
    resp = agent.invoke(tc["prompt"])

    # Collect artifacts from the module-level buffers in tools.py.
    # get_plot_results() / get_table_results() clear the buffer after reading,
    # so each test starts with a clean slate.
    plot_results = get_plot_results()
    table_results = get_table_results()

    return resp.text, resp.tool_called, plot_results, table_results


# ── Assertion check ────────────────────────────────────────────────────────────

def check(
    tc: dict[str, Any],
    text: str,
    plot_results: list,
    table_results: list,
) -> list[str]:
    """Return a list of failure strings; empty list → all assertions passed."""
    lower = text.lower()
    failures: list[str] = []

    # Text-based assertions
    for phrase in tc.get("must_contain", []):
        if phrase.lower() not in lower:
            failures.append(f"must_contain: '{phrase}' not in response")

    for phrase in tc.get("must_not_contain", []):
        if phrase.lower() in lower:
            failures.append(f"must_not_contain: forbidden string '{phrase}' found")

    # Artifact-based assertions
    if tc.get("requires_plot"):
        valid = [p for p in plot_results if _valid_plot(p)]
        if not valid:
            n = len(plot_results)
            detail = f"{n} produced but all failed validation (empty/blank)" if n else "none produced"
            failures.append(f"requires_plot: no valid PlotResult ({detail})")

    if tc.get("requires_table"):
        valid = [t for t in table_results if _valid_table(t)]
        if not valid:
            n = len(table_results)
            detail = f"{n} produced but all failed validation (empty/no rows)" if n else "none produced"
            failures.append(f"requires_table: no valid TableResult ({detail})")

    if "expected_artifact_count" in tc:
        expected = tc["expected_artifact_count"]
        actual = len(plot_results) + len(table_results)
        if actual != expected:
            failures.append(
                f"expected_artifact_count: expected {expected}, got {actual} "
                f"({len(plot_results)} plot(s), {len(table_results)} table(s))"
            )

    return failures


# ── Failure category classifier ────────────────────────────────────────────────

def classify_failure(result: CaseResult) -> str:
    """Map a failing result to a broad error category for the summary section."""
    if result.status == "error":
        return "crash"
    if result.status == "skip":
        return "skipped"
    combined = " ".join(result.failures).lower() + result.final_text.lower()
    if "requires_plot" in combined or "requires_table" in combined:
        return "missing_artifact"
    if "gene" in combined and "not found" in combined:
        return "gene_not_found"
    if ("metadata" in combined or "obs" in combined or "key" in combined) and "not found" in combined:
        return "missing_metadata_key"
    if result.category == "out_of_scope":
        return "out_of_scope"
    if "must_contain" in combined:
        return "assertion_failed"
    return "other"


# ── Markdown report ────────────────────────────────────────────────────────────

def _artifact_label(r: CaseResult) -> str:
    parts = []
    if r.plot_count:
        parts.append(f"{r.plot_count} plot{'s' if r.plot_count > 1 else ''}")
    if r.table_count:
        parts.append(f"{r.table_count} table{'s' if r.table_count > 1 else ''}")
    return ", ".join(parts) if parts else "none"


def write_report(
    results: list[CaseResult],
    output_path: Path,
    dataset_path: Path,
    model: str,
    elapsed_s: float,
) -> None:
    n_pass  = sum(1 for r in results if r.status == "pass")
    n_fail  = sum(1 for r in results if r.status == "fail")
    n_warn  = sum(1 for r in results if r.status == "warn")
    n_error = sum(1 for r in results if r.status == "error")
    n_skip  = sum(1 for r in results if r.status == "skip")

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines: list[str] = [
        "# nvwa-mvp Regression Report",
        "",
        f"**Generated:** {now}  ",
        f"**Dataset:** `{dataset_path}`  ",
        f"**Model:** `{model}`  ",
        f"**Duration:** {elapsed_s:.1f}s  ",
        "",
        "---",
        "",
        "## Summary",
        "",
        "| Result | Count |",
        "|--------|-------|",
        f"| ✅ Pass  | {n_pass} |",
        f"| ❌ Fail  | {n_fail} |",
        f"| ⚠️ Warn  | {n_warn} |",
        f"| 💥 Error | {n_error} |",
        f"| ⏭ Skip  | {n_skip} |",
        f"| **Total** | **{len(results)}** |",
        "",
        "---",
        "",
        "## Per-case Results",
        "",
        "| case_id | category | status | duration | artifacts | failure reason |",
        "|---------|----------|--------|----------|-----------|----------------|",
    ]

    for r in results:
        icon = {"pass": "✅", "fail": "❌", "warn": "⚠️", "error": "💥", "skip": "⏭"}.get(r.status, "?")
        reason = "; ".join(r.failures[:2]) if r.failures else (r.notes[:60] if r.status == "skip" else "—")
        lines.append(
            f"| {r.case_id} | {r.category} | {icon} {r.status.upper()} "
            f"| {r.duration_s:.1f}s | {_artifact_label(r)} "
            f"| {reason} |"
        )

    lines += [
        "",
        "---",
        "",
        "## Top Failure Reasons",
        "",
    ]

    failed = [r for r in results if r.status in ("fail", "error")]
    if failed:
        from collections import Counter
        cats = Counter(r.failure_category for r in failed)
        lines += [
            "| Category | Count |",
            "|----------|-------|",
        ]
        for cat, count in cats.most_common():
            lines.append(f"| `{cat}` | {count} |")
    else:
        lines.append("_No failures recorded._")

    lines += [
        "",
        "---",
        "",
        "## Detailed Results",
        "",
    ]

    for r in results:
        icon = {"pass": "✅", "fail": "❌", "warn": "⚠️", "error": "💥", "skip": "⏭"}.get(r.status, "?")
        lines += [
            f"### {icon} {r.case_id}",
            "",
            f"- **Category:** {r.category}",
            f"- **Expected status:** {r.expected_status}",
            f"- **Actual status:** {r.status.upper()}",
            f"- **Duration:** {r.duration_s:.2f}s",
            f"- **Artifacts:** {_artifact_label(r)}",
        ]
        if r.notes:
            lines.append(f"- **Notes:** {r.notes}")
        lines.append("")

        if r.failures:
            lines.append("**Failures:**")
            lines.append("")
            for f in r.failures:
                lines.append(f"- {f}")
            lines.append("")

        if r.final_text:
            preview = r.final_text[:400].replace("\n", " ").replace("|", "\\|")
            lines.append(f"**Response (truncated):** {preview}")
            lines.append("")

        lines.append("---")
        lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run nvwa-mvp in-process regression tests.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--dataset", required=True,
        help="Absolute path to the .h5ad dataset file.",
    )
    parser.add_argument(
        "--config", default=str(DEFAULT_CONFIG),
        help=f"Path to test YAML (default: {DEFAULT_CONFIG}).",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Run only the first N test cases.",
    )
    parser.add_argument(
        "--cases",
        default=None,
        help="Comma-separated list of case_ids to run, e.g. HP-04,IE-02,OOS-01",
    )
    parser.add_argument(
        "--model", default=DEFAULT_MODEL,
        help=f"OpenAI model name (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--report", default=None,
        help=(
            "Output path for the markdown report. "
            "Defaults to reports/regression_report_YYYYMMDD_HHMMSS.md "
            "(timestamped so each run is preserved)."
        ),
    )
    parser.add_argument(
        "--api-key", default=os.environ.get("OPENAI_API_KEY", ""),
        help="OpenAI API key (default: OPENAI_API_KEY env var).",
    )
    args = parser.parse_args()

    # ── Validate inputs ────────────────────────────────────────────────────────

    _check_deps()

    if not args.api_key:
        print("ERROR: OPENAI_API_KEY is not set. Pass --api-key or export the env var.",
              file=sys.stderr)
        return 2

    dataset_path = Path(args.dataset)
    if not dataset_path.exists():
        print(f"ERROR: Dataset not found: {dataset_path}", file=sys.stderr)
        return 2

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: Config not found: {config_path}", file=sys.stderr)
        return 2

    with config_path.open() as f:
        config = yaml.safe_load(f)

    test_cases: list[dict[str, Any]] = config.get("tests", [])
    if args.cases:
        allowed = {c.strip() for c in args.cases.split(",")}
        test_cases = [t for t in test_cases if t["case_id"] in allowed]
    elif args.limit:
        test_cases = test_cases[: args.limit]

    if not test_cases:
        print("No test cases found.", file=sys.stderr)
        return 1

    if args.report:
        report_path = Path(args.report)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = DEFAULT_REPORT_DIR / f"regression_report_{ts}.md"

    # ── Banner ─────────────────────────────────────────────────────────────────

    print(f"\n{'='*62}")
    print(f"  nvwa-mvp Regression Suite")
    print(f"  dataset : {dataset_path}")
    print(f"  model   : {args.model}")
    print(f"  cases   : {len(test_cases)}")
    print(f"  report  : {report_path}")
    print(f"{'='*62}\n")

    results: list[CaseResult] = []
    suite_start = time.monotonic()

    for tc in test_cases:
        case_id         = tc["case_id"]
        category        = tc.get("category", "unknown")
        expected_status = tc.get("expected_status", "PASS").upper()
        tc_notes        = tc.get("notes", "")

        # Manual / SKIP cases (WARN with no must_contain and "MANUAL" in prompt)
        is_manual = (
            expected_status == "WARN"
            and not tc.get("must_contain")
            and "MANUAL" in tc.get("prompt", "").upper()
        )

        print(f"[{case_id}] {tc.get('prompt', '')[:60]!r:<62} ", end="", flush=True)

        result = CaseResult(
            case_id=case_id,
            category=category,
            prompt=tc.get("prompt", ""),
            expected_status=expected_status,
            notes=tc_notes,
        )

        if is_manual:
            result.status = "skip"
            results.append(result)
            print("SKIP (manual)")
            continue

        # Pre-flight: dataset capability checks (cheap, no agent invocation).
        if tc.get("requires_precomputed_umap"):
            import anndata as _ad
            _adata_check = _ad.read_h5ad(dataset_path)
            if "X_umap" not in _adata_check.obsm:
                result.status = "fail"
                result.failures = [
                    "requires_precomputed_umap: X_umap not found in dataset .obsm; "
                    "this test requires a preprocessed dataset"
                ]
                result.failure_category = "missing_artifact"
                results.append(result)
                print(f"FAIL (0.0s) [none]")
                print(f"         ✗ {result.failures[0]}")
                continue

        t0 = time.monotonic()
        try:
            final_text, tool_called, plot_results, table_results = run_one(
                tc, dataset_path, args.api_key, args.model
            )
            result.duration_s = time.monotonic() - t0
            result.final_text = final_text
            result.tool_called = tool_called
            result.plot_count = sum(1 for p in plot_results if _valid_plot(p))
            result.table_count = sum(1 for t in table_results if _valid_table(t))

            failures = check(tc, final_text, plot_results, table_results)
            result.failures = failures

            if not failures:
                result.status = "pass"
            elif expected_status == "WARN":
                result.status = "warn"
            else:
                result.status = "fail"

        except Exception:
            result.duration_s = time.monotonic() - t0
            result.status = "error"
            tb = traceback.format_exc()
            result.failures = [tb.strip().splitlines()[-1]]
            result.final_text = tb

        result.failure_category = classify_failure(result)
        results.append(result)

        icon = {"pass": "PASS", "fail": "FAIL", "warn": "WARN",
                "error": "ERROR", "skip": "SKIP"}.get(result.status, result.status)
        artifact_info = _artifact_label(result)
        print(f"{icon} ({result.duration_s:.1f}s) [{artifact_info}]")
        for f in result.failures[:2]:
            print(f"         ✗ {f}")

    # ── Console summary ────────────────────────────────────────────────────────

    elapsed = time.monotonic() - suite_start
    n_pass  = sum(1 for r in results if r.status == "pass")
    n_fail  = sum(1 for r in results if r.status == "fail")
    n_warn  = sum(1 for r in results if r.status == "warn")
    n_error = sum(1 for r in results if r.status == "error")
    n_skip  = sum(1 for r in results if r.status == "skip")

    print(f"\n{'='*62}")
    print(f"  PASS={n_pass}  FAIL={n_fail}  WARN={n_warn}  "
          f"ERROR={n_error}  SKIP={n_skip}  ({elapsed:.1f}s)")
    print(f"{'='*62}")

    write_report(results, report_path, dataset_path, args.model, elapsed)
    print(f"Report → {report_path}")

    return 0 if (n_fail == 0 and n_error == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
