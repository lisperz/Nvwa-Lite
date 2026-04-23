"""One-shot import map of src/.

Walks src/, parses each .py via ast, and emits a markdown report of internal
module dependencies. Run with no args; output goes to local/reports/import_map.md.
"""

from __future__ import annotations

import ast
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src"
OUT = REPO / "local" / "reports" / "import_map.md"


def module_name(path: Path) -> str:
    rel = path.relative_to(REPO).with_suffix("")
    parts = rel.parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def resolve_relative(module: str | None, level: int, current: str) -> str | None:
    if level == 0:
        return module
    base_parts = current.split(".")[:-level]
    if module:
        base_parts.append(module)
    return ".".join(base_parts) if base_parts else None


def extract_imports(path: Path, current: str) -> list[str]:
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError:
        return []
    imports: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            resolved = resolve_relative(node.module, node.level, current)
            if resolved:
                imports.append(resolved)
    return imports


def is_internal(mod: str) -> bool:
    return mod.startswith("src.") or mod == "src"


def top_level(mod: str) -> str:
    parts = mod.split(".")
    return ".".join(parts[:2]) if len(parts) >= 2 else mod


def main() -> int:
    if not SRC.exists():
        print(f"src/ not found at {SRC}", file=sys.stderr)
        return 1

    files = sorted(p for p in SRC.rglob("*.py") if "__pycache__" not in p.parts)

    per_file: dict[str, list[str]] = {}
    top_edges: set[tuple[str, str]] = set()
    external: set[str] = set()

    for path in files:
        mod = module_name(path)
        imports = extract_imports(path, mod)
        internal = sorted({i for i in imports if is_internal(i)})
        per_file[mod] = internal
        for imp in imports:
            if is_internal(imp):
                src_tl, dst_tl = top_level(mod), top_level(imp)
                if src_tl != dst_tl:
                    top_edges.add((src_tl, dst_tl))
            else:
                external.add(imp.split(".")[0])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Import Map — src/\n")
    lines.append(f"Files scanned: {len(files)}\n")

    lines.append("## Top-level module graph\n")
    lines.append("```mermaid")
    lines.append("graph LR")
    for src_mod, dst_mod in sorted(top_edges):
        a = src_mod.replace(".", "_")
        b = dst_mod.replace(".", "_")
        lines.append(f"  {a}[{src_mod}] --> {b}[{dst_mod}]")
    lines.append("```\n")

    lines.append("## Per-file internal imports\n")
    by_pkg: dict[str, list[str]] = defaultdict(list)
    for mod in sorted(per_file):
        by_pkg[top_level(mod)].append(mod)
    for pkg in sorted(by_pkg):
        lines.append(f"### {pkg}\n")
        for mod in by_pkg[pkg]:
            deps = per_file[mod]
            if deps:
                lines.append(f"- **{mod}**")
                for d in deps:
                    lines.append(f"  - → {d}")
            else:
                lines.append(f"- **{mod}** _(no internal imports)_")
        lines.append("")

    lines.append("## External packages referenced\n")
    for pkg in sorted(external):
        lines.append(f"- {pkg}")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT.relative_to(REPO)} ({len(files)} files, {len(top_edges)} cross-module edges)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
