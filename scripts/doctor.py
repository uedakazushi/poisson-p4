#!/usr/bin/env python3
"""Report local research prerequisites without running certificate computations."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_DEPENDENCIES = {"sympy": "1.14.0", "numpy": "2.3.5", "mpmath": "1.3.0"}


def command_output(command: list[str]) -> str | None:
    try:
        result = subprocess.run(
            command, cwd=ROOT, capture_output=True, text=True, timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def inspect_environment() -> dict:
    dependencies = {}
    for name, expected in EXPECTED_DEPENDENCIES.items():
        try:
            installed = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            installed = None
        dependencies[name] = {
            "expected": expected, "installed": installed,
            "matches": installed == expected,
        }

    paths = {
        "manifest": ROOT / "sources/manifest.json",
        "v24_checks": ROOT / "legacy/v24/code/run_checks.py",
        "v24_integrity": ROOT / "legacy/v24/code/verify_bundle.py",
        "v14_checks": ROOT / "legacy/v14/code/run_checks.py",
    }
    sources = {name: path.is_file() for name, path in paths.items()}
    sources["archive_count"] = len(list((ROOT / "sources").glob("v[0-9][0-9].zip")))
    sources["materialized_version_count"] = sum(
        path.is_dir() for path in (ROOT / "legacy").glob("v[0-9][0-9]")
    )

    latex = {name: shutil.which(name) for name in ("lualatex", "latexmk", "kpsewhich", "pdftoppm")}
    latex["packages"] = {
        name: command_output(["kpsewhich", name]) if latex["kpsewhich"] else None
        for name in ("ltjsarticle.cls", "luaotfload-main.lua", "HaranoAjiMincho-Regular.otf")
    }
    latex["ready_for_japanese_lualatex"] = bool(
        latex["lualatex"] and all(latex["packages"].values())
    )

    # Never print remotes, credential-bearing URLs, environment contents, or diffs.
    git_available = shutil.which("git") is not None
    git = {"available": git_available}
    if git_available:
        status = command_output(["git", "status", "--porcelain=v1", "-z"])
        git.update(
            branch=command_output(["git", "symbolic-ref", "--quiet", "--short", "HEAD"]),
            commit=command_output(["git", "rev-parse", "--verify", "HEAD"]),
            clean=(not status) if status is not None else None,
        )

    ready = (
        sys.version_info >= (3, 11)
        and all(item["matches"] for item in dependencies.values())
        and all(sources[name] for name in paths)
    )
    return {
        "repository": str(ROOT),
        "python": {"version": sys.version.split()[0], "executable": sys.executable,
                   "supported": sys.version_info >= (3, 11)},
        "dependencies": dependencies, "sources": sources, "latex": latex,
        "git": git, "ready_for_research": ready,
        "scope": "Availability check only; no certificate computation or source hash audit.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Print structured diagnostic information.")
    parser.add_argument("--strict", action="store_true", help="Exit 1 if Python dependencies or research sources are missing.")
    args = parser.parse_args()
    report = inspect_environment()
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        print(f"Repository: {ROOT}")
        print(f"Python: {report['python']['version']} ({sys.executable})")
        for name, item in report["dependencies"].items():
            print(f"  {name}: {item['installed'] or 'missing'}; expected {item['expected']}")
        print(f"Source archives: {report['sources']['archive_count']}")
        print(f"Materialized versions: {report['sources']['materialized_version_count']}")
        for name in ("manifest", "v24_checks", "v24_integrity", "v14_checks"):
            print(f"  {name}: {'present' if report['sources'][name] else 'missing'}")
        print("Japanese LuaLaTeX: " + ("available" if report["latex"]["ready_for_japanese_lualatex"] else "incomplete (optional for certificate checks)"))
        for name, path in report["latex"]["packages"].items():
            print(f"  {name}: {'present' if path else 'missing'}")
        git = report["git"]
        print(f"Git: branch={git.get('branch') or '(unavailable)'}; commit={git.get('commit') or '(no commit)'}; clean={git.get('clean')}")
        print("Research setup: " + ("ready" if report["ready_for_research"] else "incomplete; run bash scripts/setup.sh"))
        print(report["scope"])
    return 1 if args.strict and not report["ready_for_research"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
