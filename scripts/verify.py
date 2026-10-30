#!/usr/bin/env python3
"""Run the archived exact audits without writing into immutable research sources.

The quick suite is the full v24 audit (1298 checks), not a sample. The rank2
suite reproduces the ten v14 verification groups. Integrity verifies source
archives/materialization and runs the original v24/v23 manifest checker on a
temporary copy with the legacy history symlinks expanded.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]


class VerificationFailure(RuntimeError):
    """A required archived audit or its expected result failed."""


def execute(command: list[str], cwd: Path, log: Path) -> None:
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise VerificationFailure(f"Command exited {result.returncode}; see {log}")


def load_report(path: Path) -> dict:
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise VerificationFailure(f"Missing or invalid result {path}: {exc}") from exc
    if not isinstance(report, dict) or report.get("all_passed") is not True:
        raise VerificationFailure(f"The audit did not report all_passed=true: {path}")
    return report


def require_file(path: Path) -> None:
    if not path.is_file():
        raise VerificationFailure(f"Missing {path}; run bash scripts/setup.sh or python scripts/bootstrap_sources.py first.")


def quick(output: Path, jobs: int) -> dict:
    del jobs  # The immutable v24 runner owns its own four-worker configuration.
    package = ROOT / "legacy/v24"
    script = package / "code/run_checks.py"
    require_file(script)
    report_path = output / "v24_verification.json"
    execute([sys.executable, str(script), "--output", str(report_path)], package, output / "v24.log")
    report = load_report(report_path)
    if report.get("total_checks") != 1298:
        raise VerificationFailure(f"Expected the frozen v24 1298-check audit, received {report.get('total_checks')}.")
    return {"all_passed": True, "total_checks": 1298, "report": str(report_path.relative_to(output.parent))}


def rank2(output: Path, jobs: int) -> dict:
    package = ROOT / "legacy/v14"
    script = package / "code/run_checks.py"
    require_file(script)
    generated = output / "v14_results"
    execute([sys.executable, str(script), "--output", str(generated), "--jobs", str(jobs)], package, output / "v14.log")
    report_path = generated / "verification_summary.json"
    report = load_report(report_path)
    if report.get("number_of_checks") != 10:
        raise VerificationFailure(f"Expected the frozen v14 ten verification groups, received {report.get('number_of_checks')}.")
    return {"all_passed": True, "verification_groups": 10, "report": str(report_path.relative_to(output.parent))}


def integrity(output: Path, jobs: int) -> dict:
    del jobs
    bootstrap = ROOT / "scripts/bootstrap_sources.py"
    require_file(bootstrap)
    require_file(ROOT / "sources/manifest.json")
    execute([sys.executable, str(bootstrap), "--check"], ROOT, output / "source_integrity.log")
    package = ROOT / "legacy/v24"
    require_file(package / "code/verify_bundle.py")
    # pathlib.rglob in the original checker does not recurse into symlinked
    # history. Expand it only in a temporary copy; never patch the old checker.
    with tempfile.TemporaryDirectory(prefix="poisson_p4_integrity_") as temp:
        expanded = Path(temp) / "v24"
        shutil.copytree(package, expanded, symlinks=False, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        execute([sys.executable, str(expanded / "code/verify_bundle.py")], expanded, output / "v24_manifest.log")
    report = load_report(output / "v24_manifest.log")
    return {"all_passed": True, "source_archives_and_materialization": "passed", "v24_manifest": report}


def resolve_output(value: Path | None, suite: str) -> Path:
    if value is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        result = ROOT / "reports" / f"{stamp}-{suite}"
    else:
        result = value if value.is_absolute() else ROOT / value
    result = result.resolve()
    for protected in (ROOT / "legacy", ROOT / "sources", ROOT / ".git"):
        protected = protected.resolve()
        if result == protected or protected in result.parents:
            raise VerificationFailure(f"Report output must not be inside immutable/protected directory {protected}.")
    if result == ROOT:
        raise VerificationFailure("Choose a report subdirectory, for example --output reports/manual.")
    if result.exists() and any(result.iterdir()):
        raise VerificationFailure(f"Report directory is not empty: {result}. Choose a new --output directory to preserve earlier results.")
    result.mkdir(parents=True, exist_ok=True)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--suite", choices=("quick", "rank2", "integrity", "all"), default="quick")
    parser.add_argument("--output", type=Path, help="New report directory; relative paths use the repository root. Default: reports/<UTC timestamp>-<suite>.")
    parser.add_argument("--jobs", type=int, default=4, help="Concurrent v14 verification groups (default: 4).")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    try:
        output = resolve_output(args.output, args.suite)
    except (VerificationFailure, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2

    started = time.monotonic()
    names = ("integrity", "quick", "rank2") if args.suite == "all" else (args.suite,)
    runners = {"quick": quick, "rank2": rank2, "integrity": integrity}
    results = []
    for name in names:
        print(f"RUN {name}", flush=True)
        destination = output / name
        destination.mkdir()
        suite_started = time.monotonic()
        try:
            result = {"suite": name, **runners[name](destination, args.jobs)}
        except (VerificationFailure, OSError, ValueError) as exc:
            result = {"suite": name, "all_passed": False, "error": str(exc)}
        result["seconds"] = round(time.monotonic() - suite_started, 3)
        results.append(result)
        print(f"{'PASS' if result['all_passed'] else 'FAIL'} {name}", flush=True)
        if not result["all_passed"]:
            print(result["error"], file=sys.stderr, flush=True)
            # Never execute research code after an integrity failure.
            if name == "integrity":
                break

    report = {
        "all_passed": len(results) == len(names) and all(item["all_passed"] for item in results),
        "requested_suite": args.suite, "suites": results,
        "python_version": sys.version.split()[0],
        "seconds": round(time.monotonic() - started, 3),
        "scope": "Exact archived identities and rank certificates; this does not replace proofs of geometric exhaustiveness.",
    }
    summary = output / "verification_summary.json"
    summary.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Report: {summary}", flush=True)
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
