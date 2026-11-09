#!/usr/bin/env python3
"""Verify the v16 proof certificates in disposable copies.

Run from any directory: python3 code/run_checks.py --output /tmp/v16_checks
Requires Python >=3.11; two auxiliary checks use SymPy 1.14.0.
No Singular installation or Groebner calculation is needed.
Archived code, certificates and historical editions are never overwritten.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
OLD = "previous_v15/code/normal_geometry"
GEOM = "code/geometry"
SEMI = "code/semisimple"
NIL = "code/nilpotent"
CASES = [
    ("v15_original_Jacobi", f"{OLD}/verify_rnc4_slices.py", [],
     f"{OLD}/rnc4_independent_verification.json", None),
    ("conormal_cover_28", f"{GEOM}/verify_geometry.py", [],
     f"{GEOM}/geometry_independent_verification.json", "all_passed"),
    ("semisimple_21", f"{SEMI}/verify_candidate_witnesses.py", [],
     f"{SEMI}/candidate_independent_verification.json", "all_21_passed"),
    ("semisimple_nonnormal", f"{SEMI}/verify_nonnormal_representatives.py", [],
     f"{SEMI}/nonnormal_representatives_independent_verification.json", "all_checks_passed"),
    ("semisimple_fixed_scheme", f"{SEMI}/verify_fixed_h.py", [],
     f"{SEMI}/fixed_h_independent_verification.json", "all_checks_passed"),
    ("semisimple_full_tangents", f"{SEMI}/verify_h_tangent.py", [],
     f"{SEMI}/h_tangent_independent_verification.json", "all_checks_passed"),
    ("nilpotent_7", f"{NIL}/verify_nilpotent_support.py",
     ["--source", f"{OLD}/rnc4_slices.json", "--candidates", f"{GEOM}/candidate_parameterizations.json",
      "--output", "nilpotent_support_checked.json"], "nilpotent_support_checked.json", "all_passed"),
    ("nilpotent_nonnormal", f"{NIL}/verify_geometric_obstruction.py",
     ["--geometry", GEOM, "--output", "nilpotent_geometry_checked.json"],
     "nilpotent_geometry_checked.json", "all_passed"),
    ("nilpotent_fixed_scheme_and_tangents", f"{NIL}/verify_invariant_slice.py",
     ["--source", f"{OLD}/rnc4_slices.json", "--output", "nilpotent_fixed_checked.json"],
     "nilpotent_fixed_checked.json", "all_passed"),
    ("hard_nilpotent_independent", f"{GEOM}/verify_hard_nilpotent.py",
     ["--nilpotent", NIL], f"{GEOM}/hard_nilpotent_independent_verification.json", None),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(case, output):
    name, script, arguments, report_name, flag = case
    started = time.monotonic()
    result = {"name": name, "script": script, "script_sha256": sha(ROOT / script)}
    with tempfile.TemporaryDirectory(prefix="p4_v16_check_") as tmp:
        work = Path(tmp)
        shutil.copytree(ROOT / "code", work / "code",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        # Only the direct v15 inputs are needed, not a mutable history symlink.
        shutil.copytree(ROOT / OLD, work / OLD,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        report_path = work / report_name
        report_path.unlink(missing_ok=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        env.pop("P4_POISSON_V15_DIR", None)
        env.pop("P4_POISSON_V16_GEOMETRY_DIR", None)
        proc = subprocess.run([sys.executable, str(work / script), *arguments],
                              cwd=work, env=env, capture_output=True, text=True)
        (output / f"{name}.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        try:
            if proc.returncode:
                raise RuntimeError(f"verifier exited with status {proc.returncode}")
            report = json.loads(report_path.read_text())
            if flag:
                assert report[flag] is True, f"report flag {flag} is not true"
            if name == "v15_original_Jacobi":
                assert report["rank_exact"] == 24 and report["independent_kernel_vectors"] == 36
                assert len(report["cases"]) == 2
                assert all(c["equivalent_selected_quadrics_verified"] == 94
                           and c["divergence_exact"] and c["restriction_exact"]
                           for c in report["cases"].values())
            if name == "hard_nilpotent_independent":
                assert report["source_rebuilt_from_v15"] and report["singular_source_agrees"]
                assert len(report["targets"]) == 3
            if name == "semisimple_21":
                assert len(report["checks"]) == 21
                assert sum(c.get("ideal_equalities_verified", False) for c in report["checks"]) == 19
            if name == "nilpotent_7":
                assert report["number_of_cases"] == 7
                assert report["number_of_exact_target_identities"] == 23
            if name == "conormal_cover_28":
                assert report["candidate_families_verified"] == 28
                assert report["all_branch_bases_verified_complete"]
            (output / "reports" / f"{name}.json").write_text(
                json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            result["passed"] = True
        except Exception as exc:
            result.update(passed=False, error=str(exc))
    result["seconds"] = round(time.monotonic() - started, 3)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix="p4_v16_results_"))
    (output / "reports").mkdir(parents=True, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = [pool.submit(execute, case, output) for case in CASES]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result['name']}: {'PASS' if result['passed'] else 'FAIL'}", flush=True)
    order = {case[0]: i for i, case in enumerate(CASES)}
    results.sort(key=lambda r: order[r["name"]])
    report = {"all_passed": all(r["passed"] for r in results), "number_of_checks": len(results),
              "python_version": sys.version.split()[0],
              "scope": "Exact identities and full-scheme slice certificates for the v16 normal RNC4 exclusion; not a global component classification.",
              "input_sha256": {
                  "v15_rnc4_slices": sha(ROOT / OLD / "rnc4_slices.json"),
                  "v15_equation_reduction": sha(ROOT / OLD / "rnc4_equation_reduction.json"),
                  "v16_candidates": sha(ROOT / GEOM / "candidate_parameterizations.json")},
              "checks": results}
    (output / "verification_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"all_passed": report["all_passed"], "number_of_checks": len(results),
                      "output": str(output)}, indent=2))
    if not report["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
