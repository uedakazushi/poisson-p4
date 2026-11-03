#!/usr/bin/env python3
"""Portable full135 cubic replay and exact comparison to the handoff matrix."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import shutil
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import s61_binary_hyperplane_cubic_matrix_20261102 as producer


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def portable_compile(source, target):
    flags = subprocess.check_output(
        ["pkg-config", "--cflags", "--libs", "flint"], text=True
    ).split()
    subprocess.run(
        ["c++", "-w", "-O3", "-std=c++17", str(source), *flags,
         "-lgmp", "-o", str(target)],
        check=True, capture_output=True, text=True,
    )


def run(output):
    begin = time.time()
    output.mkdir(parents=True, exist_ok=False)
    producer.compile_native = portable_compile
    result = producer.build_matrix(output / "source")
    expected_path = ROOT / "reports/S61_BINARY_HYPERPLANE_CUBIC_MODULAR_20261102.json"
    expected = json.loads(expected_path.read_text())
    for key in ("full_matrix_sha256", "schur_matrix_sha256"):
        assert result[key] == expected["source_matrix"][key]
    archive = ROOT / "research/handoff/20261102/prior_binary_probes/HYPERPLANE_CUBIC_D_DEGREE_20261102.tar.gz"
    member = "native_data/hyperplane-full-cubic-schur36-input.txt"
    with tarfile.open(archive) as tar:
        saved = tar.extractfile(member).read()
    current = (output / "source/schur").read_bytes()
    old_lines = saved.decode().splitlines()
    new_lines = current.decode().splitlines()
    assert new_lines[0] == "42 70 21"
    assert old_lines[0] == "42 70 36"
    assert len(new_lines) == len(old_lines) == 2943
    # The archived native determinant input has no pivot metadata. Its final
    # two lines select a minor; every one of its 2940 matrix entries agrees.
    assert new_lines[3:] == old_lines[1:2941]
    assert old_lines[2941] == "1"
    assert len(old_lines[2942].split()) == 72
    with tarfile.open(output / "matrices.tar.gz", "w:gz") as tar:
        for name in ("full", "schur", "weights", "row_meta"):
            tar.add(output / "source" / name, arcname=name)
    shutil.rmtree(output / "source")
    (output / "source-report.json").write_text(json.dumps(result, indent=2) + "\n")
    report = dict(
        all_passed=True,
        scope="Exact QQ reconstruction of the full135 weighted cubic matrix; no all-core rank or ambient nilpotence assertion.",
        source=result,
        portable_compile_only=True,
        expected_report=str(expected_path.relative_to(ROOT)),
        archived_member=member,
        archived_member_sha256=hashlib.sha256(saved).hexdigest(),
        archived_and_current_matrix_entry_lines_equal=2940,
        archived_header_records_target_rank=36,
        current_header_records_removed_known_unit_pivots=21,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in
                      (Path(__file__), expected_path, archive)},
        generated_sha256={name: sha(output / name) for name in
                          ("source-report.json", "matrices.tar.gz")},
        seconds=time.time() - begin,
    )
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = run(args.output)
    print(json.dumps({key: report[key] for key in
                      ("all_passed", "archived_and_current_matrix_entry_lines_equal", "seconds")}))
