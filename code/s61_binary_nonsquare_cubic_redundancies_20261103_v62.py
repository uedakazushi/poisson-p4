#!/usr/bin/env python3
"""Replay three exact known-unit row redundancies of the full135 cubic matrix."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(output):
    begin = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    certificate = ROOT / "code/data/s61_binary_nonsquare_cubic_redundancies_20261103_v62.json"
    native = ROOT / "code/s61_binary_nonsquare_cubic_redundancies_20261103_v62.cpp"
    note = ROOT / "research/notes/S61_BINARY_NONSQUARE_CUBIC_REDUNDANCIES_20261103_v62.md"
    previous = ROOT / "reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62"
    source_report = json.loads((previous / "report.json").read_text())
    data = json.loads(certificate.read_text())
    assert source_report["all_passed"]
    expected = data["original_full_matrix_sha256"]
    assert source_report["source"]["full_matrix_sha256"] == expected
    full = output / "full-input.txt"
    with tarfile.open(previous / "matrices.tar.gz") as archive:
        full.write_bytes(archive.extractfile("full").read())
    assert digest(full) == expected

    (output / "units.txt").write_text("\n".join(data["known_units"]) + "\n")
    with (output / "relations.txt").open("w") as stream:
        stream.write(f"{len(data['relations'])} 63\n")
        for relation in data["relations"]:
            stream.write(f"{relation['deleted_original_row']} {relation['coefficient_degree_bound']}\n")
            coefficients = {item["original_row"]: item["terms"] for item in relation["coefficients"]}
            for row in range(63):
                terms = coefficients.get(row, [])
                stream.write(str(len(terms)) + "".join(" " + " ".join(map(str, term)) for term in terms) + "\n")

    flags = subprocess.check_output(["pkg-config", "--cflags", "--libs", "flint"], text=True).split()
    executable = output / "native"
    command = ["c++", "-O2", "-std=c++17", str(native), *flags, "-o", str(executable)]
    compilation = subprocess.run(command, capture_output=True, text=True, check=True)
    (output / "compile.log").write_text(compilation.stdout + compilation.stderr)
    remaining = output / "remaining.txt"
    process = subprocess.run(
        [str(executable), str(full), str(output / "units.txt"), str(output / "relations.txt"),
         str(output / "normalization.txt"), str(remaining)],
        text=True, capture_output=True, check=True, timeout=60,
    )
    (output / "native.log").write_text(process.stdout + process.stderr)
    assert "EXACT_KNOWN_UNIT_ROW_NORMALIZATION_PASS 63" in process.stdout
    assert "EXACT_CUBIC_REDUNDANCY_PASS 3 REMAINING_ROWS 60" in process.stdout
    for index, relation in enumerate(data["relations"]):
        assert f"EXACT_QQ_ROW_RELATION_PASS {index} DELETED_ORIGINAL_ROW {relation['deleted_original_row']}" in process.stdout
    assert remaining.open().readline().strip() == "60 91 57"
    remaining_hash = digest(remaining)
    with remaining.open("rb") as source, gzip.open(output / "remaining.txt.gz", "wb", compresslevel=9) as target:
        shutil.copyfileobj(source, target)
    with gzip.open(output / "remaining.txt.gz", "rb") as stream:
        assert hashlib.sha256(stream.read()).hexdigest() == remaining_hash
    removed = [item["deleted_original_row"] for item in data["relations"]]
    files = [Path(__file__), native, certificate, note, previous / "report.json", previous / "matrices.tar.gz"]
    report = dict(
        all_passed=True,
        scope="Three exact QQ polynomial syzygies with known-unit coefficients; no all-point cubic rank or ambient-removal theorem.",
        coefficient_field="QQ[a,e,d] localized at a*d*F*B*C*J",
        original_matrix_shape=[63, 91],
        original_full_matrix_sha256=expected,
        deleted_original_zero_based_rows=removed,
        certificate_coefficient_degree_bounds=[item["coefficient_degree_bound"] for item in data["relations"]],
        exact_known_unit_normalization_and_reconstruction_checked=True,
        each_relation_checked_on_all_91_columns=True,
        each_deleted_coefficient_is_a_known_unit=True,
        previous_deleted_coefficients_zero_in_later_relations=True,
        remaining_matrix_shape=[60, 91],
        remaining_original_row_indices=[i for i in range(63) if i not in removed],
        remaining_matrix_sha256=remaining_hash,
        generic_rank_target_from_previous_geometry=57,
        additional_generic_row_dependencies=3,
        all_point_cubic_rank57_claimed=False,
        all_core_quartic_or_ambient_claimed=False,
        discovery_scope="Finite coefficient degree only; the replay verifies identities independently of modular discovery.",
        compile_command=command,
        input_sha256={str(path.relative_to(ROOT)): digest(path) for path in files},
        seconds=time.monotonic() - begin,
    )
    full.unlink()
    remaining.unlink()
    executable.unlink()
    report["generated_sha256"] = {path.name: digest(path) for path in output.iterdir() if path.is_file()}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = run(arguments.output)
    print("PASS", result["deleted_original_zero_based_rows"], result["remaining_matrix_shape"], result["seconds"])
