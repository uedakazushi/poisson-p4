#!/usr/bin/env python3
"""Same weight4 minor on all four special heavy-curve points; proper QQ lift."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

import numpy as np

from s61_binary_nonsquare_heavy_curve_cubic_20261103_v62 import QuadraticField, modular
from s61_binary_nonsquare_partial135_20261103_v62 import ARCHIVE, EXPECTED, MEMBER

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_report(path):
    report = json.loads(path.read_text())
    assert report["all_passed"]
    for name, expected in report["input_sha256"].items():
        assert sha(ROOT/name) == expected, name
    for name, expected in report.get("generated_sha256", {}).items():
        assert sha(path.parent/name) == expected, name
    return report


def gaussian(field, values):
    """Exact quadratic-field elimination, with vectorized integer operations."""
    A = np.array(values, dtype=np.int64)
    nr, nc, dim = A.shape
    assert dim == 2
    rank, pivots, determinant = 0, [], (1, 0)
    for j in range(nc):
        ids = np.flatnonzero(np.any(A[rank:, j, :], axis=1))
        if not len(ids):
            continue
        i = int(ids[0])+rank
        if i != rank:
            A[[i, rank]] = A[[rank, i]]
            determinant = ((-determinant[0]) % PRIME, (-determinant[1]) % PRIME)
        determinant = field.mul(determinant, tuple(map(int, A[rank, j])))
        inverse = field.inv(tuple(map(int, A[rank, j])))
        x, y = A[rank, :, 0].copy(), A[rank, :, 1].copy()
        A[rank, :, 0] = (x*inverse[0]-field.c*y*inverse[1]) % PRIME
        A[rank, :, 1] = (x*inverse[1]+y*inverse[0]-field.b*y*inverse[1]) % PRIME
        if rank+1 < nr:
            aa = A[rank+1:, j, 0].copy()[:, None]
            bb = A[rank+1:, j, 1].copy()[:, None]
            x, y = A[rank, :, 0].copy()[None, :], A[rank, :, 1].copy()[None, :]
            A[rank+1:, :, 0] = (A[rank+1:, :, 0]-aa*x+field.c*bb*y) % PRIME
            A[rank+1:, :, 1] = (A[rank+1:, :, 1]-aa*y-bb*x+field.b*bb*y) % PRIME
        pivots.append(j)
        rank += 1
        if rank == nr:
            break
    return rank, pivots, determinant if nr == nc and rank == nr else None


def run(output):
    began = time.time()
    output.mkdir(parents=True, exist_ok=False)
    heavy_path = ROOT/"reports/S61_BINARY_NONSQUARE_HEAVY_CURVE_CUBIC_20261103_v62/report.json"
    partial_path = ROOT/"reports/S61_BINARY_NONSQUARE_PARTIAL135_20261103_v62/report.json"
    heavy = check_report(heavy_path)
    partial = check_report(partial_path)
    assert heavy["special_generic_a_algebra_length"] == 4
    assert heavy["special_projective_infinity_dimension"] == 0
    assert heavy["all_vertical_exception_fibers_zero_dimensional"]
    with tarfile.open(ARCHIVE) as archive:
        raw = archive.extractfile(MEMBER).read()
    assert hashlib.sha256(raw).hexdigest() == EXPECTED
    old = json.loads(raw)
    additional_path = partial_path.parent/"additional-G6.json"
    additional = json.loads(additional_path.read_text())
    records = old["equations_sparse"] + [dict(original_row=134, weight=3, polynomial=additional["G6"])]
    variables = old["remaining_variables"]
    weights = dict(zip(variables, additional["remaining_weights"]))
    assert len(variables) == 10 and sorted(weights.values()) == [1]*7+[2]*3

    def monomials(weight, start=0):
        if weight == 0:
            yield ()
            return
        for k in range(start, len(variables)):
            v = variables[k]
            if weights[v] <= weight:
                for mon in monomials(weight-weights[v], k):
                    yield (v,)+mon

    mons = {w: list(monomials(w)) for w in range(5)}
    columns = mons[4]
    index = {mon: j for j, mon in enumerate(columns)}
    labels = [(i, shift) for i, q in enumerate(records) for shift in mons[4-q["weight"]]]
    assert len(columns) == 300 and len(labels) == 428
    common_rows, points = None, []
    m2 = shutil.which("M2")
    assert m2
    for number, point in enumerate(heavy["quadratic_Galois_pairs"]):
        b, c = point["quadratic"]
        E, D = point["e"], point["d"]
        field = QuadraticField(b, c)
        maxe = max(ex[1] for q in records for _, cs in q["polynomial"] for ex, _ in cs)
        maxd = max(ex[2] for q in records for _, cs in q["polynomial"] for ex, _ in cs)
        ep = [field.pow(E, k) for k in range(maxe+1)]
        dp = [field.pow(D, k) for k in range(maxd+1)]
        evaluated = []
        for q in records:
            row = {}
            for mon, coefficients in q["polynomial"]:
                value = (0, 0)
                assert sum(weights[j] for j in mon) == q["weight"]
                for ex, coefficient in coefficients:
                    scalar = modular(coefficient)
                    term = field.mul(ep[ex[1]], dp[ex[2]])
                    value = field.add(value, (term[0]*scalar % PRIME, term[1]*scalar % PRIME))
                if value != (0, 0):
                    row[tuple(mon)] = value
            evaluated.append(row)
        matrix = np.zeros((len(labels), len(columns), 2), dtype=np.int64)
        for k, (i, shift) in enumerate(labels):
            for mon, value in evaluated[i].items():
                matrix[k, index[tuple(sorted(mon+shift))]] = value
        if common_rows is None:
            rank, common_rows, _ = gaussian(field, matrix.transpose(1, 0, 2))
            assert rank == 300
        rank, _, determinant = gaussian(field, matrix[common_rows])
        assert rank == 300 and determinant != (0, 0)
        expressions = []
        for row in evaluated:
            expressions.append("+".join("("+str(v[0])+"+"+str(v[1])+"*kk62)"+
                                          "".join("*q"+str(j) for j in mon)
                                          for mon, v in row.items()) or "0")
        source = "\n".join([
            "F=toField((ZZ/32003)[kk62]/ideal(kk62^2+"+str(b)+"*kk62+"+str(c)+"));R=F["+
            ",".join("q"+str(v) for v in variables)+",Degrees=>{"+
            ",".join(str(weights[v]) for v in variables)+"}];",
            "I=ideal("+",".join(expressions)+");G=gb I;",
            "assert(dim(R/I)==0);assert(numColumns basis(R/I)==65);",
            "assert(apply(toList(0..5),n->hilbertFunction(n,R/I))=={1,7,23,34,0,0});",
            'print "FULL_PARTIAL135_TRANSVERSE_LENGTH65_PASS";',
        ])+"\n"
        script = output/("point"+str(number)+".m2")
        script.write_text(source)
        cp = subprocess.run([m2, "--script", str(script)], cwd=ROOT, text=True,
                            capture_output=True, timeout=120)
        (output/("point"+str(number)+".log")).write_text(cp.stdout+cp.stderr)
        cp.check_returncode()
        assert "FULL_PARTIAL135_TRANSVERSE_LENGTH65_PASS" in cp.stdout
        points.append(dict(quadratic=[b, c], e=E, d=D, common_quartic_minor_determinant=determinant,
                           rank=rank, transverse_length=65, weighted_Hilbert_function=[1,7,23,34,0,0]))
    degree = max(sum(ex) for k in common_rows for _, cs in records[labels[k][0]]["polynomial"] for ex, _ in cs)
    paths = [Path(__file__), heavy_path, partial_path, additional_path, ARCHIVE,
             ROOT/"code/s61_binary_nonsquare_heavy_curve_cubic_20261103_v62.py"]
    report = dict(
        all_passed=True, prime=PRIME, partial_source_includes_G6=True,
        variables=variables, variable_weights=[weights[v] for v in variables], cuts=[5,19,23],
        weight4_matrix_shape=[428,300], weight4_column_monomials=columns,
        common_minor_row_indices=common_rows,
        common_minor_source_labels=[[records[labels[k][0]]["original_row"], list(labels[k][1])] for k in common_rows],
        common_entry_homogenization_degree=degree, homogeneous_determinant_degree=300*degree,
        quadratic_Galois_pairs=points,
        conclusion="On every QQ heavy-defect curve, the partial135 weight4 matrix is generically surjective. Together with the quadratic/cubic bounds the actual transverse length is <=65 generically.",
        arithmetic_bridge="The same integral-ann ordinaryP4 proper model and finitefree4 special cover as the preceding cubic theorem; entrywise fixed-degree homogenized quartic determinant.",
        limitations=["Does not prove all other nonsquare defect curves are removed.",
                     "Ambient H11 removal uses the previously proved relative-free65 and nuisance-contraction theorem, not fixed-fiber nilpotence alone."],
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},
        generated_sha256={path.name:sha(path) for path in output.iterdir()},
        seconds=time.time()-began,
    )
    (output/"report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({key:report[key] for key in ("all_passed","conclusion","seconds")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    run(parser.parse_args().output)
