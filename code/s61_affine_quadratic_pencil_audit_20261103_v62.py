#!/usr/bin/env python3
"""Exact scope audit of constant-column preprocessing for all134 quadrics."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from s61_affine_normal_projection_20261102 import P, pivot_minor

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    previous = ROOT / 'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62'
    saved = json.loads((previous / 'report.json').read_text())
    assert saved['all_passed']
    for name, value in saved['input_sha256'].items():
        assert sha(ROOT / name) == value
    for name, value in saved['evidence_sha256'].items():
        assert sha(previous / name) == value
    with gzip.open(previous / 'restricted-original134.json.gz', 'rt') as stream:
        encoded = json.load(stream)
    assert len(encoded) == 134
    monomials = list(itertools.combinations_with_replacement(range(45), 2))
    positions = {mon: j for j, mon in enumerate(monomials)}
    exact = defaultdict(Fraction)
    for i, row in enumerate(encoded):
        for monomial, coefficients in row:
            if len(monomial) != 2:
                continue
            j = positions[tuple(monomial)]
            for exponent, value in coefficients:
                assert len(exponent) == 5 and sum(exponent) <= 1
                k = 0 if not sum(exponent) else exponent.index(1) + 1
                exact[k, i, j] += Fraction(value)
    C = np.zeros((6, 134, len(monomials)), dtype=np.int64)
    for (k, i, j), q in exact.items():
        assert q.denominator % P
        C[k, i, j] = q.numerator * pow(q.denominator, -1, P) % P
    variable = {j for (k, i, j), q in exact.items() if k and q}
    constant_columns = [j for j in range(len(monomials)) if j not in variable]
    assert len(constant_columns) == 25
    A = sp.Matrix(134, 25, lambda i, j: sp.Rational(exact[0, i, constant_columns[j]]))
    _, basis = A.rref()
    assert len(basis) == 16
    B = A[:, list(basis)]
    row_basis = list(B.T.rref()[1])
    Q = B[row_basis, :].inv() * A[row_basis, :]
    assert B * Q == A
    assert len(pivot_minor(C[0][:, constant_columns])[0]) == 16
    coefficient_ranks = []
    for k in range(1, 6):
        nonzero_rows = np.nonzero(np.any(C[k], axis=1))[0]
        rank = len(pivot_minor(C[k])[0])
        assert rank == len(nonzero_rows) == 64
        coefficient_ranks.append(rank)
    rows, columns, pivot_product = pivot_minor(C[0])
    assert len(rows) == len(columns) == 134 and pivot_product
    M = C[0][:, columns].copy()
    U = np.eye(134, dtype=np.int64)
    for j in range(134):
        k = j + np.nonzero(M[j:, j])[0][0]
        M[[j, k]] = M[[k, j]]
        U[[j, k]] = U[[k, j]]
        inv = pow(int(M[j, j]), -1, P)
        M[j] = M[j] * inv % P
        U[j] = U[j] * inv % P
        factors = M[:, j].copy()
        factors[j] = 0
        M = (M - factors[:, None] * M[j]) % P
        U = (U - factors[:, None] * U[j]) % P
    assert np.array_equal(M, np.eye(134, dtype=np.int64))
    N = np.array([U @ C[k][:, columns] % P for k in range(1, 6)])
    traces = [int(np.trace(matrix) % P) for matrix in N]
    assert all(traces)
    W = np.eye(134, dtype=np.int64)
    image_ranks = [134]
    for _ in range(134):
        image = np.concatenate([matrix @ W % P for matrix in N], axis=1)
        _, chosen, _ = pivot_minor(image)
        image_ranks.append(len(chosen))
        if len(chosen) == W.shape[1]:
            break
        W = image[:, chosen]
    assert image_ranks == [134, 64, 61, 61]
    np.savez_compressed(args.output / 'coefficient-pencils-modp.npz', C=C, N=N,
                        pivot_columns=columns, inverse=U,
                        quadratic_monomials=np.array(monomials), constant_columns=constant_columns)
    proof = dict(constant_column_monomials=[monomials[j] for j in constant_columns],
                 basis_columns=list(basis), row_basis=row_basis,
                 rational_span_coefficients=[[str(q) for q in row] for row in Q.tolist()])
    (args.output / 'constant-column-exact-span.json').write_text(json.dumps(proof, indent=2) + '\n')
    work = ROOT / 'work/affine-infinity-lift-v62'
    for suffix in ['sing', 'log', '-exit.json']:
        name = 'q32-coordinate-unit' + (suffix if suffix.startswith('-') else '.' + suffix)
        path = work / name
        (args.output / (name + '.gz')).write_bytes(gzip.compress(path.read_bytes(), mtime=0))
    receipt = json.loads((work / 'q32-coordinate-unit-exit.json').read_text())
    assert receipt['exit_code'] == 14 and receipt['virtual_memory_KiB'] == 3145728
    log = (work / 'q32-coordinate-unit.log').read_text()
    assert 'no more memory' in log and '_GB_DONE' not in log and '_COFACTOR_PASS' not in log
    source = (work / 'q32-coordinate-unit.sing').read_text()
    original = (previous / 'wronski-exact.sing').read_text().split('ideal I=', 1)[1].split(';\n', 1)[0].split(',')[:134]
    assert source.split('ideal I=', 1)[1].split(';\n', 1)[0].split(',') == original
    assert 'poly x32=cc;' in source
    inputs = [Path(__file__), previous / 'report.json', previous / 'restricted-original134.json.gz',
              ROOT / 'research/notes/S61_AFFINE_QUADRATIC_PENCIL_AUDIT_20261103_v62.md']
    answer = dict(id='S61_AFFINE_QUADRATIC_PENCIL_AUDIT_20261103_v62', all_checks_passed=True,
                  original_equations=134, quadratic_monomials=1035, exact_constant_point_rank=134,
                  exact_parameter_free_columns=25, exact_parameter_free_column_rank=16,
                  exact_each_parameter_coefficient_rank=coefficient_ranks,
                  selected_minor_logarithmic_derivatives_mod32003=traces,
                  selected_minor_nonconstant_over_QQ_in_each_parameter=True,
                  simultaneous_nilpotent_inverse_shortcut_for_selected_pivot=False,
                  image_flag_ranks_mod32003=image_ranks,
                  q32_algebraicity_probe=receipt, q32_algebraicity_proved=False,
                  affine_dimension_or_length_result_obtained=False,
                  limitations=['Only columns with identically parameter-free coefficients and one selected full minor are audited.',
                               'This does not rule out other constant-determinant minors, other preconditioners, or other algorithms.',
                               'The image-flag ranks are stated over F32003; no QQ equality for that flag is asserted.'],
                  input_sha256={str(path.relative_to(ROOT)): sha(path) for path in inputs},
                  evidence_sha256={path.name: sha(path) for path in args.output.iterdir() if path.is_file()},
                  seconds=time.monotonic() - started)
    (args.output / 'report.json').write_text(json.dumps(answer, indent=2) + '\n')
    print('PASS exact quadratic-pencil scope audit; parameter-free rank16; no q32 GB result;', answer['seconds'])


if __name__ == '__main__':
    main()
