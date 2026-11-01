#!/usr/bin/env python3
"""Independent strict-action and full-matrix AB-degeneration review for R0."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
import sympy as s
from s61_model import X, quotient_coefficients
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.monotonic(); h = s.symbols('h0:9')
    prior_path = ROOT/'reports/TORIC_ZERO_JET_DIMENSION_REVIEW_20261101.json'
    prior = json.loads(prior_path.read_text())
    M = s.Matrix(prior['r0_action_matrix'])[:, :6]  # Remove the Q-scaling column.
    determinants = {}
    for name, tail, rows, cols in [
        ('distinct', [1, 0, 0, 1], [2, 4, 5, 6, 7, 8], list(range(6))),
        ('double', [0, 1, 0, 0], [2, 3, 5, 6, 7], [0, 2, 3, 4, 5])]:
        determinant = s.expand(M.subs(dict(zip(h[5:], tail, strict=True))).extract(rows, cols).det())
        assert determinant != 0 and not determinant.free_symbols
        determinants[name] = str(determinant)
    T = M.subs(dict(zip(h[5:], [1, 0, 0, 0], strict=True))); rows = [2, 5, 6]
    for cs in combinations(range(6), 3):
        determinant = s.expand(T.extract(rows, cs).det())
        if determinant != 0 and not determinant.free_symbols: cols = list(cs); break
    rest = [i for i in range(6) if i not in cols]
    N = (T[:, rest]-T[:, cols]*T.extract(rows, cols).inv()*T[rows, rest]).applyfunc(s.expand)
    minors = [s.expand(N.extract(rr, cc).det()) for rr in combinations(range(9), 2)
              for cc in combinations(range(3), 2)]
    for variable, extra, power in [(h[4], [], 2), (h[3], [h[4]], 3), (h[1], [h[4], h[3]], 2)]:
        gb = s.groebner(minors+extra, *h[:5])
        assert gb.reduce(variable**power)[1] == 0
    Hzero = M.subs(dict(zip(h[2:], [0, 1, 0, 0, 0, 0, 0], strict=True)))
    value = s.expand(Hzero.extract([0, 1, 2, 3, 4], [0, 2, 3, 4, 5]).det())
    assert value != 0 and not value.free_symbols

    x0, x, y, z, w = X
    P = s.zeros(5); P[1, 2] = P[3, 4] = x0*x0
    P[2, 3] = -z*z; P[2, 4] = 2*z*w; P -= P.T
    actions = []
    for i in range(5):
        for j in range(5):
            A = s.zeros(5); A[i, j] = 1
            actions.append(s.Matrix(quotient_coefficients(A*P+P*A.T-X[j]*P.diff(X[i]))))
    action = s.Matrix.hstack(*actions); assert action.rank() == 20
    E10, E20, E21 = [s.zeros(5) for _ in range(3)]
    E10[1, 0] = 1; E20[2, 0] = 1; E21[2, 1] = 1
    basis = [s.eye(5), E10, E20, E21, s.diag(0, -1, 1, 1, -1)/2]
    kernel = s.Matrix.hstack(*[A.reshape(25, 1) for A in basis])
    assert kernel.rank() == 5 and action*kernel == s.zeros(150, 5)
    assert E21*E10-E10*E21 == E20

    # Check the entire pulled-back4x4 bivector, rather than separate blocks.
    t, d = s.symbols('t d', nonzero=True); f0, f1, f2 = s.symbols('f0 f1 f2')
    A = s.Matrix([[t, 0], [2*t, d*t**3]]); S = s.zeros(4)
    for i, ii in enumerate([0, 2]):
        for j, jj in enumerate([0, 2]): S[ii, jj] = A[i, j]
    for i, ii in enumerate([1, 3]):
        for j, jj in enumerate([1, 3]): S[ii, jj] = A.inv().T[i, j]
    assert (S*J*S.T-J).applyfunc(s.cancel) == s.zeros(4)
    R, Z = A*s.Matrix([x, z]); F = f0*x*x+f1*x*z+f2*z*z
    source = s.zeros(4); source[0, 1] = 1-R*R; source[2, 3] = 1-Z*Z
    source[1, 3] = F/A.det(); source -= source.T
    pulled = (S.inv()*source*S.inv().T).applyfunc(s.cancel)
    assert all(s.denom(v).subs(t, 0) != 0 for v in pulled)
    limit = pulled.subs(t, 0); target = J.copy()
    target[1, 2] += 6*x*x/d; target[2, 1] -= 6*x*x/d
    target[1, 3] += F; target[3, 1] -= F
    assert (limit-target).applyfunc(s.expand) == s.zeros(4)
    a, b = s.symbols('a b'); change = s.Matrix([[1, 0, 0, 0], [0, 1, -b, a], [-a, 0, 1, 0], [-b, 0, 0, 1]])
    assert change*J*change.T == J
    assert change*s.Matrix([1, 0, a, b]) == s.Matrix([1, 0, 0, 0])
    assert change*s.Matrix([0, 1, 0, 0]) == s.Matrix([0, 1, 0, 0])
    paths = [Path(__file__), prior_path, ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py',
             ROOT/'reports/TORIC_ZERO_JET_R0_20261101.json']
    return {'id': 'TORIC_ZERO_JET_R0_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'independent_strict_constant_minors': determinants,
        'triple_root_radical_forcing': ['h4', 'h3', 'h1'],
        'binary_double_strict_dimension': 4, 'binary_double_kernel_nonabelian': True,
        'full_bivector_AB_deformation_limit_checked': True,
        'x_component_translation_symplectic_normalization_checked': True,
        'conclusion': 'Reduced F0 intersect R0 is contained in closed AB.',
        'scope': 'No scheme intersection or full zero-first-jet fiber exhaustion.',
        'proof_note': 'research/notes/TORIC_ZERO_JET_R0_REVIEW_20261101.md',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); report = verify(); args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out: json.dump(report, out, indent=2); out.write('\n')
    print(json.dumps({'all_passed': True, 'seconds': report['elapsed_seconds']}))


if __name__ == '__main__': main()
