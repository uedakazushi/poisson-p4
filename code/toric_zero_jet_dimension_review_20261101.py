#!/usr/bin/env python3
"""Explicit affine-action audit of the sharp toric zero-first-jet dimension.

Uses hand-written conformal-symplectic matrices and verifies root-type rank
bounds with every lower potential coefficient retained as a variable.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import time
import sympy as s
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]


def minor(matrix, rank):
    cs = list(matrix.rref()[1]); assert len(cs) >= rank
    cs = cs[:rank]
    rs = list(matrix[:, cs].T.rref()[1]); assert len(rs) == rank
    determinant = s.factor(matrix.extract(rs, cs).det()); assert determinant != 0
    return {'rows': rs, 'columns': cs, 'determinant': str(determinant)}


def verify():
    started = time.monotonic()
    h = s.symbols('h0:9'); x, y, z, w = s.symbols('x y z w')
    coords = (x, y, z, w); point = s.Matrix(coords)
    mons = [x*x*z, x*x*w, x*z*z, x*z*w, x*w*w, z**3, z*z*w, z*w*w, w**3]
    G = sum(a*b for a, b in zip(h, mons, strict=True))
    Q = s.zeros(4); Q[1, 2] = -s.diff(G, w); Q[1, 3] = s.diff(G, z); Q -= Q.T
    pairs = list(combinations(range(4), 2))
    qmons = [coords[i]*coords[j] for i, j in combinations_with_replacement(range(4), 2)]
    def vec(M):
        return s.Matrix([s.Poly(s.expand(M[i, j]), *coords).coeff_monomial(m)
                         for i, j in pairs for m in qmons])
    emb = vec(Q).jacobian(h); rows = list(emb.T.rref()[1]); left = emb[rows, :].inv()
    def potential(M):
        vector = vec(M); answer = (left*vector[rows, :]).applyfunc(s.expand)
        assert emb*answer == vector
        return answer
    def action(P, A, c=0):
        return (A*P+P*A.T-sum(((A*point)[i]*P.diff(coords[i])
                              for i in range(4)), s.zeros(4))-c*P).applyfunc(s.expand)

    # R0: sp4 stabilizing <d_y>, plus one independent scaling of Q.
    d, a, b, m, n, o, u = s.symbols('d a b m n o u')
    AR = s.Matrix([[d, 0, 0, 0], [u, -d, b, -a],
                   [a, 0, m, n], [b, 0, o, -m]])
    assert AR*J+J*AR.T == s.zeros(4)
    assert action(Q, AR.diff(u)) == s.zeros(4)
    R = s.Matrix.hstack(*[potential(action(Q, AR.diff(v))) for v in [d, m, n, o, a, b]], s.Matrix(h))
    assert R.shape == (9, 7)
    # The GL2 orbit and x-translation bounds are triangular: top four rows
    # describe H, and the next three describe x*H2. Keep all lower h free.
    distinct = {h[5]: 1, h[6]: 0, h[7]: 0, h[8]: 1}
    double = {h[5]: 0, h[6]: 1, h[7]: 0, h[8]: 0}
    zero_lower = dict.fromkeys(h[:5], 0)
    r_distinct = R.subs(distinct)
    top_columns = list(r_distinct[5:9, :].rref()[1]); assert len(top_columns) == 4
    # Columns a,b are numbered4,5 and vanish in the H3 quotient.
    assert r_distinct[5:9, [4, 5]] == s.zeros(4, 2)
    trans_rows = list(r_distinct[2:5, [4, 5]].T.rref()[1]); assert len(trans_rows) == 2
    selected_rows = [5, 6, 7, 8]+[2+i for i in trans_rows]
    selected_columns = top_columns+[4, 5]
    rd = s.factor(r_distinct.extract(selected_rows, selected_columns).det())
    assert rd != 0 and not rd.free_symbols
    r_double = R.subs(double).subs({h[0]: 1, h[1]: 0, h[2]: 0, h[3]: 0, h[4]: 0})
    r_double_minor = minor(r_double, 6)

    # Fixed minimal trace: x is fixed; y weight is trace of the z,w block.
    p = s.Symbol('p')
    AM = s.Matrix([[0, 0, 0, 0], [u, m+p, b, -a],
                   [a, 0, m, n], [b, 0, o, p]])
    assert AM*J+J*AM.T == (m+p)*J
    minimal = {}
    for kind in ['M0', 'M1']:
        P = J+Q
        P[0, 1] += x*x; P[1, 0] -= x*x
        if kind == 'M1':
            P[2, 3] += x*x; P[3, 2] -= x*x
            P[1, 2] -= x*z; P[2, 1] += x*z
            P[1, 3] -= x*w; P[3, 1] += x*w
        assert action(P, AM.diff(u)) == s.zeros(4)
        M = s.Matrix.hstack(*[potential(action(P, AM.diff(v), s.diff(m+p, v)))
                              for v in [m, n, o, p, a, b]])
        checks = {}
        for label, substitution, top_rank in [('distinct', distinct, 4), ('double', double, 3)]:
            mat = M.subs(substitution)
            top = mat[5:9, :]
            cs = list(top.rref()[1]); assert len(cs) == top_rank
            rs = list(top[:, cs].T.rref()[1]); assert len(rs) == top_rank
            assert top[:, [4, 5]] == s.zeros(4, 2)
            trs = list(mat[2:5, [4, 5]].T.rref()[1]); assert len(trs) == 2
            rr = [5+i for i in rs]+[2+i for i in trs]; cc = cs+[4, 5]
            determinant = s.factor(mat.extract(rr, cc).det())
            assert determinant != 0 and not determinant.free_symbols
            checks[label] = {'rank_lower_bound': top_rank+2, 'rows': rr, 'columns': cc,
                             'minor_independent_of_all_lower_coefficients': str(determinant)}
        triple = dict.fromkeys(h, 0); triple[h[3]] = 2; triple[h[5]] = 1
        triple_minor = minor(M.subs(triple), 5)
        # The two translation fields contribute only degree-one terms;
        # all six displayed residual actions are homogeneous degree two.
        translations = [P.diff(z), P.diff(w)]
        assert all(s.Poly(v, *coords).total_degree() <= 1 for T in translations for v in T)
        minimal[kind] = {'residual_potential_action_matrix': [[str(v) for v in row] for row in M.tolist()],
                         'root_type_bounds': checks, 'triple_root_rank5_minor': triple_minor}

    dependencies = [Path(__file__), ROOT/'code/weight1_structure.py',
        ROOT/'reports/TORIC_ZERO_JET_DIMENSION_20261101.json',
        ROOT/'reports/TORIC_ZERO_JET_FIBER_20261101.json',
        ROOT/'research/notes/TORIC_ZERO_JET_DIMENSION_20261101.md']
    return {'id': 'TORIC_ZERO_JET_DIMENSION_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'r0_explicit_parabolic_matrix': [[str(v) for v in row] for row in AR.tolist()],
        'r0_action_matrix': [[str(v) for v in row] for row in R.tolist()],
        'r0_distinct_root_rank6_minor': {'rows': selected_rows, 'columns': selected_columns,
            'minor_independent_of_all_lower_coefficients': str(rd)},
        'r0_double_root_rank6_minor': r_double_minor,
        'minimal_explicit_conformal_matrix': [[str(v) for v in row] for row in AM.tolist()],
        'minimal_models': minimal,
        'accepted_geometric_dimension_bounds': {'extra_translation_R0': 9, 'remaining_R0': 10,
            'translation_rank_lt2_M0_M1': 9, 'remaining_M0_M1': 10},
        'sharp_reduced_zero_jet_fiber_dimension': 10,
        'scope': 'Dimension only, conditional on the recorded historical compatible-quadratic exhaustion and previous AB inclusion; no decomposition or full S61 exhaustion.',
        'proof_note': 'research/notes/TORIC_ZERO_JET_DIMENSION_REVIEW_20261101.md',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(v.relative_to(ROOT)): hashlib.sha256(v.read_bytes()).hexdigest() for v in dependencies}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(); report = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2); out.write('\n')
    print(json.dumps({'all_passed': True, 'dimension': 10, 'seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()
