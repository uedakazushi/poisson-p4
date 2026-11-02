#!/usr/bin/env python3
"""Exact product-pair reduction and generic cubic-product rank on the K=0 chart."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

import numpy as np
import sympy as s
from sympy.polys.rings import ring
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank, determinant_mod

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'research/notes/S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102.md'
INPUT = ROOT / 'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'


def verify():
    start = time.time()
    data = json.loads(INPUT.read_text())
    assert data['all_passed']
    for path, expected in data['input_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected
    R, a, e, d = ring('a,e,d', s.QQ)
    P, *v = ring('v0,v1,v2,v3,v4,v5,v6,v7', R)
    quadrics = [P.from_dict({tuple(mon): R.from_dict({tuple(mm): s.QQ(c) for mm, c in terms})
                           for mon, terms in row}) for row in data['quadratic_rows_sparse']]
    assert len(quadrics) == 6
    Z, W = 3*v[0] - v[1] + v[6], v[3]
    A = -(d*e+8)*quadrics[1] + 8*quadrics[2] + d*quadrics[3]
    B = (-10*a-3*e**2+3*e)*quadrics[1] + 3*(e-1)*quadrics[3] + 10*a*quadrics[4]
    linear = A.exquo(W)
    assert max(map(sum, linear)) == 1
    assert B == -a*Z*linear
    assert a*Z*A + W*B == 0
    # Retaining q0,q1,q3,q5 and replacing q2,q4 by A,B has determinant 80*a.
    change = s.eye(6)
    change[2, :] = s.Matrix([[0, -(d*e+8).as_expr(), 8, d.as_expr(), 0, 0]])
    change[4, :] = s.Matrix([[0, (-10*a-3*e**2+3*e).as_expr(), 0, 3*(e-1).as_expr(), 10*a.as_expr(), 0]])
    assert s.expand(change.det()) == 80*a.as_expr()
    light = [0, 1, 2, 3, 4, 6, 7]
    products = [v[j]*q for q in quadrics for j in light]
    monomials = sorted(set().union(*(q.keys() for q in products)))
    assert len(products) == 42 and len(monomials) == 84
    # This polynomial relation has a unit coefficient 10*a on W*q4 on the a-open.
    relation = [P.zero for _ in quadrics]
    relation[1] = -a*(d*e+8)*Z + (-10*a-3*e**2+3*e)*W
    relation[2] = 8*a*Z
    relation[3] = a*d*Z + 3*(e-1)*W
    relation[4] = 10*a*W
    assert sum((f*q for f, q in zip(relation, quadrics)), P.zero) == 0
    prime = 32003
    point = (1, 2, 3)
    def evaluate(poly):
        value = s.Rational(poly.evaluate(0, point[0]).evaluate(0, point[1]).evaluate(0, point[2]))
        assert int(value.q) % prime
        return int(value.p) % prime * pow(int(value.q), -1, prime) % prime
    matrix = np.array([[evaluate(q.get(mon, R.zero)) for mon in monomials] for q in products], dtype=np.int64)
    rank, selected_rows = elimination_rank(matrix.T, prime)
    _, selected_columns = elimination_rank(matrix[selected_rows], prime)
    minor = determinant_mod(matrix[selected_rows][:, selected_columns], prime)
    assert rank == 41 and minor
    factors = [a, 500*a+90*e**2+3*e-24,
               50*a*d+18*d*e**2-9*d*e-2*d+60*e-40,
               100*a+18*e**2-9*e-2]
    assert all(evaluate(f) for f in factors)
    paths = [Path(__file__), NOTE, INPUT,
             ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py']
    return {
        'id': 'S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102',
        'all_passed': True,
        'quadratic_source_rows': data['quadratic_rows'],
        'quadric_indexing': 'q0,...,q5 in the stored source-row order',
        'Z': str(Z.as_expr()), 'W': str(W.as_expr()),
        'linear_form': str(linear.as_expr()),
        'change_determinant': '80*a',
        'product_pair_identities': ['A=W*ell', 'B=-a*Z*ell'],
        'syzygy_coefficients': [str(f.as_expr()) for f in relation],
        'product_matrix_shape': [42, 84],
        'generic_product_rank': 41,
        'rank_witness': {'prime': prime, 'point': point, 'rows': selected_rows,
                         'columns': selected_columns, 'minor': minor},
        'scope': ['The six-quadratic ideal is (q0,q1,q3,q5,W*ell,Z*ell) on a*F*B*C nonzero.',
                  'The cubic-product rank is at most 41 at every such point and exactly 41 generically.',
                  'No assertion that ell,Z,W are independent at every point.',
                  'No complete cubic or quartic rank coverage, length65, or ambient saturation removal.'],
        'proof_note': str(NOTE.relative_to(ROOT)),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'seconds': time.time()-start,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ['all_passed', 'generic_product_rank', 'seconds']}))
