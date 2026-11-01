#!/usr/bin/env python3
"""Independent linear-algebra review of the null-correlation dual curve.

Builds the primitive tensor from the closed/primitive constraint kernel,
without the audited S61 primitive_beta or Pfaffian helper implementations.
"""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import sympy as s

from weight1_structure import constraints

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    j = s.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0],
                  [0, 0, 0, 1], [0, 0, -1, 0]])
    y = s.Matrix(s.symbols('y0:4'))
    b = s.symbols('b0:16')
    basis = s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    vector = basis*s.Matrix(b)
    pairs = list(combinations(range(4), 2))
    beta = s.zeros(4)
    for h, (i, k) in enumerate(pairs):
        beta[i, k] = sum(vector[4*h+a]*y[a] for a in range(4))
        beta[k, i] = -beta[i, k]
    assert beta[0, 1]+beta[2, 3] == 0
    alpha = -beta*y
    derivative = s.Matrix(4, 4, lambda i, k: s.diff(alpha[k], y[i])-s.diff(alpha[i], y[k]))
    assert (derivative-3*beta).applyfunc(s.expand) == s.zeros(4)
    t = s.Symbol('t')
    phi = s.zeros(5)
    phi[:4, :4] = beta+t*j
    phi[:4, 4] = j*y
    phi[4, :4] = -(j*y).T
    pfaffians = []
    for omitted in range(5):
        i, k, l, m = [a for a in range(5) if a != omitted]
        pfaffians.append(s.expand((-1)**omitted*(phi[i, k]*phi[l, m]
                           -phi[i, l]*phi[k, m]+phi[i, m]*phi[k, l])))
    pf = s.expand(beta[0, 1]*beta[2, 3]-beta[0, 2]*beta[1, 3]+beta[0, 3]*beta[1, 2])
    expected = (-j*(beta-t*j)*y).col_join(s.Matrix([t*t+pf]))
    assert (s.Matrix(pfaffians)-expected).applyfunc(s.expand) == s.zeros(5, 1)
    assert ((j*beta)**2+pf*s.eye(4)).applyfunc(s.expand) == s.zeros(4)
    assert ((j*beta-t*s.eye(4))*j*(beta-t*j)*y+(pf+t*t)*y).applyfunc(s.expand) == s.zeros(4, 1)
    assert pfaffians[4].subs(dict.fromkeys(y, 0)) == t*t

    zz = s.symbols('z0:5')
    z = s.zeros(4)
    for (i, k), value in zip(pairs, [*zz, -zz[0]], strict=True):
        z[i, k], z[k, i] = value, -value
    ell = s.Matrix([sum(z[i, k]*s.diff(beta[i, k], y[a]) for i, k in pairs)
                    for a in range(4)])
    graph = z.row_join(-j*ell).col_join((j*ell).T.row_join(s.zeros(1)))
    annihilation = s.expand(sum(graph[i, k]*phi[i, k] for i, k in combinations(range(5), 2)))
    assert annihilation == 0
    canonical = (j*z*j).row_join(ell).col_join((-ell.T).row_join(s.zeros(1)))
    congruence = s.diag(j, s.ones(1))
    assert (congruence*canonical*congruence.T+graph).applyfunc(s.expand) == s.zeros(5)
    # Uniform injectivity uses the graph projection and the last column Jy.
    assert s.Matrix([graph[i, k] for i, k in pairs]).jacobian(zz).rank() == 5
    assert s.Matrix([*phi[:4, 4], phi[0, 1]]).jacobian([*y, t]).det() == 1
    paths = [Path(__file__), ROOT/'code/weight1_structure.py',
             ROOT/'code/null_correlation_duality_20261101.py',
             ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md']
    report = dict(id='NULL_CORRELATION_DUALITY_REVIEW_20261101', all_passed=True,
                  field='QQ[b0,...,b15,y0,...,y3,t,z0,...,z4]',
                  independent_construction='Closed/primitive constraint kernel; explicit complementary-index Pfaffians; no audited primitive_beta or Pfaffian helpers',
                  Cartan_identity=True, all_five_Pfaffians=True,
                  chartwise_scheme_inverse_identity=True, projection_center_absent=True,
                  universal_annihilator_identity=True, canonical_congruence=True,
                  both_parameter_spaces_injective_dimension=5,
                  primary_read='Shinder-Zhang, DOI 10.1112/blms.12339, Proposition 2.8 and its hypotheses/proof',
                  accepted_conclusion='For Delta != 0, the null-correlation zero scheme is the orthogonal elliptic quintic Pic^3(C_beta), with the same Jacobian; not a claimed canonical equality of generic torsors.',
                  limitations=['The bundle cohomology and chartwise scheme argument are reviewed in the note, not inferred from point sampling.',
                               'No new Fisher invariant contraction or full-S61/affine-infinity conclusion.',
                               'The old G6 equivariance proof is reviewed as a dependency; its 36-dimensional syzygy computation is not repeated.'],
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'duality_review': 'accepted'}))


if __name__ == '__main__':
    main()
