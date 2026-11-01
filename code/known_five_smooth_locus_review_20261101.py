#!/usr/bin/env python3
"""Independent universal dual-Fisher recognition and MP basepoint evaluation."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61, X, coefficients

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    definition = json.loads((ROOT/'reports/TORIC_DEGREE5_COVARIANT_20261101.json').read_text())
    u = s.symbols(' '.join(definition['beta_variables']))
    tensor = s.Matrix(definition['weight_basis_24_by16'])*s.Matrix(definition['b_in_u'])*s.Matrix(u)
    j = s.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
    yy = s.Matrix(s.symbols('y0:4'))
    t = s.Symbol('t')
    coordinates = [*yy, t]
    pairs = list(combinations(range(4), 2))
    beta = s.zeros(4)
    for h, (a, b) in enumerate(pairs):
        beta[a, b] = sum(tensor[4*h+k]*yy[k] for k in range(4))
        beta[b, a] = -beta[a, b]
    phi = (beta+t*j).row_join(j*yy).col_join((-(j*yy).T).row_join(s.zeros(1)))
    pf = []
    for omitted in range(5):
        a, b, c, d = [k for k in range(5) if k != omitted]
        pf.append(s.expand((-1)**omitted*(phi[a, b]*phi[c, d]
                        -phi[a, c]*phi[b, d]+phi[a, d]*phi[b, c])))
    pf = s.Matrix(pf)
    assert (phi*pf).applyfunc(s.expand) == s.zeros(5, 1)
    derivatives = [pf.diff(z) for z in coordinates]
    matrices = [phi.diff(z) for z in coordinates]
    omega = s.zeros(5)
    for a, b in combinations(range(5), 2):
        c, d, e = [k for k in range(5) if k not in (a, b)]
        permutation = [a, b, c, d, e]
        sign = (-1)**sum(permutation[i] > permutation[k] for i in range(5) for k in range(i+1, 5))
        omega[a, b] = s.expand(sign*((derivatives[c].T*matrices[d]*derivatives[e])[0]
                         +(derivatives[e].T*matrices[c]*derivatives[d])[0]
                         +(derivatives[d].T*matrices[e]*derivatives[c])[0]))
        omega[b, a] = -omega[a, b]
    order = [4, 0, 1, 2, 3]
    P = omega.extract(order, order).subs(dict(zip([t, *yy], X)), simultaneous=True)
    assert P.subs(dict(zip(X, [1, 0, 0, 0, 0]))) == s.diag(s.zeros(1), -5*j)
    model = S61()
    normal = model.normalize(-P/5)
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    two = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 2]
    assert all(s.expand(normal['parameters'][i]-2*v) == 0 for i, v in zip(one, u))
    q51 = json.loads((ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json').read_text())
    substitution = dict(zip(u, [2*v for v in u]))
    assert all(s.expand(s.sympify(q).subs(substitution, simultaneous=True)-normal['parameters'][i]) == 0
               for i, q in zip(two, q51['quadratic_section_q']))
    assert all(s.expand(a-b) == 0 for a, b in zip(coefficients(model.reconstruct(normal['parameters'])), normal['coefficients']))
    # MP evaluation: the wedge of five projected primitive Z's is nonzero.
    Z = s.eye(6)[:, :5]
    Z[5, 0] = -1
    covector = s.Matrix([(-1)**i*Z[[k for k in range(6) if k != i], :].det() for i in range(6)])
    assert covector.T*Z == s.zeros(1, 5)
    wedge_pairing = s.zeros(6)
    for a, b, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        wedge_pairing[a, b] = wedge_pairing[b, a] = sign
    value = wedge_pairing.inv()*covector
    assert value == s.Matrix([-1, 0, 0, 0, 0, -1])
    paths = [Path(__file__), ROOT/'code/s61_model.py',
             ROOT/'reports/TORIC_DEGREE5_COVARIANT_20261101.json',
             ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',
             ROOT/'research/notes/KNOWN_FIVE_SMOOTH_LOCUS_20261101.md']
    report = dict(id='KNOWN_FIVE_SMOOTH_LOCUS_REVIEW_20261101', all_passed=True,
                  field='QQ[all16 primitive variables]',
                  independence='Constraint-basis beta and hand-coded complementary Pfaffians/cyclic Fisher formula; audited primitive_beta and pfaffOmega helpers not used.',
                  universal_normalized_Q51_2beta=True, full150_reconstruction=True,
                  MP_projected_primitive_basis=encode(Z), MP_basepoint_nonzero_bivector=list(map(str, value)),
                  primary_review=['Polishchuk 2204.09152v2 Theorem A and Introduction',
                                  'Markarian-Polishchuk 2310.18759v3 exterior construction, Theorem A, section3.2'],
                  accepted_conclusion='The reduced known three-component union is finite etale degree5 over the entire Delta!=0 open, with no additional open restriction.',
                  limitations=['No smoothness, reducedness, isolation, or exhaustion of the ambient full S61 fibers is inferred.',
                               'Q52 polynomial identification uses the established generic reconstruction and the polynomial MP construction reviewed in the note.'],
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'known_five_over_Delta_review': 'accepted'}))


def encode(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


if __name__ == '__main__':
    main()
