#!/usr/bin/env python3
"""Direct dual-number substitution in all135 equations for the old incidence."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM

from s61_marked_sigma_point_lift_20261103_v62 import load_model, coerce
from s61_model import X
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import sp_basis

ROOT = Path(__file__).resolve().parents[1]


def dm(rows, domain):
    return DomainMatrix.from_rep(DDM(rows, (len(rows), len(rows[0])), domain))


def verify(model_folder, output):
    source = ROOT/'reports/S61_MARKED_SIGMA_OLD_DOUBLE_QUOTIENT_20261103_v62_FINAL/report.json'
    report = json.loads(source.read_text())
    L = load_model(model_folder)
    K = s.QQ.algebraic_field(s.sqrt(-2), alias='h')
    h = s.Symbol('h')
    z = s.symbols('z0:9')
    R = K.poly_ring(*z)
    def parse(value):
        result = R.zero
        for mon, c in s.Poly(s.sympify(value), *z, h, domain=s.QQ).terms():
            term = R.convert(coerce(K, c)*K.unit**mon[-1], K)
            for i, e in enumerate(mon[:-1]):
                term *= R.gens[i]**e
            result += term
        return result
    H = [[parse(v) for v in row] for row in report['H']]
    beta = [sum((v*R.gens[i] for i, v in enumerate(row)), R.zero) for row in H]
    delta_beta = list(map(parse, report['delta_beta']))
    eta_q = list(map(parse, report['eta_q']))
    normal = [K.zero]*13
    normal[9], normal[10] = K.unit, K.one
    point = L.lift(K.unit, normal, domain=K)
    assert point['all135_zero'] and not point['t']
    q = point['q']
    assert not eta_q[25]
    parameters = dict(zip(L.one, zip(beta, delta_beta)))
    parameters.update(dict(zip(L.two, [(R.convert(v, K), d) for v, d in zip(q, eta_q)])))
    def multiply(a, b):
        return a[0]*b[0], a[0]*b[1]+a[1]*b[0]
    for equation in L.original:
        constant, epsilon = R.zero, R.zero
        for mon, n, d in equation:
            term = coerce(R, s.Rational(n, d)), R.zero
            for i in mon:
                term = multiply(term, parameters[i])
            for _ in range(2-sum(i in L.two for i in mon)):
                term = multiply(term, (R.zero, R.zero))  # t is identically zero
            constant += term[0]
            epsilon += term[1]
        assert not constant and not epsilon
    geom = coords()
    y = s.Matrix(X[1:])
    physical = geom['E']*s.Matrix([K.to_sympy(v) for v in q])
    Q = s.zeros(4)
    for k, (i, j) in enumerate(pairs):
        Q[i, j] = sum(physical[10*k+t]*y[a]*y[b] for t, (a, b) in enumerate(m2))
    Q -= Q.T
    orbit = []
    for a in sp_basis():
        derivative = a*Q+Q*a.T-sum(((a*y)[i]*Q.diff(y[i]) for i in range(4)), s.zeros(4))
        orbit.append([K.from_sympy(s.Poly(derivative[i, j], *y).coeff_monomial(y[a]*y[b]))
                      for i, j in pairs for a, b in m2])
    # Include q itself to work in physical projective coefficient space.
    base = dm(list(map(list, zip(list(map(K.from_sympy, physical)), *orbit))), K)
    assert base.rank() == 9  # eight projective orbit directions + radial direction
    def evaluate(poly, values):
        return sum((c*s.prod(values[i]**e for i, e in enumerate(mon) if e)
                    for mon, c in poly.items()), K.zero)
    H_constant = dm([[evaluate(v, [K.zero]*9) for v in row] for row in H], K)
    hrows = list(H_constant.transpose().rref()[1])
    actual_z = H_constant.extract(hrows, range(9)).inv()*dm([[point['beta'][i]] for i in hrows], K)
    samples = [dict(name='generic-integer-H', values=[K.convert(i+1) for i in range(9)]),
               dict(name='actual-central-beta', values=[row[0] for row in actual_z.to_list()])]
    witnesses = []
    for sample in samples:
        vector = [evaluate(v, sample['values']) for v in eta_q]
        extra = geom['E']*s.Matrix([K.to_sympy(v) for v in vector])
        rank = base.hstack(dm([[K.from_sympy(v)] for v in extra], K)).rank()
        assert rank == 10
        augmented = base.hstack(dm([[K.from_sympy(v)] for v in extra], K))
        residue_ranks = []
        for root in [3083, (-3083) % 32003]:
            ground = s.GF(32003)
            def residue(value):
                result = 0
                for c in value.to_list():
                    assert int(c.denominator) % 32003
                    result = (result*root+int(c.numerator)*pow(int(c.denominator), -1, 32003)) % 32003
                return ground.convert(result)
            rb = dm([[residue(v) for v in row] for row in base.to_list()], ground).rank()
            ra = dm([[residue(v) for v in row] for row in augmented.to_list()], ground).rank()
            assert (rb, ra) == (9, 10)
            residue_ranks.append(dict(h_root=root, radial_and_orbit_rank=rb, augmented_rank=ra))
        actual_beta = [evaluate(v, sample['values']) for v in beta]
        if sample['name'] == 'actual-central-beta':
            assert actual_beta == point['beta']
        witnesses.append(dict(name=sample['name'], augmented_physical_rank=rank,
                              normal_to_projective_q_orbit=True, residue_ranks=residue_ranks))
    for p, digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest
    paths = [Path(__file__), source, ROOT/'code/s61_marked_sigma_old_double_quotient_20261103_v62.py',
             ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_GEOMETRY_20261103_v62_REVIEW.json',
             ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py']
    result = dict(id='S61_MARKED_SIGMA_OLD_DOUBLE_QUOTIENT_20261103_v62_REVIEW', all_passed=True,
                  direct_entire_H_dual_number_substitution_all135=True,
                  q25_preserved=True, t_identically_zero=True,
                  independent_physical_projective_orbit_rank=8, normal_witnesses=witnesses,
                  same_incidence_generic_length_lower_bound=2,
                  ambient_equality_requires_separate_transverse_upper=True,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')
    print('PASS ALL135 DUAL NUMBERS; NORMAL WITNESSES', witnesses)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-folder', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.model_folder, args.output)
