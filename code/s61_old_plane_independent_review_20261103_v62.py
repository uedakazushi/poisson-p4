#!/usr/bin/env python3
"""Independent bounded audit of the old marked-square plane incidence data.

Use the original135 equations to check the reported entire H, reconstruct the
primitive jet action from tensors, and prove rank bounds using covectors plus
one exact lower witness.  Do not repeat the author's generic fraction-free
rank computation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM

from s61_marked_sigma_point_lift_20261103_v62 import load_model, coerce
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61, X
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import COORDS, sp_basis, infinitesimal_action

ROOT = Path(__file__).resolve().parents[1]


def dm(rows, domain):
    return DomainMatrix.from_rep(DDM(rows, (len(rows), len(rows[0])), domain))


def verify(output, model):
    started = time.monotonic()
    producer = ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_GEOMETRY_20261103_v62.json'
    evidence = json.loads(producer.read_text())
    K = s.QQ.algebraic_field(s.sqrt(-2), alias='h')
    hsymbol = s.Symbol('h')
    parse = lambda value: K.from_sympy(s.sympify(value).subs(hsymbol, s.sqrt(-2)))
    lifter = load_model(model)
    one, two = lifter.one, lifter.two
    geom = coords()

    # Rebuild the stored primitive representation from actual tensor pullbacks.
    z, beta_variables, beta_form = primitive_beta(S61())
    beta_tensor = s.Matrix([beta_form[i, j].diff(z[k]) for (i, j), k in COORDS])
    B = beta_tensor.jacobian(beta_variables)
    beta_rows = list(B.T.rref()[1])
    left = B[beta_rows, :].inv()
    primitive_actions = []
    for a in sp_basis():
        action = (left*infinitesimal_action(beta_tensor, a)[beta_rows, :]).applyfunc(s.expand)
        assert B*action == infinitesimal_action(beta_tensor, a)
        primitive_actions.append(action.jacobian(beta_variables))
    stored = s.Matrix(json.loads((ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json').read_text())['primitive_action']).applyfunc(s.sympify)
    assert s.Matrix.hstack(*(a*s.Matrix(beta_variables) for a in primitive_actions)) == stored
    results = []
    for case in evidence['cases']:
        H = dm([[parse(v) for v in row] for row in case['H']], K)
        assert H.rank() == case['H_dimension']
        q = list(map(parse, case['q']))
        zz = s.symbols(f'z0:{H.shape[1]}')
        R = K.poly_ring(*zz)
        Hp = H.convert_to(R)
        beta = Hp*dm([[R.gens[i]] for i in range(H.shape[1])], R)
        beta_values = [row[0] for row in beta.to_list()]
        parameters = dict(zip(one, beta_values)) | dict(zip(two, [R.convert(v, K) for v in q]))
        # At t=0 exactly the q-degree-two terms survive.  This checks every
        # equation directly over the entire reported H, including G6.
        for equation in lifter.original:
            total = R.zero
            for mon, n, d in equation:
                if sum(i in two for i in mon) != 2:
                    continue
                term = coerce(R, s.Rational(n, d))
                for i in mon:
                    term *= parameters[i]
                total += term
            assert not total

        # Independently reconstruct the full highest linear beta matrix.
        M = [[K.zero]*16 for _ in range(135)]
        for row, equation in enumerate(lifter.original):
            for mon, n, d in equation:
                beta_indices = [one.index(i) for i in mon if i in one]
                if sum(i in two for i in mon) != 2 or len(beta_indices) != 1:
                    continue
                value = coerce(K, s.Rational(n, d))
                for i in mon:
                    if i in two:
                        value *= q[two.index(i)]
                M[row][beta_indices[0]] += value
        assert dm(M, K).rank() + H.shape[1] == 16
        assert (dm(M, K)*H).is_zero_matrix

        # Projective q orbit rank from physical sixty coefficients, without
        # using the author's Pi-projection or its q25 chart correction.
        physical = geom['E']*s.Matrix([K.to_sympy(v) for v in q])
        y = s.Matrix(X[1:])
        tensor = s.zeros(4)
        for k, (i, j) in enumerate(pairs):
            tensor[i, j] = sum(physical[k*10+t]*y[a]*y[b] for t, (a, b) in enumerate(m2))
        tensor -= tensor.T
        orbit_columns = []
        for a in sp_basis():
            derivative = a*tensor+tensor*a.T-sum(((a*y)[i]*tensor.diff(y[i]) for i in range(4)), s.zeros(4))
            orbit_columns.append([K.from_sympy(s.Poly(derivative[i, j], *y).coeff_monomial(y[a]*y[b]))
                                  for i, j in pairs for a, b in m2])
        physical_values = list(map(K.from_sympy, physical))
        projective_rank = dm(list(map(list, zip(physical_values, *orbit_columns))), K).rank()-1
        assert projective_rank == case['orbit_dimension']

        # Full polynomial covectors give a generic upper bound; a single exact
        # point proves the matching lower bound for H + all ten beta actions.
        action_columns = []
        for action in primitive_actions:
            action_columns.append([sum((coerce(R, c)*b for c, b in zip(row, beta_values)), R.zero)
                                   for row in action.tolist()])
        Ap = dm(list(map(list, zip(*action_columns))), R)
        def parse_polynomial(value):
            out = R.zero
            for mon, c in s.Poly(s.sympify(value), *zz, hsymbol, domain=s.QQ).terms():
                term = R.convert(coerce(K, c)*K.unit**mon[-1], K)
                for i, e in enumerate(mon[:-1]):
                    term *= R.gens[i]**e
                out += term
            return out
        covector_rows = [[parse_polynomial(v) for v in row]
                         for row in case['annihilating_covectors']]
        if covector_rows:
            covectors = dm(covector_rows, R)
            assert (covectors*Hp).is_zero_matrix and (covectors*Ap).is_zero_matrix
            upper = 16-len(covector_rows)
        else:
            upper = 16
        target = case['generic_beta_differential_rank']
        assert upper == target
        witness = None
        for seed in range(1, 12):
            values = [K.convert(seed+i) for i in range(H.shape[1])]
            def evaluate(poly):
                return sum((c*s.prod(values[i]**e for i, e in enumerate(mon) if e)
                            for mon, c in poly.items()), K.zero)
            action_numeric = dm([[evaluate(v) for v in row] for row in Ap.to_list()], K)
            rank = H.hstack(action_numeric).rank()
            if covector_rows:
                assert dm([[evaluate(v) for v in row] for row in covectors.to_list()], K).rank() == len(covector_rows)
            if rank == target:
                witness = seed
                break
        assert witness is not None
        results.append(dict(case=case['case'], entire_H_satisfies_original135=True,
                            highest_kernel_dimension=H.shape[1], projective_q_orbit_dimension=projective_rank,
                            generic_beta_rank=target, independent_upper_covectors=len(covector_rows),
                            exact_lower_witness_seed=witness))
        print('PASS', results[-1], flush=True)
    # The binary discriminant and the finite SL2 double-root frame are exact.
    Z, W, u, v = s.symbols('Z W u v')
    G = -Z**3/3 + s.sqrt(-2)*u*Z*W**2/4 + v*W**3/4
    discriminant = s.discriminant(s.expand(12*G).subs(W, 1), Z)
    assert s.simplify(discriminant + 432*(2*s.sqrt(-2)*u**3+9*v**2)) == 0
    assert s.simplify(G.subs({u:s.sqrt(-2), v:2*s.sqrt(-2)/3})
                      + (Z-s.sqrt(-2)*W/2)**2*(Z+s.sqrt(-2)*W)/3) == 0
    frame = s.Matrix([[1,-s.sqrt(-2)/2], [2/(3*s.sqrt(-2)),s.Rational(2,3)]])
    assert s.simplify(frame.det()) == 1
    paths = [Path(__file__), producer, ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
             ROOT/'code/primitive_jet_elliptic_curve_20261101.py', ROOT/'code/weight1_structure.py',
             ROOT/'code/toric_closed_jet_coordinates_20261102.py']
    report = dict(id='S61_MARKED_SIGMA_OLD_PLANE_GEOMETRY_20261103_v62_REVIEW', all_passed=True,
                  primitive_beta_action_reconstructed=True, cases=results,
                  binary_discriminant_and_SL2_frame_checked=True,
                  frozen_model_input_sha256=hashlib.sha256(lifter.input_path.read_bytes()).hexdigest(),
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  ambient_nilpotence_asserted=False, seconds=time.monotonic()-started)
    with output.open('x') as stream:
        json.dump(report, stream, indent=2);stream.write('\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--model-folder', type=Path, required=True)
    args = parser.parse_args()
    print('COMPLETE', verify(args.output, args.model_folder)['seconds'])
