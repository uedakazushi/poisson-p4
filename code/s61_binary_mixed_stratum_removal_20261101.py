#!/usr/bin/env python3
"""Relative C3 binary/mixed junction: exact smooth-source checks and proof dependencies."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s
from s61_model import S61, X, J as J0
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_minimal_infinity_local_data_20261101 import calculate

ROOT = Path(__file__).resolve().parents[1]


def verify():
    start = time.time()
    stems = [
        'S61_BINARY_MIXED_CONNECTION', 'S61_BINARY_MIXED_FORMAL',
        'S61_BINARY_MIXED_FORMAL_REVIEW', 'S61_BINARY_MIXED_RELATIVE_DOUBLE_LINE',
        'S61_BINARY_UNIFORM_RANK', 'S61_BINARY_UNIFORM_QUADRATIC',
        'S61_BINARY_UNIFORM_CUBIC',
    ]
    paths = [ROOT / 'reports' / (stem + '_20261101.json') for stem in stems]
    reports = {stem: json.loads(path.read_text()) for stem, path in zip(stems, paths)}
    assert all(r['all_passed'] for r in reports.values())
    old = reports['S61_BINARY_MIXED_CONNECTION']
    formal = reports['S61_BINARY_MIXED_FORMAL']
    double = reports['S61_BINARY_MIXED_RELATIVE_DOUBLE_LINE']
    assert formal['associated_graded_ideal_identified']
    assert formal['completed_homogenizer_square_zero']
    assert formal['completed_homogenizer_nonzero']
    assert formal['initial_ideal_degree'] == 65
    assert old['joint_image_dimension'] == 16
    assert old['point']['characteristic_zero_transverse_length_upper_bound'] == 67
    assert double['flat_relative_double_line']
    assert double['all_135_actual_double_line_identities']
    assert reports['S61_BINARY_UNIFORM_CUBIC']['normal_length_lower_bound'] == 65

    # Reconstruct the actual binary family, keeping beta symbolic before the point check.
    model = S61()
    z, uu, beta = primitive_beta(model)
    sample = {s.Symbol(k): s.Rational(v) for k, v in old['point']['sample'].items()}
    a, c, d, h, k = s.symbols('a c d h k')
    variables = [a, c, d, h, k]
    v1 = s.Matrix([a, 1, 0, c])
    v2 = s.Matrix([-c, 0, 1, d])
    origin = dict.fromkeys(variables, s.S.Zero)
    section = s.Matrix([
        (v1.T * beta.subs(dict(zip(z, v1))) * v2)[0],
        (v1.T * beta.subs(dict(zip(z, v2))) * v2)[0],
    ]).applyfunc(s.expand)
    assert section.subs(origin).subs(sample) == s.zeros(2, 1)
    vertex_jac = section.jacobian([a, c, d]).subs(origin).subs(sample)
    vertex_pivots = list(vertex_jac.to_DM().rref()[1])
    assert len(vertex_pivots) == 2
    vertex_minor = vertex_jac[:, vertex_pivots].det()
    direction = vertex_jac.nullspace()[0]

    two = [i for i, wt in enumerate(model.reduced['parameter_weights']) if wt == 2]
    fixed = two.index(29)
    chart = [i for i in range(45) if i != fixed]
    q0 = s.zeros(45, 1)
    q0[fixed] = 1
    xx, yy, zz, ww = X[1:]
    ell1, ell2 = xx-a*yy+c*zz, ww-c*yy-d*zz
    binary = (ell1*ell2+h*ell1**2+k*ell2**2) * (v1*v2.T-v2*v1.T)
    derivatives = []
    for var in variables:
        P = X[0]**2 * J0
        P[1:, 1:] += binary.diff(var).subs(origin)
        q = s.Matrix([model.normalize(P)['parameters'][i] for i in two])
        q -= q0*q[fixed]
        derivatives.append(q[chart, :].col_join(s.zeros(1, 1)))
    tangent = s.Matrix.hstack(
        s.Matrix.hstack(*derivatives[:3])*direction, derivatives[3], derivatives[4])
    free = formal['point']['tangent_free_rows'] if 'tangent_free_rows' in formal['point'] else [19,20,22,23,24,27,32,33,34,36,37,43]
    T = s.Matrix(formal['point']['tangent_basis']).applyfunc(s.Rational)
    assert T[free, :] == s.eye(12)
    coordinates = tangent[free, :]
    assert T*coordinates == tangent
    assert coordinates == s.Matrix(old['point']['actual_binary_family_tangent']).applyfunc(s.Rational)
    base = [0, 1, 3]
    base_minor = coordinates[base, :].det()
    assert base_minor != 0

    # Q0 and Q1 are in the same projective quadratic-form fiber of the vertex family.
    P = X[0]**2 * J0
    P[1:, 1:] += binary.subs({a: 0, c: 0, d: 0, h: 1, k: 1})
    q1 = s.Matrix([model.normalize(P)['parameters'][i] for i in two])
    assert q1 == s.Matrix(old['point']['same_beta_regular_point_q']).applyfunc(s.Rational)
    p61 = s.zeros(61, 1)
    for symbol, value in sample.items():
        p61[int(str(symbol)[1:])-1] = value
    for i, val in zip(two, q1):
        p61[i] = val
    local = calculate(model, True, p61)
    assert local['leading_jacobian_rank'] == 34
    assert local['homogenized_jacobian_rank'] == 35
    assert local['quadratic_span_rank'] == 5
    assert local['common_hessian_rank'] == 7
    assert local['homogenizer_zero_on_tangent']

    # Check the already certified smooth base-change minor and exact source point.
    normal = s.Matrix(double['beta_normal_matrix']).applyfunc(s.Rational)
    columns = double['beta_minor_columns']
    normal_minor = normal[:, columns].det()
    assert normal_minor == -s.Rational(1133614671336, 625)
    assert normal_minor == s.Rational(double['beta_minor'])
    assert s.Rational(old['beta_dominance_minor']) != 0
    enc = lambda M: [[str(v) for v in row] for row in M.tolist()]
    paths += [
        Path(__file__), ROOT/'code/s61_model.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'code/s61_minimal_infinity_local_data_20261101.py',
        ROOT/'research/notes/S61_BINARY_MIXED_FORMAL_20261101.md',
        ROOT/'research/notes/S61_BINARY_MIXED_FORMAL_REVIEW_20261101.md',
        ROOT/'research/notes/S61_BINARY_MIXED_RELATIVE_DOUBLE_LINE_20261101.md',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        ROOT/'research/notes/S61_BINARY_MIXED_STRATUM_REMOVAL_20261101.md',
    ]
    return dict(
        id='S61_BINARY_MIXED_STRATUM_REMOVAL_20261101', all_passed=True, field='QQ',
        vertex_section=list(map(str, section)), vertex_jacobian=enc(vertex_jac),
        vertex_minor_columns=vertex_pivots, vertex_minor=str(vertex_minor),
        actual_binary_family_tangent=enc(coordinates), binary_base_columns=base,
        binary_base_minor=str(base_minor),
        same_vertex_projective_quadratic_fiber=True,
        same_component_regular_point=dict(leading_rank=34, full_rank=35,
            quadratic_rank=5, common_hessian_radical_dimension=3,
            full_raw_local_data=local),
        smooth_beta_base_change_minor=str(normal_minor),
        relative_beta_dimension=9, relative_binary_dimension=12,
        actual_flat_double_line_dimension=10,
        point_kernel_torsion_free_rank=65,
        relative_kernel_torsion_free_rank=65,
        algebraic_proof_dependencies=[
            'The specified point theorem proves the actual kernel N0 finite torsion-free of rank65.',
            'The actual closed relative double line is beta-flat, so Nrel tensor k equals N0.',
            'The relative binary family is smooth of dimension3 over the nine beta parameters.',
            'Its general point has normal length at least65 by the universal quadratic/cubic certificates.',
            'Equal-rank torsion-free specialization lemma gives Nrel torsion-free over the twelve-variable regular base.',
            'Generic binary length65 gives t^2=0 generically; the double line has t=0; torsion-freeness gives t^2=0 everywhere.',
            'The smooth Sp4 base change is faithfully flat; homogenized equivariance descends t^2=0 to the ambient local ring.',
            'A nonempty annihilator open on the irreducible16-dimensional stratum removes its entire sufficiently general beta fiber.',
        ],
        theorem_conclusions=dict(relative_homogenizer_square_zero=True,
            ambient_local_homogenizer_square_zero=True,
            entire_general_beta_C3_stratum_removed=True),
        limitations=[
            'The ambient t^2 identity is asserted on a nonempty open, not on every special point of the stratum.',
            'The entire sufficiently general beta fiber of this sixteen-dimensional C3 stratum is removed.',
            'Other binary strata, flag points, affine finiteness and full S61 exhaustion remain separate.',
            'No freeness of Nrel or Srel and no rank67 claim is made.',
        ],
        proof_note='research/notes/S61_BINARY_MIXED_STRATUM_REMOVAL_20261101.md',
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({k: report[k] for k in ['all_passed', 'seconds']}))
