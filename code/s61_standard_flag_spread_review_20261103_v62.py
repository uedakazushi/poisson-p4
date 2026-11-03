#!/usr/bin/env python3
"""Audit full formal integrality, Borel action and D5 data for the standard flag."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, J, sp_basis, infinitesimal_action

ROOT = Path(__file__).resolve().parents[1]
P = 32003


def verify(output):
    source_path = ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    d5_path = ROOT/'reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json'
    report_path = ROOT/'reports/S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.json'
    source, d5, report = [json.loads(p.read_text()) for p in [source_path, d5_path, report_path]]
    b = s.symbols('b:15')
    substitution = {b[0]:1,b[1]:1,b[3]:1,b[12]:0}
    entries = set()
    for key in ['jacobian', 'tangent', 'left_kernel']:
        entries.update(v for row in source['linear'][key] for v in row)
    entries.update(v for row in source['quadratic']['quadratic_matrix'] for v in row)
    coefficient_count = 0
    for text in entries:
        numerator, denominator = s.fraction(s.cancel(s.sympify(text)))
        # Audit BEFORE specializing b0,b1,b3,b12, so terms vanishing there
        # cannot hide a bad p denominator in the full formal neighborhood.
        for polynomial in [numerator, denominator]:
            coefficients = s.Poly(polynomial, *b, domain=s.QQ).coeffs()
            assert all(int(c.q) % P for c in coefficients)
            coefficient_count += len(coefficients)
        value = s.cancel(denominator.subs(substitution))
        assert value.is_Rational and int(value.p) % P and int(value.q) % P
    z, u, beta = primitive_beta(S61())
    tensor = s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS])
    embedding = tensor.jacobian(u)
    rows = list(embedding.T.rref()[1])
    left = embedding[rows,:].inv()
    H = s.Matrix(source['linear']['beta_kernel']).applyfunc(s.sympify)
    hrows = list(H.T.rref()[1])
    hleft = H[hrows,:].inv()
    generator = s.Matrix(source['frame_change']['unipotent_generator']).applyfunc(s.Rational)
    g = s.eye(4)+generator
    assert generator**2 == s.zeros(4) and g.T*J*g == J
    shifted = g.T*beta.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g
    shifted_tensor = s.Matrix([s.expand(shifted[i,j]).diff(z[k]) for (i,j),k in COORDS])
    shifted_u = (left*shifted_tensor[rows,:]).subs(dict(zip(u,H*s.Matrix(b))), simultaneous=True).applyfunc(s.expand)
    shifted_b = (hleft*shifted_u[hrows,:]).applyfunc(s.expand)
    assert (H*shifted_b-shifted_u).applyfunc(s.expand) == s.zeros(16,1)
    assert shifted_b[3] == b[3]+b[0]
    assert all(shifted_b[i] == b[i] for i in [0,1,12])
    normal = []
    for a in sp_basis():
        action = left*infinitesimal_action(tensor,a)[rows,:]
        normal.append(s.expand(action[0].subs(dict(zip(u,H*s.Matrix(b))), simultaneous=True)))
    expected_normal = [s.expand(s.sympify(v).subs(dict(zip(u,H*s.Matrix(b))), simultaneous=True))
                       for v in d5['incidence_normal']]
    assert normal == expected_normal and normal[5] == -b[0]
    # The full D5 matrix is regular in all beta parameters, with only u5
    # inverted.  Verify arithmetic coefficients before imposing u5=1.
    for text in {v for row in d5['kernel'] for v in row}:
        numerator, denominator = s.fraction(s.cancel(s.sympify(text)))
        terms = s.Poly(denominator, s.Symbol('u5'), domain=s.QQ).terms()
        assert len(terms) == 1
        coefficient = terms[0][1]
        assert int(coefficient.p) % P and int(coefficient.q) % P
        assert all(int(c.q) % P for c in s.Poly(numerator,*u,domain=s.QQ).coeffs())
    witnesses = [s.sympify(source['linear']['jacobian_minor']).subs(substitution),
                 s.sympify(source['quadratic']['quadratic_minor']).subs(substitution),
                 s.sympify(d5['kernel_minor']).subs(s.Symbol('u5'),1),
                 s.sympify(d5['radial_C_normal_complement_minor']).subs(s.Symbol('u5'),1), -s.S.One]
    residues = []
    for value in witnesses:
        value = s.Rational(s.cancel(value))
        residue = int(value.p)*pow(int(value.q),-1,P) % P
        assert residue
        residues.append(residue)
    for path, digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    paths = [Path(__file__), source_path, d5_path, report_path,
             ROOT/'research/notes/S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.md']
    result = dict(id='S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62_REVIEW', all_passed=True,
                  universal_entries_checked_before_specialization=len(entries),
                  universal_numerator_and_denominator_coefficients_checked=coefficient_count,
                  no_bad_p_denominator_hidden_by_specialization=True,
                  finite_Borel_beta_action_independently_reconstructed=True,
                  flag_incidence_beta_smooth_component_independently_reconstructed=True,
                  D5_full_beta_integrality_checked=True, five_witness_residues=residues,
                  relative_closed_D5_and_transverse_upper5_argument_accepted=True,
                  all_h_standard_rho_axis_accepted=True, other_points_exhausted=False,
                  input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})
    with output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS', len(entries), 'universal entries;', coefficient_count, 'coefficients;', residues)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.output)
