#!/usr/bin/env python3
"""A polynomial E2-root coordinate on the reduced toric S61 component.

The tensor contraction is defined on every S61 point. Its identification with
an E2 root is proved on the toric component, not on the full S61 scheme.
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
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61
from q51_universal_section_20261101 import affine, exterior_derivative
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from toric_degree3_scalar_20261101 import scalar_D3, symmetric_cubic_pair
from weight1_structure import J, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]


def quadratic_tensor_coefficients(eta, y, coordinates):
    """Return eta_ij,ab with eta_ij(y)=sum_ab eta_ij,ab*y_a*y_b."""
    return s.Matrix([s.Poly(eta[i,j],*y).coeff_monomial(y[a]*y[b])/(1 if a==b else 2)
                     for i,j,a,b in coordinates])


def tensor_pairing_matrix(coordinates):
    """Full ordered-index contraction eta_ij,ab xi_kl,cd Jik Jjl Jac Jbd."""
    result = s.zeros(len(coordinates))
    mate = {i:next((k,J[i,k]) for k in range(4) if J[i,k]) for i in range(4)}
    for row,(i,j,a,b) in enumerate(coordinates):
        k,ei = mate[i]; l,ej = mate[j]
        c,ea = mate[a]; d,eb = mate[b]
        sign = 1
        if k>l:
            k,l = l,k
            sign = -1
        c,d = sorted((c,d))
        column = coordinates.index((k,l,c,d))
        result[row,column] = 2*(1 if a==b else 2)*ei*ej*ea*eb*sign
    assert result==result.T
    return result


def verify():
    started = time.monotonic()
    model = S61()
    z,u,beta = primitive_beta(model)
    y = s.Matrix(z)
    two = [i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    q = s.symbols('delta_q0:45')
    values = [q[two.index(i)] if i in two else 0 for i in range(61)]
    quadratic_bivector = affine(model.reconstruct(values),y)
    eta_delta = -J*quadratic_bivector*J
    coordinates = [(i,j,a,b) for i,j in combinations(range(4),2)
                   for a,b in combinations_with_replacement(range(4),2)]
    embedding = quadratic_tensor_coefficients(eta_delta,y,coordinates).jacobian(q)
    full_pairing = tensor_pairing_matrix(coordinates)
    metric = embedding.T*full_pairing*embedding
    metric_determinant = metric.det()
    assert metric_determinant != 0

    xfield = -J*beta*y
    derivative = xfield.jacobian(y)
    eta1 = (derivative.T*J*derivative).applyfunc(s.expand)
    pfaffian = s.expand(beta[0,1]*beta[2,3]-beta[0,2]*beta[1,3]+beta[0,3]*beta[1,2])
    eta2 = exterior_derivative(-pfaffian*J*y,y)
    tensors = [quadratic_tensor_coefficients(eta,y,coordinates) for eta in (eta1,eta2)]
    section_path = ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json'
    section = json.loads(section_path.read_text())
    for tensor,stored in zip(tensors,section['covariants_in_q_coordinates'][:2],strict=True):
        covariant = s.Matrix(list(map(s.sympify,stored)))
        assert (embedding*covariant-tensor).applyfunc(s.expand)==s.zeros(60,1)
    functional = ((-18*tensors[0]-9*tensors[1]).T*full_pairing*embedding).applyfunc(s.expand)
    assert all(s.Poly(value,*u).total_degree()==2 for value in functional if value)

    plane_path = ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json'
    plane = json.loads(plane_path.read_text())
    assert list(map(str,u))==plane['beta_variables']
    nmatrix = s.Matrix(plane['q_linear_matrix_N'])
    center = s.Matrix(list(map(s.sympify,plane['q_center'])))
    offset = s.expand((functional*center)[0])
    c = s.symbols('c0:20')
    d = scalar_D3([beta.diff(v) for v in z])
    scalar_row = s.Matrix([symmetric_cubic_pair(d,s.Matrix(c))]).jacobian(c)
    residual = (functional*nmatrix-scalar_row).applyfunc(s.expand)
    assert all(s.Poly(value,*u).total_degree()==3 for value in residual if value)

    u_path = ROOT/'reports/TORIC_DEGREE5_COVARIANT_20261101.json'
    slice_path = ROOT/'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json'
    v_path = ROOT/'reports/TORIC_JACOBIAN_TWO_TORSION_20261101.json'
    u_report = json.loads(u_path.read_text())
    slice_report = json.loads(slice_path.read_text())
    v_report = json.loads(v_path.read_text())
    slice_matrix = s.Matrix(slice_report['slice_matrix_b'])
    assert slice_matrix==s.Matrix(v_report['slice_matrix_b'])
    basis = s.Matrix(u_report['weight_basis_24_by16'])
    assert basis==s.Matrix(slice_report['weight_basis_24_by16'])
    point = slice_matrix*s.Matrix([1,0,0,0,0,0])
    rows = list(basis.T.rref()[1])
    inverse = basis[rows,:].inv()
    tangent = s.Matrix.hstack(*[inverse*infinitesimal_action(basis*point,a)[rows,0]
                               for a in sp_basis()],slice_matrix)
    tangent_determinant = tangent.det()
    assert tangent_determinant != 0

    polynomial_ring,*variables = ring('s0,s1,s2,s3,s4,s5',QQ)
    source_map = s.Matrix(u_report['b_in_u']).inv()*slice_matrix
    images = [sum((polynomial_ring(QQ.from_sympy(source_map[i,j]))*variables[j]
                   for j in range(6) if source_map[i,j]),polynomial_ring.zero) for i in range(16)]
    powers = [[value**degree for degree in range(4)] for value in images]
    def restrict(expression):
        answer = polynomial_ring.zero
        for monomial,coefficient in s.Poly(expression,*u).terms():
            term = polynomial_ring(QQ.from_sympy(coefficient))
            for i,degree in enumerate(monomial):
                if degree:
                    term *= powers[i][degree]
            answer += term
        return answer
    def decode(encoded):
        return polynomial_ring.from_dict({tuple(mon):QQ.from_sympy(s.Rational(co)) for mon,co in encoded})
    restricted = [restrict(value) for value in residual]
    U = list(map(decode,slice_report['C5_sparse']))
    V = list(map(decode,v_report['C9_sparse']))
    assert not sum((a*b for a,b in zip(restricted,U,strict=True)),polynomial_ring.zero)
    assert not sum((a*b for a,b in zip(restricted,V,strict=True)),polynomial_ring.zero)

    paths = [Path(__file__),section_path,plane_path,u_path,slice_path,v_path,
             ROOT/'code/s61_model.py',ROOT/'code/q51_universal_section_20261101.py',
             ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
             ROOT/'code/toric_degree3_scalar_20261101.py',ROOT/'code/weight1_structure.py']
    paths += [model.source/name for name in ('symplectic_slice_certificate.json',
              'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json')]
    return {'id':'TORIC_POLYNOMIAL_E2_COORDINATE_20261101','all_passed':True,
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'first_jet_variables':list(map(str,u)),'weight2_u_indices_zero_based':two,
        'tensor_coordinates':[list(item) for item in coordinates],
        'tensor_pairing_convention':'sum_all_ordered_indices eta_ij,ab xi_kl,cd J_ik J_jl J_ac J_bd; symmetric coefficient is half the off-diagonal polynomial coefficient',
        'K45_pairing_matrix':[[str(v) for v in row] for row in metric.tolist()],
        'K45_pairing_determinant':str(metric_determinant),
        'two_covariant_form_identities_checked':True,'contraction_constants':[-18,-9],
        'lambda_coefficients_q':list(map(str,functional)),
        'q_center':list(map(str,center)),'lambda_center':str(offset),
        'polynomial_X_formula':'sum_i lambda_i(beta)*(q_i-q_center_i(beta))',
        'polynomial_weight':4,'polynomial_degree_in_q':1,
        'residual_row_lambda_N_minus_D':list(map(str,residual)),
        'residual_full16_monomial_count':sum(len(s.Poly(f,*u).terms()) for f in residual),
        'slice_matrix_b':slice_report['slice_matrix_b'],
        'dominance_tangent_determinant':str(tangent_determinant),
        'residual_U_contraction_zero_on_dominant_slice':True,
        'residual_V_contraction_zero_on_dominant_slice':True,
        'extension_argument':'Both scalar functions are invariant on the generic toric cover. Their equality on the dominant six-slice cover extends by Sp4 covariance and density, then by reduced closure.',
        'scope':'The polynomial is defined on S61, but its cubic E2 equation is asserted only on the reduced T01 closure. No full-fiber containment or exhaustion.',
        'proof_note':'research/notes/TORIC_POLYNOMIAL_E2_COORDINATE_20261101.md',
        'elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    report = verify()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as handle:
        json.dump(report,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'weight':4,'degree_in_q':1,
                      'seconds':report['elapsed_seconds'],'output':str(args.output)}))


if __name__=='__main__':
    main()
