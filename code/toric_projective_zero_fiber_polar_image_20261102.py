#!/usr/bin/env python3
"""Prove that the discriminant polar directions dominate projective AB."""
import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from ab_ideal_high_generator_vanishing_20261102 import BASIS, MONOMIALS, PAIRS, bivector_action
from s61_model import S61, X
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]
STEM = 'TORIC_PROJECTIVE_ZERO_FIBER_POLAR_IMAGE_20261102'
P = 32749


def read(name):
    result = json.loads((ROOT/'reports'/(name+'.json')).read_text())
    assert result['all_passed']
    return result


def reduce_coefficient(value):
    c = Fraction(str(value))
    assert c.denominator % P
    return c.numerator * pow(c.denominator,-1,P) % P


def reduce_matrix(matrix):
    return np.array([[reduce_coefficient(x) for x in row] for row in matrix.tolist()],dtype=np.int64)


def value_and_gradient(polynomial, point):
    value = np.zeros(7,dtype=np.int64)
    assert all(point)
    for exponents, coefficient in polynomial:
        m = reduce_coefficient(coefficient)
        for s,power in zip(point,exponents,strict=True):
            m = m*pow(int(s),power,P) % P
        value[0] = (value[0]+m) % P
        for i,power in enumerate(exponents):
            value[1+i] = (value[1+i]+m*power*pow(int(point[i]),-1,P)) % P
    return value


def pivots(matrix):
    reduced = matrix.copy() % P
    rows = list(range(len(reduced)))
    row_indices = []
    columns = []
    for i in range(reduced.shape[1]):
        j = len(columns)
        if j == len(reduced):
            break
        possible = np.flatnonzero(reduced[j:,i])
        if not len(possible):
            continue
        k = j+int(possible[0])
        reduced[[j,k]] = reduced[[k,j]]
        rows[j],rows[k] = rows[k],rows[j]
        row_indices.append(rows[j])
        columns.append(i)
        reduced[j] = reduced[j]*pow(int(reduced[j,i]),-1,P) % P
        reduced[j+1:] = (reduced[j+1:]-reduced[j+1:,i,None]*reduced[j]) % P
    return row_indices,columns


def raw_quadratic_coordinates():
    model = S61()
    q = sp.symbols('q0:45')
    q_iter = iter(q)
    representative = model.reconstruct([next(q_iter) if weight == 2 else 0
                                        for weight in model.reduced['parameter_weights']])
    affine = representative.subs(X[0],1)
    y = X[1:]
    rows = []
    for i,j in PAIRS:
        comp = affine[1+i,1+j] - y[i]*affine[0,1+j] + y[j]*affine[0,1+i]
        poly = sp.Poly(comp,*y)
        rows.extend(poly.coeff_monomial(y[r]*y[t]) for r,t in MONOMIALS)
    assert len(rows) == len(BASIS) == 60
    E = sp.Matrix(rows).jacobian(q)
    assert E.rank() == 45
    return reduce_matrix(E)


def verify():
    start = time.monotonic()
    assert sp.isprime(P)
    dependencies = [
        'TORIC_ZERO_JET_ASSOCIATIVE_LIFTS_20261102',
        'TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_20261102',
        'TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_REVIEW_20261102',
        'PRIMITIVE_CURVE_DOMINANT_SLICE_20261101','PRIMITIVE_CURVE_DISCRIMINANT_20261101',
        'TORIC_DEGREE5_COVARIANT_20261101','TORIC_JACOBIAN_TWO_TORSION_20261101',
        'TORIC_UNIVERSAL_THREE_PLANE_20261101','TORIC_C9_POLYNOMIAL_20261101',
        'AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102','AB_INCIDENCE_20261101']
    reports = {name:read(name) for name in dependencies}
    lifts = reports[dependencies[0]]
    assert lifts['exact_row_identities'] == len(lifts['associative_lifts']) == 274
    for i,record in enumerate(lifts['associative_lifts']):
        assert i == record['quadratic_index']
        assert all(len(term) == 4 and 0 <= term[0] < 45 and 0 <= term[1] <= term[2] < 16
                   for term in record['beta2_z_terms'])
        assert all(len(term) == 5 and all(0<=j<16 for j in term[:4])
                   for term in record['beta4_terms'])
    slice_data = reports['PRIMITIVE_CURVE_DOMINANT_SLICE_20261101']
    delta = reports['PRIMITIVE_CURVE_DISCRIMINANT_20261101']
    U = reports['TORIC_DEGREE5_COVARIANT_20261101']
    V = reports['TORIC_JACOBIAN_TWO_TORSION_20261101']
    plane = reports['TORIC_UNIVERSAL_THREE_PLANE_20261101']
    assert slice_data['slice_matrix_b'] == V['slice_matrix_b']
    assert slice_data['weight_basis_24_by16'] == U['weight_basis_24_by16']
    assert U['beta_variables'] == plane['beta_variables']
    assert U['c_variables'] == plane['c_variables']
    base = np.array(delta['slice_line_base'])
    line = np.array(delta['slice_line_direction'])
    t0 = 4577
    s = (base+t0*line) % P
    c4,c6 = [value_and_gradient(slice_data[k],s) for k in ['c4_sparse','c6_sparse']]
    u = np.array([value_and_gradient(f,s) for f in slice_data['C5_sparse']])
    v = np.array([value_and_gradient(f,s) for f in V['C9_sparse']])
    assert u.shape == v.shape == (20,7)
    C4,C6 = int(c4[0]),int(c6[0])
    assert C4 != 0 and (C4**3-C6**2) % P == 0
    discriminant_gradient = (3*pow(C4,2,P)*c4[1:]-2*C6*c6[1:]) % P
    slope = int(discriminant_gradient@line) % P
    assert slope == 19085
    # The cone representative gamma=N(c4 V-3c6 U). Columns 1..6 are derivatives on the slice.
    w = (C4*v[:,0]-3*C6*u[:,0]) % P
    dw = (v[:,0,None]*c4[1:]+C4*v[:,1:]-3*u[:,0,None]*c6[1:]-3*C6*u[:,1:]) % P
    parameters = list(map(sp.Symbol,plane['beta_variables']))
    N = sp.Matrix(plane['q_linear_matrix_N']).applyfunc(sp.sympify)
    old_from_slice = sp.Matrix(U['b_in_u']).inv()*sp.Matrix(slice_data['slice_matrix_b'])
    N_derivatives = np.array([reduce_matrix(N.diff(b)) for b in parameters])
    dN = np.einsum('frj,fn->nrj',N_derivatives,reduce_matrix(old_from_slice)) % P
    N0 = np.einsum('nrj,n->rj',dN,s) % P
    gamma = N0@w % P
    derivative = (N0@dw+np.einsum('nrj,j->rn',dN,w)) % P
    # A polynomial/rational identity on the whole discriminant must hold at this smooth chart.
    for relation in reports['AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102']['quadratic_relations']:
        assert sum(reduce_coefficient(c)*int(gamma[i])*int(gamma[j]) for i,j,c in relation) % P == 0
    pivot = int(np.flatnonzero(discriminant_gradient)[0])
    tangent = np.eye(6,dtype=np.int64)
    tangent[pivot] = -discriminant_gradient * pow(int(discriminant_gradient[pivot]),-1,P) % P
    tangent = np.delete(tangent,pivot,axis=1)
    assert np.all(discriminant_gradient@tangent % P == 0)
    E = raw_quadratic_coordinates()
    raw = E@gamma % P
    orbit = []
    for a,b in combinations_with_replacement(range(4),2):
        symmetric = sp.zeros(4)
        symmetric[a,b] = symmetric[b,a] = 1
        A = J*symmetric
        assert A*J+J*A.T == sp.zeros(4)
        orbit.append(reduce_matrix(bivector_action(A))@raw % P)
    joint = np.column_stack([raw,*orbit,E@derivative@tangent % P])
    selected,columns = pivots(joint)
    assert len(selected) == len(columns) == 10
    minor = [list(map(int,joint[i,columns])) for i in selected]
    determinant = int(sp.Matrix(minor).det()) % P
    assert determinant
    orbit_rank = len(pivots(joint[:,1:11])[0])
    orbit_and_radial = len(pivots(joint[:,:11])[0])
    assert (orbit_rank,orbit_and_radial) == (8,9)
    # Derivatives of the six raw quadratic coefficients determine constant translation invariance.
    hess = np.zeros((6,4,4),dtype=np.int64)
    for index,((i,j),(a,b)) in enumerate(BASIS):
        h = PAIRS.index((i,j))
        hess[h,a,b] += raw[index]
        hess[h,b,a] += raw[index]
    divergence = np.zeros((4,4),dtype=np.int64)
    for h,(i,j) in enumerate(PAIRS):
        divergence[j] += hess[h,i]
        divergence[i] -= hess[h,j]
    assert len(pivots(hess.reshape((24,4)))[0]) == 2
    assert len(pivots(divergence)[0]) == 2 and np.all(divergence@divergence % P == 0)
    files = [Path(__file__),ROOT/'research/notes'/(STEM+'.md')]
    files.extend(ROOT/'reports'/(n+'.json') for n in dependencies)
    files.extend(ROOT/'research/notes'/(n+'.md') for n in dependencies)
    files.extend(ROOT/'code'/n for n in ['s61_model.py','weight1_structure.py',
                                        'ab_ideal_high_generator_vanishing_20261102.py'])
    files.extend(ROOT/'legacy/v02/data'/n for n in ['symplectic_slice_certificate.json',
                 'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    return dict(id=STEM,all_passed=True,
                polar_numerator='N(beta)*(c4(beta)*V(beta)-3*c6(beta)*U(beta))',
                numerator_degree=18,projective_image_closure='P(AB)',projective_image_dimension=9,
                generic_fiber_dimension_from_affine_discriminant=6,
                projective_zero_fiber='Proj C[h,q]/I_AB(q)',
                whole_boundary_closed_subscheme_of='W16 x P(AB)',
                probe_field=P,line_base=base.tolist(),line_direction=line.tolist(),line_t=t0,
                slice_parameters=s.tolist(),c4=C4,c6=C6,
                gradient_c4_cubed_minus_c6_squared=discriminant_gradient.tolist(),
                line_simple_root_derivative=slope,old_q45_direction=gamma.tolist(),
                tangent_basis=tangent.tolist(),joint_minor_rows=selected,joint_minor_columns=columns,
                joint_minor=minor,joint_minor_determinant_residue=determinant,
                joint_rank_at_reduction=10,orbit_rank_at_reduction=orbit_rank,
                orbit_and_radial_rank_at_reduction=orbit_and_radial,
                translation_kernel_dimension_at_reduction=2,modular_rank_at_reduction=2,
                all274_relations_verified_at_reduction=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in files},seconds=time.monotonic()-start)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS',report['joint_minor_determinant_residue'],round(report['seconds'],3))
