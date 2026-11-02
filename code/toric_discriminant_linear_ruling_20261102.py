#!/usr/bin/env python3
"""Compute the discriminant ruling and the critical first jet in an adapted symplectic frame."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61

ROOT = Path(__file__).resolve().parents[1]
STEM = 'TORIC_DISCRIMINANT_LINEAR_RULING_20261102'


def ambient_first_derivative(model, coefficients, beta_positions):
    """Differentiate the original equations once in beta and restrict to beta=0."""
    result = sp.zeros(len(model.reduced['equations']), len(beta_positions))
    index = {position: column for column, position in enumerate(beta_positions)}
    for row, equation in enumerate(model.reduced['equations']):
        for monomial, numerator, denominator in equation:
            beta_indices = [i for i in monomial if i in index]
            if len(beta_indices) == 1:
                coefficient = sp.Rational(numerator, denominator)
                coefficient *= sp.prod(coefficients[i] for i in monomial if i not in index)
                result[row, index[beta_indices[0]]] += coefficient
    return result


def verify():
    model = S61()
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    position_variables, u, beta = primitive_beta(model)
    u4,u5,u12,u13,u14,u25,u26,u27,u34,u35,u36,u37,u48,u49,u50,u51 = u
    A,B,C,D,E = parameters = sp.symbols('A B C D E')

    # In the adapted frame Qxy=-Ax², Qzw=-Bz², Qyw=Cx²+Dxz+Ez².
    coefficients = [0]*61
    for index, value in [(0,-A),(37,C),(39,D),(44,E),(58,-B)]:
        coefficients[index] = value
    from s61_model import X
    tensor = model.reconstruct(coefficients).subs(X[0],1)
    y = X[1:]
    quad = sp.zeros(4)
    quad[0,1] = 1-A*y[0]**2
    quad[2,3] = 1-B*y[2]**2
    quad[1,3] = C*y[0]**2+D*y[0]*y[2]+E*y[2]**2
    quad = quad-quad.T
    for i in range(4):
        for j in range(i+1,4):
            actual = tensor[i+1,j+1]-y[i]*tensor[0,j+1]+y[j]*tensor[0,i+1]
            assert sp.expand(actual-quad[i,j]) == 0
    jacobian = ambient_first_derivative(model, coefficients, one)
    fixed = [u4,u5,u12,u26,u27,u49-2*u14]
    recovery = [2*A*u13+3*D*u14, 2*A*u35-3*C*u14,
                2*B*u25-D*u51, B*(2*u37+u48)+E*u51]
    expected = sp.Matrix(fixed+recovery).jacobian(u)
    assert expected.to_DM().rank() == 10
    assert jacobian.to_DM().rank() == 10
    assert jacobian.col_join(expected).to_DM().rank() == 10
    recovery_matrix = sp.Matrix(recovery).jacobian(parameters)
    determinant = sp.factor(recovery_matrix[:,[0,1,2,4]].det())
    assert determinant != 0
    representative = sp.Matrix([-3*u14/(2*u13),u51/(2*u25),-u35/u13,
                                1,-(2*u37+u48)/(2*u25)])
    assert (recovery_matrix*representative).applyfunc(sp.cancel) == sp.zeros(4,1)

    p,l,m = sp.symbols('p l m')
    left, right = sp.Matrix([p,1,l,0]), sp.Matrix([l,0,m,1])
    zero = {p:0,l:0,m:0}
    sections = sp.Matrix([(left.T*beta.subs(dict(zip(position_variables, point)))*right)[0]
                          for point in (left, right)])
    sections = sections.applyfunc(sp.expand)
    derivative = sections.jacobian((p,l,m)).subs(zero)
    conditions = list(sections.subs(zero))+list(derivative[:,0])+list(derivative[:,2])
    restricted = sp.Matrix(conditions).jacobian(u)
    assert restricted.to_DM().rank() == restricted.col_join(sp.Matrix(fixed).jacobian(u)).to_DM().rank() == 6
    specialize = {u4:0, u5:0, u12:0, u26:0, u27:0, u49:2*u14}
    sections = sections.subs(specialize)
    assert sections.jacobian((p,l,m)).subs(zero) == sp.Matrix([[0,3*u14,0],[0,u51,0]])
    transverse = u51*sections[0]-3*u14*sections[1]
    hessian = sp.hessian(transverse,(p,m)).subs(zero)
    sample = {value:i+1 for i,value in enumerate(u)}
    hessian_witness = hessian.subs(sample).det()
    assert hessian_witness == -78265
    sources = [Path(__file__), ROOT/'research/notes'/f'{STEM}.md',
               ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
               ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/s61_model.py']
    for name in ['TORIC_PROJECTIVE_ZERO_FIBER_POLAR_IMAGE_20261102',
                 'TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_20261102',
                 'TORIC_ZERO_FIBER_GENERIC_REDUCED_20261102',
                 'AB_INCIDENCE_20261101','PRIMITIVE_CURVE_DISCRIMINANT_20261101']:
        path = ROOT/'reports'/f'{name}.json'
        assert json.loads(path.read_text())['all_passed']
        sources.append(path)
    return dict(id=STEM, all_passed=True, field='QQ(A,B,C,D,E)',
                beta_variable_order=list(map(str,u)), fixed_critical_conditions=list(map(str,fixed)),
                recovery_equations=list(map(str,recovery)), recovery_minor=str(determinant),
                unique_projective_coefficients_normalized_D_equals_one=list(map(str,representative)),
                ambient_beta_jacobian_shape=[134,16], ambient_beta_jacobian_generic_rank=10,
                projective_ruling_dimension=5, polar_graph_birational=True,
                Hessian_kernel_plane=[[str(x) for x in row] for row in hessian.tolist()],
                ordinary_node_hessian_value=str(hessian_witness),
                generic_critical_planes=1, general_canonical_discriminant_curve='irreducible rational quintic with one ordinary node',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS', report['recovery_minor'], report['ordinary_node_hessian_value'])
