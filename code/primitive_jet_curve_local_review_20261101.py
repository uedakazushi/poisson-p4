#!/usr/bin/env python3
"""Independent tensor/Pfaffian and local infinity-dimension certificate review.

The curve construction is rebuilt from closed/primitive tensor constraints,
without importing the original curve or flag-obstruction implementations.
The local homogeneous chart is rebuilt directly from the archived equations.
Geometric generic smoothness and dimension arguments are in the review note.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import sympy as s
from s61_model import S61, X
from elliptic_weight1_dominance import projective_affine

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT/'reports/S61_FLAG_LOCAL_DIMENSION_20261101.json'
EXTRA = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'


def pf4(a):
    return s.expand(a[0,1]*a[2,3]-a[0,2]*a[1,3]+a[0,3]*a[1,2])


def curve_checks():
    pairs = list(combinations(range(4), 2))
    coeff = s.symbols('b0:24')
    y = s.symbols('y0:4')
    beta = s.zeros(4)
    for h, (i, j) in enumerate(pairs):
        beta[i,j] = sum(coeff[4*h+a]*y[a] for a in range(4))
        beta[j,i] = -beta[i,j]
    relations = [s.diff(beta[j,k],y[i])-s.diff(beta[i,k],y[j])+s.diff(beta[i,j],y[k])
                 for i,j,k in combinations(range(4),3)]
    relations += [s.diff(beta[0,1]+beta[2,3],v) for v in y]
    constraint = s.Matrix(relations).jacobian(coeff)
    basis = s.Matrix.hstack(*constraint.nullspace())
    assert constraint.rank() == 8 and basis.shape == (24,16)
    a,b,c,d,e = p = s.symbols('a b c d e')
    plucker = s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,-a],[-c,-e,a,0]])
    dual = s.Matrix([[0,-a,-e,d],[a,0,c,-b],[e,-c,0,a],[-d,b,-a,0]])
    quadric = -a*a-b*e+c*d
    assert plucker*dual == -quadric*s.eye(4)
    cov = s.Matrix([sum(s.diff(beta[i,j],y[h])*plucker[i,j] for i,j in pairs)
                    for h in range(4)])
    quadrics = (plucker*cov).applyfunc(s.expand)
    monomials = [m for m in combinations_with_replacement(range(5),2) if m != (0,0)]
    coefficients = s.Matrix([s.Poly(f.subs(a*a,-b*e+c*d),*p).coeff_monomial(p[i]*p[j])
                             for f in quadrics for i,j in monomials]).jacobian(coeff)
    assert (coefficients*basis).rank() == 16
    v, w = s.Matrix([0,-1,0,0]), s.Matrix([0,0,0,-1])
    direction = v*w.T-w*v.T
    flag_cov = s.Matrix([sum(s.diff(beta[i,j],y[h])*direction[i,j] for i,j in pairs)
                         for h in range(4)])
    evaluation_full = flag_cov.jacobian(coeff)*basis
    evaluation_plane = s.Matrix.vstack(v.T*evaluation_full,w.T*evaluation_full)
    assert evaluation_full.rank() == 4 and evaluation_plane.rank() == 2
    # Pfaffian identity for arbitrary covector: stronger than a sampled beta.
    arbitrary = s.Matrix(s.symbols('ell0:4'))
    matrix = dual.row_join(arbitrary).col_join((-arbitrary.T).row_join(s.zeros(1)))
    signed = s.Matrix([(-1)**i*pf4(matrix.minor_submatrix(i,i)) for i in range(5)])
    expected = (plucker*arbitrary).col_join(s.Matrix([quadric]))
    assert (signed-expected).applyfunc(s.expand) == s.zeros(5,1)
    assert (matrix*signed).applyfunc(s.expand) == s.zeros(5,1)
    return {'tensor_constraint_rank':8, 'closed_primitive_dimension':16,
            'section_map_rank':16, 'plane_evaluation_rank':2,
            'full_covector_evaluation_rank':4,
            'universal_Pfaffian_and_syzygy_identities':True}


def local_checks():
    model = S61()
    stored = json.loads(LOCAL.read_text())
    weights = model.reduced['parameter_weights']
    one = [i for i,w in enumerate(weights) if w == 1]
    two = [i for i,w in enumerate(weights) if w == 2]
    value = {int(i)-1:s.Rational(v) for i,v in stored['beta'].items()}
    k = two.index(37)
    columns = [i for i in range(45) if i != k]
    chart_index = {two[j]:i for i,j in enumerate(columns)}
    equations = model.reduced['equations']+[json.loads(EXTRA.read_text())['universal_new_equation']]
    chart = []
    for equation in equations:
        polynomial = {}
        for mon,n,d in equation:
            coefficient = s.Rational(n,d)
            qcount, newmon = 0, []
            for index in mon:
                if weights[index] == 1:
                    coefficient *= value[index]
                else:
                    qcount += 1
                    if index != 37:
                        newmon.append(chart_index[index])
            newmon += [44]*(2-qcount)
            key = tuple(sorted(newmon))
            polynomial[key] = polynomial.get(key,0)+coefficient
        assert polynomial.get((),0) == 0
        chart.append(polynomial)
    jacobian = s.Matrix([[f.get((j,),0) for j in range(45)] for f in chart])
    assert jacobian == s.Matrix(stored['homogenized_jacobian']).applyfunc(s.Rational)
    tangent = s.Matrix(stored['tangent_basis']).applyfunc(s.Rational)
    left = s.Matrix(stored['cokernel_rows']).applyfunc(s.Rational)
    assert len(jacobian.to_DM().rref()[1]) == 38
    assert tangent.rank() == 7 and tangent[44,:] == s.zeros(1,7)
    assert jacobian*tangent == s.zeros(135,7)
    assert left*jacobian == s.zeros(10,45)
    monomials = list(combinations_with_replacement(range(7),2))
    quadratic_rows = []
    for polynomial in chart:
        row = {}
        for mon,coefficient in polynomial.items():
            if len(mon) != 2:
                continue
            for i in range(7):
                for j in range(7):
                    key = tuple(sorted((i,j)))
                    row[key] = row.get(key,0)+coefficient*tangent[mon[0],i]*tangent[mon[1],j]
        quadratic_rows.append([row.get(mon,0) for mon in monomials])
    obstructions = (left.to_DM()*s.Matrix(quadratic_rows).to_DM()).to_Matrix()
    assert obstructions == s.Matrix(stored['quadratic_obstruction_basis']).applyfunc(s.Rational)
    assert obstructions.rank() == 10
    hessians = []
    for row in obstructions.tolist():
        f = dict(zip(monomials,row,strict=True))
        hessians += [[f[tuple(sorted((i,j)))]*(2 if i == j else 1)
                      for j in range(7)] for i in range(7)]
    hessian = s.Matrix(hessians)
    assert hessian.rank() == 4
    radical = s.Matrix(stored['quadratic_common_radical']).applyfunc(s.Rational)
    assert radical.rank() == 3 and hessian*radical == s.zeros(70,3)
    # Reconstruct the first jet independently and check the smooth flag lower bound.
    params = [value.get(i,s.S.Zero) for i in range(61)]
    affine_y = s.symbols('y0:4')
    homogeneous = model.reconstruct(params)
    affine = projective_affine(homogeneous,X,affine_y)
    origin = dict.fromkeys(affine_y,0)
    linear = sum((affine.diff(x).subs(origin)*x for x in affine_y),s.zeros(4))
    standard = s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]))
    alpha = s.Matrix([s.diff(linear[0,1]+linear[2,3],x)/5 for x in affine_y])
    yy = s.Matrix(affine_y)
    vv = standard*alpha
    primitive = -standard*(linear-2*alpha.dot(yy)*standard+yy*vv.T-vv*yy.T)*standard
    v = s.Matrix([0,-1,0,0])
    wprime = s.Matrix([0,0,1,0])
    incidence_derivative = (v.T*primitive.subs(dict(zip(affine_y,v,strict=True)))*wprime)[0]
    assert incidence_derivative != 0
    return {'homogeneous_chart_rebuilt_from_all135_equations':True,
            'jacobian_rank':38, 'tangent_dimension':7, 'tangent_t_coordinate_zero':True,
            'independent_quadratic_obstruction_count':10, 'common_hessian_rank':4,
            'common_radical_dimension':3,
            'flag_incidence_derivative':str(incidence_derivative),
            'local_dimension_bounds':[3,3],
            'scope':'One saved rational first jet and one point of the unsaturated135-generator homogenization only.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = {'id':'PRIMITIVE_JET_CURVE_LOCAL_REVIEW_20261101',
              'generated_utc':datetime.now(timezone.utc).isoformat(),
              'all_passed':True, 'curve':curve_checks(), 'local':local_checks(),
              'proof_note':'research/notes/PRIMITIVE_JET_CURVE_LOCAL_REVIEW_20261101.md',
              'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__).resolve(),LOCAL,EXTRA,
                                        ROOT/'code/s61_model.py']}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS: independent curve/Pfaffian construction and local unsaturated dimension3 review.')


if __name__ == '__main__':
    main()
