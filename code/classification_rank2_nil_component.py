#!/usr/bin/env python3
"""Whole ten-parameter (3,1), c=0, n7!=0 source and one full tangent test."""
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from classification_double_quadric import full_tangent_ranks, check_affine_poisson


def run():
    start = time.monotonic()
    source_path = Path('reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json')
    source = json.loads(source_path.read_text())['cases'][0]
    nn = s.symbols('n0:15')
    xx4 = s.symbols('x1:5')
    local = {str(t): t for t in (*nn, *xx4)}
    n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14 = nn
    sub = {n11:0,n9:0,n14:-1-2*n10,
           n13:-s.Rational(2,3)*(3*n10+1)*n6/n7,
           n12:-(2*n10+1)*n5/n7}
    parameters = (n0,n1,n2,n3,n4,n8,n10,n5,n6,n7)
    for expression in source['all_distinct_jacobi_coefficients']:
        assert s.factor(s.sympify(expression, locals=local).subs(sub)) == 0
    brackets = [s.factor(s.sympify(f, locals=local).subs(sub)) for f in source['full_six_brackets']]
    pfaffian = s.factor(brackets[0]*brackets[5]-brackets[1]*brackets[4]+brackets[2]*brackets[3])
    expected = -xx4[2]**2*(n5+s.Rational(2,3)*n6*xx4[2]+n7*xx4[3])/2
    assert s.factor(pfaffian-expected) == 0
    h = s.symbols('h')
    xx = (h,*xx4)
    P = s.zeros(5)
    for (i,j), f in zip(itertools.combinations(range(1,5),2), brackets):
        P[i,j] = s.expand(6*sum(a*h**(2-sum(e))*s.prod(x**v for x,v in zip(xx4,e))
                               for e,a in s.Poly(f,*xx4).terms()))
        P[j,i] = -P[i,j]
    point_sub = {p:1 for p in parameters}
    point = P.subs(point_sub)
    check_affine_poisson(point, xx)
    directions = [P.diff(p).subs(point_sub) for p in parameters]
    for tensor in (point,*directions):
        assert all(c.q==1 for f in tensor for c in s.Poly(f,*xx).coeffs())
    result = full_tangent_ranks(point, xx, directions)
    result.update(
        status='PASS', source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_report_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        full_source_parameters=list(map(str,parameters)),
        full_source_substitution={str(k):str(v) for k,v in sub.items()},
        all_parameter_jacobi_verified=True, all_parameter_pfaffian=str(pfaffian),
        universal_brackets={f'{i},{j}':str(P[i,j]) for i,j in itertools.combinations(range(1,5),2)},
        witness={str(k):str(v) for k,v in point_sub.items()},
        witness_brackets={f'{i},{j}':str(point[i,j]) for i,j in itertools.combinations(range(1,5),2)},
        localization='n7 != 0',
        complements='Other n7=0 branches and other modular types are retained separately.',
        seconds=round(time.monotonic()-start,4),
    )
    return result


if __name__ == '__main__':
    report=run()
    out=Path('reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_20261107.json')
    out.write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({k:report[k] for k in ('status','matrix_ranks_mod_prime',
        'projective_tangent_dimension_upper','projective_family_dimension_lower','seconds')}))
