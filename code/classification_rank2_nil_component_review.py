#!/usr/bin/env python3
"""Independent full-source and full-ambient review of the nil(3,1) component."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as S

import classification_rank3_mixed_component_review as independent

ROOT=Path(__file__).resolve().parents[1]
PAIRS4=list(itertools.combinations(range(4),2))


def sparse(p, coordinates):
    return {ij:dict(S.Poly(p[ij], *coordinates).terms())
            for ij in independent.PAIRS if p[ij]!=0}


def skew(entries, size=4):
    p=S.zeros(size)
    for ij,entry in zip(itertools.combinations(range(size),2),entries):
        p[ij]=entry
        p[ij[::-1]]=-entry
    return p


def normalize(expr, parameters):
    p=S.Poly(S.expand(expr),*parameters)
    return S.expand(p.as_expr()/p.LC())


def main(output):
    began=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    signal.alarm(90)
    source_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
    component_path=ROOT/'reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_20261107.json'
    source=json.loads(source_path.read_text())['cases'][0]
    component=json.loads(component_path.read_text())
    h,*xx=S.symbols('h x1 x2 x3 x4')
    coords=(h,*xx)
    nn=S.symbols('n0:15')
    env={str(v):v for v in (*coords,*nn)}
    entries=[S.sympify(f,locals=env) for f in source['full_six_brackets']]
    p15=skew(entries)

    # Reconstruct the full affine-linear 90-variable system without producer helpers.
    mons=[S.Integer(1),*xx,*[xx[i]*xx[j] for i,j in itertools.combinations_with_replacement(range(4),2)]]
    unknowns=S.symbols('a0:90')
    raw=skew([sum(unknowns[15*i+j]*m for j,m in enumerate(mons)) for i in range(6)])
    N=S.zeros(4)
    N[0,1]=1
    N[1,2]=1
    field=N*S.Matrix(xx)
    eq=[sum(S.diff(raw[i,j],xx[i]) for i in range(4))-field[j] for j in range(4)]
    lie=sum((field[k]*raw.diff(xx[k]) for k in range(4)),S.zeros(4))-N*raw-raw*N.T
    eq.extend(lie[ij] for ij in PAIRS4)
    top=skew([sum(unknowns[15*i+j]*mons[j] for j in range(5,15)) for i in range(6)])
    eq.extend(top[i,j]*field[k]-top[i,k]*field[j]+top[j,k]*field[i]
              for i,j,k in itertools.combinations(range(4),3))
    coeffs=[c for expression in eq for c in S.Poly(S.expand(expression),*xx).coeffs() if c!=0]
    linear,rhs=S.linear_eq_to_matrix(coeffs,unknowns)
    source_vector=S.Matrix([S.Poly(p15[ij],*xx).coeff_monomial(mon) for ij in PAIRS4 for mon in mons])
    assert all(S.expand(v)==0 for v in linear*source_vector-rhs)
    linear_cert=independent.minor_certificate(linear,75)
    source15_cert=independent.minor_certificate(source_vector.jacobian(nn),15)

    # Recover necessity of all six residuals from independent raw Jacobi coefficients.
    homogeneous15=S.zeros(5)
    for (i,j),entry in zip(itertools.combinations(range(1,5),2),entries):
        f=S.expand(h*h*entry.subs({x:x/h for x in xx},simultaneous=True))
        homogeneous15[i,j]=f
        homogeneous15[j,i]=-f
    all_raw=independent.jacobi(sparse(homogeneous15,coords),raw=True)
    actual_coefficients={normalize(c,nn) for polynomial in all_raw.values()
                         for c in polynomial.values() if S.expand(c)!=0}
    residuals=[S.sympify(v,locals=env) for v in source['all_distinct_jacobi_coefficients']]
    assert len(residuals)==6
    assert all(normalize(v,nn) in actual_coefficients for v in residuals)
    sub={nn[11]:0,nn[9]:0,nn[14]:-1-2*nn[10],
         nn[13]:-S.Rational(2,3)*(3*nn[10]+1)*nn[6]/nn[7],
         nn[12]:-(2*nn[10]+1)*nn[5]/nn[7]}
    for v in residuals:
        assert S.cancel(v.subs(sub,simultaneous=True))==0
    # Triangular necessity on n7 != 0: n11,n9,n14,n13,n12, respectively.
    triangle=[nn[11],nn[9],nn[14],nn[13],nn[12]]
    residual_order=[1,2,0,3,4]
    cumulative={}
    for var,row in zip(triangle,residual_order):
        relation=S.cancel(residuals[row].subs(cumulative,simultaneous=True))
        answer=S.solve(relation,var)
        assert len(answer)==1 and S.cancel(answer[0]-sub[var])==0
        cumulative[var]=sub[var]
    family=(6*homogeneous15.subs(sub,simultaneous=True)).applyfunc(S.expand)
    assert all(S.cancel(c)==0 for polynomial in independent.jacobi(sparse(family,coords),raw=True).values()
               for c in polynomial.values())
    pf=S.factor(family[1,2]*family[3,4]-family[1,3]*family[2,4]+family[1,4]*family[2,3])
    expected=-18*h*xx[2]**2*(nn[5]*h+S.Rational(2,3)*nn[6]*xx[2]+nn[7]*xx[3])
    assert S.cancel(pf-expected)==0
    for ij in independent.PAIRS:
        key=f'{ij[0]},{ij[1]}'
        assert S.cancel(family[ij]-S.sympify(component['universal_brackets'].get(key,'0'),locals=env))==0

    parameters=[nn[i] for i in (0,1,2,3,4,8,10,5,6,7)]
    point={v:1 for v in parameters}
    actual=family.subs(point).applyfunc(S.expand)
    raw=sparse(actual,coords)
    assert independent.jacobi(raw)==S.zeros(350,1)
    cols=[]
    for ij,mon in independent.COLUMNS:
        altered={key:dict(val) for key,val in raw.items()}
        altered.setdefault(ij,{})[mon]=altered.get(ij,{}).get(mon,0)+1
        cols.append(independent.jacobi(altered))
    derivative=S.Matrix.hstack(*cols)
    ee=S.Matrix(coords)
    vertical=[]
    orbit=[]
    for i,j in itertools.product(range(5),repeat=2):
        A=S.zeros(5)
        A[i,j]=1
        v=A*ee
        vertical.append(independent.vector(sparse(ee*v.T-v*ee.T,coords)))
        orbit.append(independent.vector(sparse(A*actual+actual*A.T-coords[j]*actual.diff(coords[i]),coords)))
    vm=S.Matrix.hstack(*vertical)
    group=S.Matrix.hstack(*vertical,*orbit,independent.vector(raw))
    source10=S.Matrix.hstack(*[independent.vector(sparse(family.diff(v).subs(point),coords)) for v in parameters])
    full=group.row_join(source10)
    assert derivative*full==S.zeros(350,full.cols)
    certificates={label:independent.minor_certificate(mat,rank) for label,mat,rank in
                  [('D',derivative,100),('vertical',vm,24),('group',group,47),('full_actual_source',full,50)]}
    paths=[Path(__file__).resolve(),Path(independent.__file__).resolve(),source_path,component_path]
    result=dict(all_passed=True,task_id='E_multiple_hyperplane',producer_arithmetic_imported=False,
        full90_linear_system=dict(certificate=linear_cert,nullity=15,source15=source15_cert,
            all_source_coefficients_satisfy_affine_system=True),
        all_six_jacobi_residuals_independently_necessary=True,
        n7_nonzero_source_triangular_elimination=True,
        universal_full10_raw_jacobi_zero=True, universal_scaled_pfaffian=str(pf),
        all350_projective_jacobi_zero=True,all61_actual_columns_annihilated=True,
        certificates=certificates,exact_field='QQ',
        source10_sparse_columns=[{str(i):str(source10[i,j]) for i in range(150) if source10[i,j]}
                                 for j in range(source10.cols)],
        matrix_sha256={label:hashlib.sha256(json.dumps(mat.tolist(),default=str).encode()).hexdigest()
                       for label,mat in [('D',derivative),('group',group),('source10',source10)]},
        projective_tangent_dimension=25,projective_actual_family_rank=25,
        smooth_full_component_proved=True,full_n7_nonzero_source_covered=True,
        n7_zero_complement_retained=True,global_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.monotonic()-began)
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    signal.alarm(0)
    print(json.dumps(dict(all_passed=True,D=100,source=50,linear_nullity=15,seconds=result['seconds'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
