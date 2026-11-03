#!/usr/bin/env python3
"""Actual full three-parameter Poisson family; canonical GL image rank26.

The former auxiliary N^2 equation is NOT imposed: every Jacobi equation
vanishes for all three parameters, including n0*n2!=0.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
X=s.Matrix(s.symbols('h x y z w'));PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in mon))for i,j in PAIRS for mon in MONS])
def matrix(M):return [[str(v)for v in row]for row in M.tolist()]
def canonical(P):
    div=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])
    return (P-(X*div.T-div*X.T)/5).applyfunc(s.expand)
def lie(P,i,j):
    A=s.zeros(5);A[i,j]=1;field=A*X
    return (A*P+P*A.T-sum((field[k]*P.diff(X[k])for k in range(5)),s.zeros(5))).applyfunc(s.expand)
def build(out):
    start=time.monotonic();ns=s.symbols('n0:3');h,x,y,z,w=X
    parent=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz'
    data=json.loads(gzip.decompress(parent.read_bytes()));point=list(map(s.sympify,data['exact99_polynomial_graph']))
    av=s.Matrix([x,y,z,w]);pairs=list(itertools.combinations(range(4),2));mons=list(map(s.sympify,data['coefficient_monomials']))
    affine=s.Matrix([[s.sympify(v)for v in row]for row in data['positive_P']])
    for k,(i,j)in enumerate(pairs):
        delta=sum(c*m for c,m in zip(point[15*k:15*k+15],mons));affine[i,j]+=delta;affine[j,i]-=delta
    affine=affine.applyfunc(s.expand);raw=s.zeros(5)
    all140=[];mon3=[s.prod(av[i]for i in ij)for d in range(4)for ij in itertools.combinations_with_replacement(range(4),d)]
    for i,j,k in itertools.combinations(range(4),3):
        expression=s.expand(sum(affine[i,b]*s.diff(affine[j,k],av[b])+affine[j,b]*s.diff(affine[k,i],av[b])+affine[k,b]*s.diff(affine[i,j],av[b])for b in range(4)))
        polynomial=s.Poly(expression,*av);all140 +=[polynomial.coeff_monomial(m)for m in mon3]
    assert len(all140)==140 and not any(all140)
    for i,j in pairs:
        polynomial=s.Poly(affine[i,j],*av)
        value=sum(c*h**(2-sum(exponent))*s.prod(X[k+1]**a for k,a in enumerate(exponent))for exponent,c in polynomial.terms())
        raw[i+1,j+1]=value;raw[j+1,i+1]=-value
    canon=canonical(raw);coefficients=vector(canon);assert len(coefficients)==150
    all350=[]
    for quad in itertools.combinations(range(5),4):
        expression=0
        for pos,e in enumerate(quad):
            i,j,k=quad[:pos]+quad[pos+1:]
            jac=sum(canon[i,b]*s.diff(canon[j,k],X[b])+canon[j,b]*s.diff(canon[k,i],X[b])+canon[k,b]*s.diff(canon[i,j],X[b])for b in range(5))
            expression+=(-1)**pos*X[e]*jac
        polynomial=s.Poly(s.expand(expression),*X)
        all350 +=[polynomial.coeff_monomial(s.prod(X[i]for i in mon))for mon in itertools.combinations_with_replacement(range(5),4)]
    assert len(all350)==350 and not any(all350)
    assert s.Matrix([s.expand(sum(s.diff(canon[i,j],X[i])for i in range(5)))for j in range(5)])==s.zeros(5,1)
    # No auxiliary incidence equations restrict the Poisson family.
    assert s.sympify(data['all209_exact_full_pullbacks'][179])==s.Rational(4,3)*ns[0]**2*ns[2]
    pf=s.factor(affine[0,1]*affine[2,3]-affine[0,2]*affine[1,3]+affine[0,3]*affine[1,2])
    assert pf!=0
    rank_records={}
    for point_n in [(0,0,0),(1,1,1)]:
        subs=dict(zip(ns,point_n));P=canon.subs(subs);v=coefficients.subs(subs)
        O=s.Matrix.hstack(*[vector(lie(P,i,j))for i,j in itertools.product(range(5),repeat=2)])
        M=O.row_join(v).row_join(coefficients.jacobian(ns).subs(subs));columns=list(M.rref()[1]);rows=list(M[:,columns].T.rref()[1])
        assert len(columns)==26
        det=s.Rational(DomainMatrix.from_Matrix(M[rows,columns]).convert_to(s.QQ).det());assert det
        rank_records[str(point_n)]=dict(full150_by29_matrix=matrix(M),rank=26,minor_rows=rows,minor_columns=columns,determinant=str(det))
    # Generic anticanonical divisor is3H+Q, with Q smooth and H tangent.
    homogeneous_pf=s.expand(h**4*pf.subs({x:x/h,y:y/h,z:z/h,w:w/h},simultaneous=True))
    quad=s.cancel(homogeneous_pf/h**2);assert s.Poly(quad,*X).total_degree()==2
    quad_matrix=s.hessian(quad,X)/2;generic=dict(zip(ns,[1,1,1]));assert quad_matrix.subs(generic).det()!=0
    assert quad_matrix[1:,1:].subs(generic).rank()==3
    out.mkdir(parents=True,exist_ok=False)
    certificate=dict(parameters=list(map(str,ns)),affine_coordinates=list(map(str,av)),full_affine_P=matrix(affine),full_raw_homogeneous_P=matrix(raw),
      full_canonical_homogeneous_P=matrix(canon),canonical150=list(map(str,coefficients)),all140_affine_Jacobi=list(map(str,all140)),all350_projective_Jacobi=list(map(str,all350)),
      auxiliary_modular_square_residual_not_imposed='4*n0^2*n2/3',rank_certificates=rank_records,
      anticanonical_quintic='h^3*Q',Q=str(quad),Q_symmetric_matrix=matrix(quad_matrix),Q_generic_det=str(quad_matrix.subs(generic).det()),H_section_generic_rank=3)
    archive=out/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63.md',parent,
      ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/report.json',ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63/report.json',
      ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json']
    result=dict(id='Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63',all_passed=True,field='QQ[n0,n1,n2]',
      actual_full_three_parameter_Poisson_family=True,auxiliary_N_squared_equation_explicitly_not_imposed=True,
      all140_affine_and_all350_projective_equations_zero_exactly=True,full150_canonical_coefficients_saved=True,
      full_GL25_P_scalar_three_parameter_image_rank=26,irreducible_projective_image_dimension_lower_bound=25,
      generic_Poisson_rank4=True,positive_point_occurs_at_parameter_origin=True,
      generic_actual_divisor='3H+Q, Q smooth, H tangent to Q',
      known_component_distinctness='All other known rank4 components have dimension<=24 except T01,T39,T40,CA,CB; those five closed components exclude the positive parameter origin by the pinned packets. Rank2 components cannot contain generic rank4 family.',
      at_least_one_new_full_P4_component=True,global_known_component_lower_bound_after_independent_audit=58,
      exact_new_component_dimension_or_generic_reducedness_claimed=False,
      firstjet_dominance_or_general_fiber_exhaustion_claimed=False,
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},certificate_sha256=sha(archive),seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ACTUAL3-PARAM FAMILY FULL140/350/150 +RANK26 NEW-COMPONENT LOWER PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);build(parser.parse_args().output)
