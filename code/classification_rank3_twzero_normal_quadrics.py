#!/usr/bin/env python3
"""True three-normal quadratic obstruction at the actual49-dimensional source."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse,vector

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_RANK3_TWZERO_GENERIC_REVIEW_20261107.json'


def matrix(encoded,x):
    P=s.zeros(5)
    for ij in itertools.combinations(range(5),2):
        P[ij]=s.sympify(encoded.get(str(ij),'0'));P[ij[::-1]]=-P[ij]
    return P


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();r=json.loads((ROOT/SOURCE).read_text())
    x=s.symbols('h u v z w');P=matrix(r['point_brackets'],x);ps=sparse(P,x)
    assert jacobi(ps)==s.zeros(350,1)
    dcols=[]
    for ij,e in COLUMNS:
        variation={key:dict(value)for key,value in ps.items()}
        variation.setdefault(ij,{})[e]=variation.get(ij,{}).get(e,0)+1
        dcols.append(jacobi(variation))
    D=s.Matrix.hstack(*dcols)
    cols=[];X=s.Matrix(x)
    for a,b in itertools.product(range(5),repeat=2):
        field=s.zeros(5,1);field[a]=x[b]
        cols.append(vector(sparse(X*field.T-field*X.T,x)))
    for a,b in itertools.product(range(5),repeat=2):
        E=s.zeros(5);E[a,b]=1
        cols.append(vector(sparse(E*P+P*E.T-x[b]*P.diff(x[a]),x)))
    cols.append(vector(ps))
    cols += [vector(sparse(matrix(enc,x),x)) for enc in r['all11_actual_parameter_directions']]
    F=s.Matrix.hstack(*cols)
    assert D*F==s.zeros(350,F.cols)
    dm=DomainMatrix.from_Matrix(D).convert_to(s.QQ)
    rr,pivots=dm.rref();rr=rr.to_Matrix();assert len(pivots)==98
    free=[i for i in range(150)if i not in pivots];kernel=[]
    for j in free:
        col=s.zeros(150,1);col[j]=1
        for i,p in enumerate(pivots):col[p]=-rr[i,j]
        denominator=s.ilcm(1,*[s.denom(v)for v in col]);kernel.append(col*denominator)
    K=s.Matrix.hstack(*kernel);assert D*K==s.zeros(350,52)
    allpivots=list(DomainMatrix.from_Matrix(F.row_join(K)).convert_to(s.QQ).rref()[1])
    fpivots=[j for j in allpivots if j<F.cols]
    normal_indices=[j-F.cols for j in allpivots if j>=F.cols]
    assert len(fpivots)==49 and len(normal_indices)==3
    R=[];rs=[]
    for j in normal_indices:
        data={}
        for (ij,e),v in zip(COLUMNS,K[:,j]):
            if v:data.setdefault(ij,{})[e]=int(v)
        rs.append(data)
        value=s.zeros(5)
        for ij,poly in data.items():
            value[ij]=sum(c*s.prod(xx**a for xx,a in zip(x,e))for e,c in poly.items())
            value[ij[::-1]]=-value[ij]
        R.append(value)
    q=[];normal_pairs=list(itertools.combinations_with_replacement(range(3),2))
    squares=[jacobi(data)for data in rs]
    for i,j in normal_pairs:
        if i==j:q.append(squares[i])
        else:q.append(jacobi(sparse(R[i]+R[j],x))-squares[i]-squares[j])
    Q=s.Matrix.hstack(*q)
    dind=D[:,list(pivots)]
    rows=list(DomainMatrix.from_Matrix(dind.T).convert_to(s.QQ).rref()[1])
    block=DomainMatrix.from_Matrix(D.extract(rows,list(pivots))).convert_to(s.QQ)
    rhs=DomainMatrix.from_Matrix(Q.extract(rows,range(6))).convert_to(s.QQ)
    correction=block.inv().matmul(rhs).to_Matrix()
    obstruction=Q-dind*correction
    assert obstruction.extract(rows,range(6))==s.zeros(98,6)
    reduced,qpivots=DomainMatrix.from_Matrix(obstruction).convert_to(s.QQ).rref()
    reduced=reduced.to_Matrix();qrank=len(qpivots)
    n=s.symbols('n0:3');mons=[n[i]*n[j]for i,j in normal_pairs]
    polynomials=[s.expand(sum(reduced[i,j]*mons[j]for j in range(6)))for i in range(qrank)]
    gb=s.groebner(polynomials,*n,order='grevlex') if polynomials else None
    zero_dimensional=bool(gb and gb.is_zero_dimensional)
    report=dict(task_id='E_quadric3',all_passed=True,
        actual_source_tangent_rank=49,full_D_rank=98,full_kernel_dimension=52,
        normal_count=3,normal_kernel_indices=normal_indices,
        source_independent_column_indices=fpivots,
        normals=[{str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]!=0}for A in R],
        D_times_all_normals_zero=True,quadratic_obstruction_rank=qrank,
        normal_monomial_order=list(map(str,mons)),
        quadratic_initial_generators=list(map(str,polynomials)),
        grevlex_basis=list(map(str,gb.polys))if gb else [],
        quadratic_support_zero_dimensional=zero_dimensional,
        local_projective_dimension_upper_if_zero_dimensional=24 if zero_dimensional else None,
        source_componenthood_proved=zero_dimensional,
        full_initial_ideal_and_nil_thickness_not_claimed=True,
        eliminated_D_rows=rows,eliminated_D_columns=list(pivots),
        coefficient_matrix_sha256=hashlib.sha256(json.dumps(obstruction.tolist(),default=str).encode()).hexdigest(),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in
            [SOURCE,'code/classification_rank3_twzero_normal_quadrics.py','code/classification_rank3_mixed_component_review.py']},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(quadratic_rank=qrank,zero_dimensional=zero_dimensional,
                         quadrics=list(map(str,polynomials)),seconds=report['seconds'])))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)
