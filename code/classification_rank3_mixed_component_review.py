#!/usr/bin/env python3
"""Independent full-ambient certificate for the mixed rank-three orbit.

No producer algebra, monomial enumerator, Jacobi, orbit, or rank helper is
imported. Coefficients are assembled by sparse integer monomial arithmetic.
"""
import argparse
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
TRIPLES=list(itertools.combinations(range(5),3))
QUADS=list(itertools.combinations(range(5),4))


def monomials(n):
    # Reverse lexicographic exponent lists, obtained without the producer's
    # combinations-with-replacement monomial generator.
    return sorted((e for e in itertools.product(range(n+1),repeat=5) if sum(e)==n),reverse=True)


M2=monomials(2);M4=monomials(4)
COLUMNS=list(itertools.product(PAIRS,M2));ROWS=list(itertools.product(QUADS,M4))


def entry(A,i,j):
    return A.get((i,j),{}) if i<j else {e:-v for e,v in A.get((j,i),{}).items()} if i>j else {}


def jacobi(A,raw=False):
    cubic={}
    for i,j,k in TRIPLES:
        f=defaultdict(int)
        for a,b,c in ((i,j,k),(j,k,i),(k,i,j)):
            for ell in range(5):
                for e,v in entry(A,a,ell).items():
                    for ee,vv in entry(A,b,c).items():
                        if ee[ell]:
                            exponent=tuple(e[q]+ee[q]-(q==ell) for q in range(5))
                            f[exponent]+=v*vv*ee[ell]
        cubic[(i,j,k)]={e:v for e,v in f.items() if v}
    if raw:return cubic
    fourth={}
    for inds in QUADS:
        f=defaultdict(int)
        for at,ell in enumerate(inds):
            for e,v in cubic[inds[:at]+inds[at+1:]].items():
                f[tuple(e[q]+(q==ell) for q in range(5))]+=(-1)**at*v
        fourth[inds]=f
    return s.Matrix([fourth[inds].get(e,0) for inds,e in ROWS])


def sparse(A,x):
    return {ij:{e:int(v) for e,v in s.Poly(A[ij],*x).terms() if v}
            for ij in PAIRS if A[ij]!=0}


def vector(A):
    return s.Matrix([A.get(ij,{}).get(e,0) for ij,e in COLUMNS])


def minor_certificate(M,rank):
    dm=DomainMatrix.from_Matrix(M).convert_to(s.QQ)
    cols=list(dm.rref()[1]);assert len(cols)==rank
    rows=list(dm.extract(range(M.rows),cols).transpose().rref()[1])
    det=dm.extract(rows,cols).det();assert det
    return dict(shape=list(M.shape),rank_QQ=rank,rows=rows,columns=cols,determinant=str(det))


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();h,u,v,z,w=x=s.symbols('h u v z w')
    P=s.zeros(5)
    for ij,f in zip(itertools.combinations(range(1,5),2),
                    [-h*v,3*h*z-2*u*w,2*h*w+2*u*u,
                     2*h*w+2*u*u+2*w*v,-2*u*v,4*w*w+6*u*z]):
        P[ij]=f;P[ij[::-1]]=-f
    R=s.zeros(5)
    for ij,f in {(0,1):-2*h*u+2*u*v,(0,2):-4*h*v+6*v*v,
                 (0,3):4*h*z-2*v*z,(0,4):2*u*u,
                 (1,2):-u*v,(1,3):3*u*z,(2,3):11*v*z}.items():
        R[ij]=f;R[ij[::-1]]=-f
    ps,rs=sparse(P,x),sparse(R,x)
    assert all(not v for v in jacobi(ps,raw=True).values())
    assert jacobi(ps)==s.zeros(350,1)
    columns=[]
    for ij,e in COLUMNS:
        variation={key:dict(value) for key,value in ps.items()}
        variation.setdefault(ij,{})[e]=variation.get(ij,{}).get(e,0)+1
        # A single coordinate-monomial bivector has zero self-Jacobi.
        columns.append(jacobi(variation))
    D=s.Matrix.hstack(*columns)
    # Actual Euler-vertical translations and actual finite GL5 action derivatives.
    V=[];O=[];X=s.Matrix(x)
    for a,b in itertools.product(range(5),repeat=2):
        field=s.zeros(5,1);field[a]=x[b]
        vertical=X*field.T-field*X.T
        V.append(vector(sparse(vertical,x)))
        E=s.zeros(5);E[a,b]=1
        action=E*P+P*E.T-x[b]*P.diff(x[a])
        O.append(vector(sparse(action,x)))
    VM=s.Matrix.hstack(*V);M=VM.row_join(s.Matrix.hstack(*O)).row_join(vector(ps))
    N=M.row_join(vector(rs))
    assert D*N==s.zeros(350,N.cols)
    certificates={name:minor_certificate(A,rank) for name,A,rank in
                  [('full350_D',D,102),('vertical',VM,24),('orbit_vertical_scale',M,47),
                   ('orbit_vertical_scale_plus_R',N,48)]}
    left_data={7:-18,14:18,18:-9,20:9,46:-5,73:36,78:18,80:-36,91:-9,117:4}
    left=s.zeros(1,350)
    for i,c in left_data.items():left[i]=c
    assert left*D==s.zeros(1,150)
    obstruction=(left*jacobi(rs))[0];assert obstruction==-140
    pf=s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
    Q=h*w+u*u;assert s.expand(pf-4*Q*Q)==0
    qmatrix=s.hessian(Q,x)/2;assert qmatrix.rank()==3
    frozen=ROOT/'work/rank3-mixed-support-obstruction.json'
    source=json.loads(frozen.read_text())
    assert source['left_cokernel_witness']=={str(i):str(c) for i,c in left_data.items()}
    report=dict(task_id='E_quadric3',all_passed=True,
        independent_polynomial_and_action_implementation=True,
        raw_Jacobi_all_ten_zero=True,full350_projective_Jacobi_zero=True,
        all_actual48_tangent_directions_killed_by_D=True,certificates=certificates,
        normal_R={str(ij):str(R[ij]) for ij in PAIRS if R[ij]!=0},
        left_cokernel=[dict(row=i,quadruple=ROWS[i][0],exponents=ROWS[i][1],coefficient=c)
                       for i,c in left_data.items()],
        left_D_zero=True,left_quadratic_obstruction=str(obstruction),
        true_formal_one_normal_support_dimension=0,
        actual_projective_orbit_dimension=22,full_projective_local_dimension=22,
        projective_orbit_closure_is_full_component=True,
        divisor=str(4*h*Q**2),quadric_rank=3,
        exact_nil_thickness_not_claimed=True,global_exhaustion_not_claimed=True,
        input_sha256={str(frozen.relative_to(ROOT)):hashlib.sha256(frozen.read_bytes()).hexdigest(),
                      'code/classification_rank3_mixed_component_review.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        matrix_sha256={name:hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest()
                       for name,A in [('D',D),('orbit_vertical_scale',M),('plus_R',N)]},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,component_dimension=22,obstruction=str(obstruction),seconds=report['seconds'])))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)
