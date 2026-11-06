#!/usr/bin/env python3
"""Check one complete quadratic curve at the nonzero-modulus nil source."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,ROWS,jacobi,sparse
ROOT=Path(__file__).resolve().parents[1]

def main(output,fixed_hyperplanes=False):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    start=time.monotonic();paths=['reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_NONZERO_MODULUS_20261107.json','reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_NONZERO_MODULUS_NORMAL_20261107.json']
    source,normal=[json.loads((ROOT/p).read_text())for p in paths];x=s.symbols('h x1:5');t=s.Symbol('t')
    def decode(data):
        A=s.zeros(5)
        for key,f in data.items():
            i,j=map(int,key.split(','));A[i,j]=s.sympify(f);A[j,i]=-A[i,j]
        return A
    def encode(A):return {f'{i},{j}':str(s.expand(A[i,j]))for i,j in itertools.combinations(range(5),2)if A[i,j]}
    P=decode(source['branches'][0]['witness_brackets']);R=decode(normal['normals'][0]);ps=sparse(P,x)
    if fixed_hyperplanes:
        A=s.zeros(5);A[0,4]=A[1,4]=1;A[1,3]=-1;Y=A*s.Matrix(x)
        S=R+A*P+P*A.T-sum((Y[i]*P.diff(x[i])for i in range(5)),s.zeros(5))
        Z=s.Matrix([s.cancel(S[0,j]/x[0])for j in range(5)])
        assert all(s.denom(z)==1 for z in Z)
        E=s.Matrix(x);S=(S-E*Z.T+Z*E.T).applyfunc(s.expand)
        assert all(S[0,j]==0 for j in range(5))
        assert all(s.rem(S[1,j],x[1],x[1])==0 for j in range(2,5))
        Q=P+t*S
        full={str(ijk):s.factor(sum(Q[i,l]*s.diff(Q[j,k],x[l])+Q[j,l]*s.diff(Q[k,i],x[l])+Q[k,l]*s.diff(Q[i,j],x[l])for l in range(5)))for ijk in itertools.combinations(range(5),3)for i,j,k in [ijk]}
        V=s.Matrix([s.cancel(Q[1,j]/x[1])for j in range(2,5)]);B=V.jacobian(x[2:5]);L=s.Symbol('L')
        Dq=x[0]*(Q[1,2]*Q[3,4]-Q[1,3]*Q[2,4]+Q[1,4]*Q[2,3])
        result=dict(task_id='E_reduced_nonnormal',geometric_hyperplane_displacements=['x4','x4-x3'],
            infinitesimal_PGL_matrix=[[str(z)for z in row]for row in A.tolist()],Euler_gauge=list(map(str,Z)),
            fixed_hyperplanes_normal=encode(S),family=encode(Q),
            all_parameter_strict_Jacobi_zero=all(f==0 for f in full.values()),
            nonzero_strict_Jacobi={k:str(f)for k,f in full.items()if f},
            cofactor=list(map(str,V)),cofactor_B=[[str(z)for z in row]for row in B.tolist()],
            cofactor_B_characteristic_polynomial=str(s.factor(B.charpoly(L).as_expr())),
            actual_Euler_Pfaffian=str(s.factor(Dq)),producer_source=Path(__file__).read_text(),
            input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
        with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
        print(json.dumps({k:result[k]for k in ['all_parameter_strict_Jacobi_zero','nonzero_strict_Jacobi','cofactor_B_characteristic_polynomial','seconds']}));return
    columns=[]
    for ij,e in COLUMNS:
        variation={key:dict(value)for key,value in ps.items()};variation.setdefault(ij,{})[e]=variation.get(ij,{}).get(e,0)+1
        columns.append(jacobi(variation))
    D=s.Matrix.hstack(*columns);dm=DomainMatrix.from_Matrix(D).convert_to(s.QQ)
    pivots=list(dm.rref()[1]);assert len(pivots)==102
    rows=list(DomainMatrix.from_Matrix(D[:,pivots].T).convert_to(s.QQ).rref()[1])
    square=jacobi(sparse(R,x));correction=-DomainMatrix.from_Matrix(D.extract(rows,pivots)).convert_to(s.QQ).inv().matmul(DomainMatrix.from_Matrix(square.extract(rows,[0])).convert_to(s.QQ)).to_Matrix()
    assert D[:,pivots]*correction+square==s.zeros(350,1)
    T=s.zeros(5)
    for index,value in zip(pivots,correction):
        ij,e=COLUMNS[index];T[ij]+=value*s.prod(z**k for z,k in zip(x,e))
    for i,j in itertools.combinations(range(5),2):T[j,i]=-T[i,j]
    Q=P+t*R+t*t*T
    J={ijk:s.expand(sum(Q[i,l]*s.diff(Q[j,k],x[l])+Q[j,l]*s.diff(Q[k,i],x[l])+Q[k,l]*s.diff(Q[i,j],x[l])for l in range(5)))for ijk in itertools.combinations(range(5),3)for i,j,k in [ijk]}
    full={str(inds):s.factor(sum((-1)**k*x[inds[k]]*J[inds[:k]+inds[k+1:]]for k in range(4)))for inds in itertools.combinations(range(5),4)}
    cubic=s.zeros(350,1)
    for row,(inds,e) in enumerate(ROWS):
        cubic[row]=s.Poly(s.expand(full[str(inds)]).coeff(t,3),*x).coeff_monomial(s.prod(z**k for z,k in zip(x,e)))
    inv=DomainMatrix.from_Matrix(D.extract(rows,pivots)).convert_to(s.QQ).inv().to_Matrix()
    remainder=cubic-D[:,pivots]*inv*cubic.extract(rows,[0])
    obstruction=None
    if any(remainder):
        index=next(i for i,a in enumerate(remainder)if a)
        functional=s.zeros(1,350);functional[index]=1
        coefficients=D[index,pivots]*inv
        for k,row in enumerate(rows):functional[row]-=coefficients[k]
        assert functional*D==s.zeros(1,150)
        assert (functional*cubic)[0]==remainder[index]
        obstruction=dict(row=index,nonzero_coefficients={str(i):str(a)for i,a in enumerate(functional)if a},
                         pairing=str(remainder[index]))
    Dq=0
    for i in range(5):
        a,b,c,d=[j for j in range(5)if j!=i];Dq+=(-1)**i*x[i]*(Q[a,b]*Q[c,d]-Q[a,c]*Q[b,d]+Q[a,d]*Q[b,c])
    result=dict(task_id='E_reduced_nonnormal',single_quadratic_candidate=True,
        exact_second_order_cancellation=True,all_parameter_projective_Jacobi_zero=all(f==0 for f in full.values()),
        family=encode(Q),second_order_correction=encode(T),nonzero_projective_Jacobi={k:str(f)for k,f in full.items()if f},
        actual_Euler_Pfaffian=str(s.factor(Dq/5)),
        genuine_103_slice_regular_columns=pivots,genuine_103_slice_equation_rows=rows,
        projected_cubic_nonzero=obstruction is not None,cubic_cokernel_functional=obstruction,
        no_third_or_higher_order_correction_computed=True,
        independent_membership_not_claimed=True,producer_source=Path(__file__).read_text(),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k]for k in ['all_parameter_projective_Jacobi_zero','nonzero_projective_Jacobi','seconds']}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--fixed-hyperplanes',action='store_true');args=ap.parse_args();main(args.output,args.fixed_hyperplanes)
