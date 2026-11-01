#!/usr/bin/env python3
"""Independent parabolic cyclicity, shear, and Bott checks for AB degree two."""
import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1]
    report_path=root/'reports/AB_ORIGIN_DEGREE_TWO_20261101.json'
    main=json.loads(report_path.read_text());assert main['all_passed']
    # Reconstruct the seven-dimensional parabolic in q=(x,z), p=(y,w) order.
    exponents=[(i,j,k,6-i-j-k) for i in range(7) for j in range(7-i) for k in range(7-i-j)]
    index={m:i for i,m in enumerate(exponents)};assert len(index)==84
    generators=[]
    for i in range(2):
        for j in range(2):
            A=s.zeros(4);A[i,j]=1;A[2+j,2+i]=-1;generators.append(A)
    for i,j in combinations_with_replacement(range(2),2):
        A=s.zeros(4);A[2+i,j]=1;A[2+j,i]=1;generators.append(A)
    def act(A,vector):
        out={}
        for m,c in vector.items():
            for i in range(4):
                if not m[i]:continue
                for j in range(4):
                    if A[i,j]:
                        n=list(m);n[i]-=1;n[j]+=1;n=tuple(n)
                        out[n]=out.get(n,s.S.Zero)+c*m[i]*A[i,j]
        return {m:c for m,c in out.items() if c}
    pivots={};queue=[]
    def insert(row):
        row=dict(row)
        while row:
            mon=min(row)
            if mon not in pivots:
                c=row[mon];row={m:v/c for m,v in row.items()};pivots[mon]=row;queue.append(row);return
            c=row[mon]
            for m,v in pivots[mon].items():
                row[m]=row.get(m,s.S.Zero)-c*v
                if not row[m]:del row[m]
    insert({(0,0,6,0):s.S.One})
    cursor=0
    while cursor<len(queue):
        vector=queue[cursor];cursor+=1
        for A in generators:insert(act(A,vector))
    assert len(pivots)==84
    # Independent tensor-block shear calculation, without the main code's
    # matrix-to-core inverse or pushforward helper.
    x,z,u=s.symbols('x z u');h=s.symbols('h0:5');a=s.symbols('a0:3');g=s.symbols('g0:3')
    H=sum(h[i]*x**(4-i)*z**i for i in range(5));S=a[0]*x*x+a[1]*x*z+a[2]*z*z
    F=s.Matrix([x*S/4+s.diff(H,z),z*S/4-s.diff(H,x)])
    A=F.jacobian([x,z]);G=s.Matrix([[0,g[0]*x*x+g[1]*x*z+g[2]*z*z],[-g[0]*x*x-g[1]*x*z-g[2]*z*z,0]])
    C=s.diag(0,u);block=s.zeros(4);block[:2,2:]=A;block[2:,:2]=-A.T;block[2:,2:]=G
    shear=s.eye(4);shear[2:,:2]=C;transformed=shear*block*shear.T
    assert transformed[:2,2:]==A
    shift=s.Poly(s.expand(transformed[2,3]-G[0,1]),x,z).coeff_monomial(x*x)
    assert shift==12*u*h[0]
    core=[(-2,0),(-1,-1),(0,-2),(-3,1),(-2,0),(-1,-1),(0,-2),(1,-3)]
    vertical=[(-3,-1),(-2,-2),(-1,-3)]
    cross=Counter(tuple(a[i]+b[i] for i in range(2)) for a in core for b in vertical)
    assert cross[(-6,0)]==1
    assert all(sum(v)==-2 for v in core) and all(sum(v)==-4 for v in vertical)
    # Apply the primary theorem to Sigma^(a,b) U* using lambda=(-b,-a).
    table=[];piece_dims=[[0,0] for _ in range(3)]
    for piece,items in enumerate(main['Levi_degree_two_decomposition']):
        for a,b,multiplicity in items:
            weight=(a+2,b+1)
            if 0 in weight or abs(weight[0])==abs(weight[1]):
                table.append([piece,a,b,multiplicity,None]);continue
            degree=sum(v<0 for v in weight)+(weight[0]+weight[1]<0)
            ordered=sorted(map(abs,weight),reverse=True);p,q=ordered[0]-2,ordered[1]-1
            dim=(p-q+1)*(q+1)*(p+2)*(p+q+3)//6
            assert degree in (0,1);piece_dims[piece][degree]+=multiplicity*dim
            table.append([piece,a,b,multiplicity,degree,p,q,dim])
    assert piece_dims==[[84,84],[486,0],[275,0]]
    assert main['triangular_pivot_product_mod_prime']==11854
    assert len(main['product_indices'])==len(main['pivot_monomials'])==761
    paths=[Path(__file__),report_path,root/'code/ab_origin_degree_two_20261101.py',root/'research/notes/AB_ORIGIN_DEGREE_TWO_20261101.md',root/'research/notes/AB_INCIDENCE_20261101.md']
    return dict(id='AB_ORIGIN_DEGREE_TWO_REVIEW_20261101',all_passed=True,
        independent_parabolic_cyclic_span_rank=84,parabolic_generator_count=7,
        independent_shear_coefficient=str(shift),unique_candidate_weight=[-6,0],
        cohomology_piece_dimensions=piece_dims,independent_Bott_table=table,
        accepted_dimension=761,main_script_rerun_verified=True,
        primary_source='https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9880&wshow=paper',
        primary_source_checked='Fonarev, Theorem 3.2, complete statement in primary HTML',
        proof_note='research/notes/AB_ORIGIN_DEGREE_TWO_REVIEW_20261101.md',
        limitations=['Degree two only; no all-degree generation or normality at the origin.',
                    'The main modular coefficient minor was rerun, not independently regenerated here.'],
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','independent_parabolic_cyclic_span_rank','accepted_dimension','seconds']}))


if __name__=='__main__':main()
