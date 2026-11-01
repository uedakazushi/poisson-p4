#!/usr/bin/env python3
"""Exact product-pair coordinates for the seven projected deeper-binary quadrics.

This is an algebraic reduction on the previously treated factor complement;
it proves no new removal or full-fiber exhaustion statement.
"""
import argparse,hashlib,json,time
from itertools import combinations_with_replacement
from pathlib import Path
import numpy as np
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'code/certificates/s61_binary_deeper_cubic_20261101.json'
SYZ=ROOT/'research/archive/binary_deeper_pid_20261101/extra_syzygy.json'
QUAD=ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
NOTE=ROOT/'research/notes/S61_BINARY_DEEPER_PRODUCT_PAIR_20261102.md'


def verify():
    started=time.time();certificate=json.loads(CERT.read_text())
    F=s.QQ.frac_field('e','d');R,*v=ring('v0:7',F)
    keep=[0,1,2,3,4,6,8];old=list(combinations_with_replacement(range(9),2));Q=[];indices=[]
    for i,row in enumerate(certificate['quadratic_rref']['rref']):
        if sum(c!='0' for c in row)==1:continue
        q=R.zero
        for (a,b),c in zip(old,row):
            if c!='0':
                assert a in keep and b in keep
                q+=F.from_sympy(s.sympify(c))*v[keep.index(a)]*v[keep.index(b)]
        Q.append(q);indices.append(i)
    assert indices==[2,4,7,12,14,18,20]
    source=json.loads(SYZ.read_text());M=[[F.zero]*7 for _ in range(7)]
    for (i,j),c in zip(source['sources'],source['coefficients']):M[i][j]=F.from_sympy(s.sympify(c))
    assert sum((sum((R(c)*v[j] for j,c in enumerate(row)),R.zero)*q for row,q in zip(M,Q)),R.zero)==0
    L=3*v[0]-v[1]+v[5];N=(Q[6]-Q[5]).exquo(L)
    assert max(map(sum,N))==1 and N.get((0,0,0,1,0,0,0))==-F.one
    scalar=M[6][3];aa=[row[5]/scalar for row in M];bb=[row[3]/scalar for row in M]
    assert bb==[F.zero]*5+[-F.one,F.one] and aa[1]==F.one
    assert all(M[i][j]/scalar==aa[i]*[3,-1,0,0,0,1,0][j]+bb[i]*(j==3)
               for i in range(7) for j in range(7))
    A=sum((a*q for a,q in zip(aa,Q)),R.zero)
    assert A==-v[3]*N and Q[6]-Q[5]==L*N
    selected=[0,2,3,4,5]
    row_change=[[F.one if i==j else F.zero for j in range(7)] for i in selected]+[bb,[-a for a in aa]]
    row_det=DomainMatrix(row_change,(7,7),F).det();assert row_det in [F.one,-F.one]

    S,*x=ring('x0:7',F)
    inverse={0:x[0],1:3*x[0]+x[3]-x[5],2:x[1],4:x[2],5:x[3],6:x[4]}
    H=sum((c*inverse[m.index(1)] for m,c in (N+v[3]).items()),S.zero)
    inverse[3]=H-x[6]
    coordinate_matrix=[[inverse[i].get(tuple(int(j==k) for j in range(7)),F.zero)
                        for k in range(7)] for i in range(7)]
    coordinate_det=DomainMatrix(coordinate_matrix,(7,7),F).det()
    assert coordinate_det in [F.one,-F.one]
    def transform(q):
        result=S.zero
        for mon,c in q.items():
            term=S(c)
            for j,power in enumerate(mon):
                if power:term*=inverse[j]**power
            result+=term
        return result
    assert transform(L)==x[5] and transform(N)==x[6]
    NL=x[5]*x[6];N2=x[6]**2-x[6]*H
    new=[transform(Q[i]) for i in selected]
    square=(0,0,0,0,0,0,2);mixed=(0,0,0,0,0,1,1)
    for i in range(5):
        new[i]-=new[i].get(square,F.zero)*N2
        new[i]-=new[i].get(mixed,F.zero)*NL
    N2-=N2.get(mixed,F.zero)*NL
    new.extend([NL,N2]);H_reduced=H-H.get((0,0,0,0,0,1,0),F.zero)*x[5]
    assert N2==x[6]**2-x[6]*H_reduced
    assert all(m[5]==m[6]==0 for m in H_reduced)
    assert all(q.get(square,F.zero)==q.get(mixed,F.zero)==0 for q in new[:5])
    known=[F.field.ring.from_expr(s.sympify(z)) for z in json.loads(QUAD.read_text())['irreducible_minor_factors']]
    coefficients=aa+bb+[c for q in list(inverse.values())+new for c in q.values()]
    for coefficient in coefficients:
        denominator=coefficient.denom
        for factor in known:
            while denominator.rem(factor)==0:denominator=denominator.exquo(factor)
        assert max(map(sum,denominator))==0

    # A single nonzero minor supplies the generic lower bound; the displayed
    # exact nonzero syzygy supplies the upper bound over QQ(e,d).
    prime=32003;point={s.Symbol('e'):2,s.Symbol('d'):3}
    def evaluate(c):
        c=s.Rational(F.to_sympy(c).subs(point))
        assert int(c.q)%prime
        return int(c.p)%prime*pow(int(c.q),-1,prime)%prime
    mons=list(combinations_with_replacement(range(7),3))
    lookup={tuple(mon.count(i) for i in range(7)):j for j,mon in enumerate(mons)}
    matrix=np.zeros((49,84),dtype=np.int64)
    for i,q in enumerate(Q):
        for j in range(7):
            for mon,c in (q*v[j]).items():matrix[7*i+j,lookup[mon]]=evaluate(c)
    rank,rows=elimination_rank(matrix.T,prime);assert rank==48
    _,cols=elimination_rank(matrix[rows,:],prime)
    determinant=determinant_mod(matrix[rows,:][:,cols],prime);assert determinant
    paths=[Path(__file__),CERT,SYZ,QUAD,NOTE,ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py']
    return dict(id='S61_BINARY_DEEPER_PRODUCT_PAIR_20261102',all_passed=True,
        coefficient_field='QQ(e,d)',old_light_indices=keep,old_quadratic_rows=indices,
        L=str(L.as_expr()),N=str(N.as_expr()),normalized_syzygy_L_coefficients=[str(F.to_sympy(a)) for a in aa],
        normalized_syzygy_v3_coefficients=[str(F.to_sympy(a)) for a in bb],
        quadratic_basis_determinant=str(F.to_sympy(row_det)),coordinate_determinant=str(F.to_sympy(coordinate_det)),
        new_to_old_linear=[str(inverse[i].as_expr()) for i in range(7)],
        new_quadrics=[str(q.as_expr()) for q in new],H=str(H_reduced.as_expr()),
        denominator_coefficients_checked=len(coefficients),all_denominators_supported_on_13_previously_treated_factors=True,
        cubic_multiplication=dict(source_dimension=49,target_dimension=84,generic_rank=48,generic_kernel_dimension=1,
            exact_rank_two_linear_syzygy_checked=True,prime=prime,point=[2,3],minor_rows=rows,minor_columns=cols,minor_mod_prime=determinant),
        scope='Exact change of coordinates and quadratic generators on the complement of the thirteen treated core factors.',
        limitations=['No new core-curve removal','No fourth-degree full-rank certificate on the entire core','No ambient nilpotence or full-fiber exhaustion'],
        proof_note=str(NOTE.relative_to(ROOT)),input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))
