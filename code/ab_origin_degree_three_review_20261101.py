#!/usr/bin/env python3
"""Independent low-degree character and parabolic checks for AB degree three."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from collections import Counter
from itertools import combinations, combinations_with_replacement
import sympy as s


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1]
    main_path=root/'reports/AB_ORIGIN_DEGREE_THREE_20261101.json';main=json.loads(main_path.read_text());assert main['all_passed']
    generators=[]
    for i in range(2):
        for j in range(2):
            A=s.zeros(4);A[i,j]=1;A[2+j,2+i]=-1;generators.append(A)
    for i,j in combinations_with_replacement(range(2),2):
        A=s.zeros(4);A[2+i,j]=1;A[2+j,i]=1;generators.append(A)
    def action(A,vector):
        result={}
        def add(mon,ij,value):
            key=(mon,ij);result[key]=result.get(key,s.S.Zero)+value
            if not result[key]:del result[key]
        for (mon,ij),c in vector.items():
            for i in range(4):
                if not mon[i]:continue
                for j in range(4):
                    if A[i,j]:
                        n=list(mon);n[i]-=1;n[j]+=1;add(tuple(n),ij,c*mon[i]*A[i,j])
            if ij is not None:
                i,j=ij
                for k in range(4):
                    if k!=j and A[i,k]:add(mon,tuple(sorted((k,j))),c*A[i,k]*(1 if k<j else -1))
                    if k!=i and A[j,k]:add(mon,tuple(sorted((i,k))),c*A[j,k]*(1 if i<k else -1))
        return result
    def cyclic(wedge):
        cache={};queue=[]
        def insert(row):
            row=dict(row)
            while row:
                pivot=min(row)
                if pivot not in cache:
                    c=row[pivot];row={k:v/c for k,v in row.items()};cache[pivot]=row;queue.append(row);return
                c=row[pivot]
                for k,v in cache[pivot].items():
                    row[k]=row.get(k,s.S.Zero)-c*v
                    if not row[k]:del row[k]
        insert({((0,0,8,0),wedge):s.S.One});cursor=0
        while cursor<len(queue):
            row=queue[cursor];cursor+=1
            for A in generators:insert(action(A,row))
        return len(cache)
    assert cyclic(None)==165
    cartan_cyclic_rank=cyclic((2,3));assert cartan_cyclic_rank==429
    # Compute the core characters by relation spans, not a Groebner basis.
    a=s.symbols('s0:3');h=s.symbols('h0:5');variables=(*a,*h)
    N=s.Matrix([[h[2]/6+a[1]/48,-h[1]/4-a[0]/48,h[0]],
                [h[3]/2+a[2]/24,-h[2]/3,h[1]/2-a[0]/24],
                [h[4],-h[3]/4+a[2]/48,h[2]/6-a[1]/48]])
    relations=[s.expand(N.extract(i,j).det()) for i in combinations(range(3),2) for j in combinations(range(3),2)]
    weights=[(2,0),(1,1),(0,2),(3,-1),(2,0),(1,1),(0,2),(-1,3)]
    chars=[];relation_ranks=[]
    for degree in range(4):
        inds=list(combinations_with_replacement(range(8),degree))
        monomials=[s.prod(variables[i] for i in ii) for ii in inds]
        wts=[tuple(sum(weights[i][j] for i in ii) for j in range(2)) for ii in inds]
        polys=[]
        if degree>=2:
            for f in relations:
                for mul in combinations_with_replacement(range(8),degree-2):polys.append(s.Poly(f*s.prod(variables[i] for i in mul),*variables))
        character=Counter(wts);total=0
        for wt in list(character):
            select=[i for i,w in enumerate(wts) if w==wt]
            matrix=s.Matrix([[f.coeff_monomial(monomials[i]) for i in select] for f in polys])
            rank=len(matrix.to_DM().rref()[1]) if polys else 0
            character[wt]-=rank;total+=rank
        assert sum(character.values())==(degree+1)**3
        chars.append(character);relation_ranks.append(total)
    vertical=[(3,1),(2,2),(1,3)];hco=[Counter(),Counter()];dims=[0,0]
    for vd in range(4):
        cv=Counter(tuple(sum(vertical[i][k] for i in mon) for k in range(2)) for mon in combinations_with_replacement(range(3),vd))
        prod=Counter()
        for aa,n in chars[3-vd].items():
            for bb,m in cv.items():prod[(aa[0]+bb[0],aa[1]+bb[1])]+=n*m
        prod=Counter({k:v for k,v in prod.items() if v})
        while prod:
            aa,bb=max(prod);mult=prod[aa,bb]
            for k in range(aa-bb+1):
                wt=(aa-k,bb+k);prod[wt]-=mult;assert prod[wt]>=0
                if not prod[wt]:del prod[wt]
            u,v=aa+2,bb+1
            if not u or not v or abs(u)==abs(v):continue
            deg=int(u<0)+int(v<0)+int(u+v<0);assert deg in (0,1)
            x,y=sorted((abs(u),abs(v)),reverse=True);x-=2;y-=1
            dim=(x-y+1)*(y+1)*(x+2)*(x+y+3)//6
            hco[deg][x,y]+=mult;dims[deg]+=mult*dim
    assert dims==[7285,594] and hco[1]==Counter({(9,1):1,(8,0):1})
    assert hco[0][9,1]==1 and hco[0][8,0]==2
    # Independently apply the two elementary shear substitutions to the only
    # possible Sym8 image, and test the resulting relations on the core chart.
    g0,g1,g2,u=s.symbols('g0 g1 g2 u');ca,cb,cc=s.symbols('A B C')
    f=ca*h[0]**2*g1+cb*h[0]*h[1]*g0+cc*h[0]*a[0]*g0
    levi=s.diff(f.subs({h[1]:h[1]-4*u*h[0],g1:g1-2*u*g0},simultaneous=True),u).subs(u,0).expand()
    radical=s.diff(f.subs({g0:g0+12*u*h[0],g1:g1+u*(6*h[1]-a[0]/2)},simultaneous=True),u).subs(u,0).expand()
    assert s.expand(levi+(2*ca+4*cb)*h[0]**2*g0)==0
    assert s.expand(radical-(6*ca+12*cb)*h[0]**2*h[1]-(-ca/2+12*cc)*h[0]**2*a[0])==0
    assert s.diff(h[0]*(g0+12*u*h[0])**2,u,2)==288*h[0]**3
    aa,bb,dd,ee=s.symbols('aa bb dd ee');chart=s.Matrix([1,aa,bb])*s.Matrix([[-aa*dd-bb*ee,dd,ee]])
    E=N.reshape(9,1).jacobian(variables);rows=list(E.T.to_DM().rref()[1]);values=E[rows,:].inv()*chart.reshape(9,1)[rows,:]
    sub=dict(zip(variables,values));polyvariables=(aa,bb,dd,ee,g0,g1,g2)
    candidates=[s.Poly(q.subs(sub,simultaneous=True),*polyvariables) for q in [h[0]**2*g1,h[0]*h[1]*g0,h[0]*a[0]*g0]]
    monomials=sorted(set().union(*(set(q.monoms()) for q in candidates)))
    assert s.Matrix([[q.coeff_monomial(mon) for q in candidates] for mon in monomials]).rank()==3
    coeff=[]
    for q in [levi,radical]:coeff.extend(s.Poly(q.subs(sub,simultaneous=True),*polyvariables).coeffs())
    constraints=s.Matrix(coeff).jacobian([ca,cb,cc]);assert constraints.rank()==2
    assert constraints*s.Matrix([1,-s.Rational(1,2),s.Rational(1,24)])==s.zeros(constraints.rows,1)
    paths=[Path(__file__),main_path,root/'code/ab_origin_degree_three_20261101.py',root/'research/notes/AB_ORIGIN_DEGREE_THREE_20261101.md',root/'research/notes/AB_INCIDENCE_20261101.md']
    return dict(id='AB_ORIGIN_DEGREE_THREE_REVIEW_20261101',all_passed=True,
        independent_Sym8_cyclic_rank=165,independent_partition91_cyclic_rank=429,
        independent_core_relation_ranks=relation_ranks,independent_associated_graded_cohomology_dimensions=dims,
        independent_H0_multiplicities={str(k):v for k,v in hco[0].items()},independent_H1_multiplicities={str(k):v for k,v in hco[1].items()},
        fiber_candidate_dimension=3,fiber_annihilator_constraint_rank=2,accepted_degree_three_dimension=6691,
        proof_note='research/notes/AB_ORIGIN_DEGREE_THREE_REVIEW_20261101.md',
        limitations=['The6691 modular coefficient minor was rerun, not independently regenerated here.','No all-degree normality conclusion.'],
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args();r=verify()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','independent_partition91_cyclic_rank','accepted_degree_three_dimension','seconds']}))


if __name__=='__main__':main()
