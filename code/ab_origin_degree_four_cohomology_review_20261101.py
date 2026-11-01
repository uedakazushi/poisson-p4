#!/usr/bin/env python3
"""Independent rank-one character, Bott, and hand-derived P-map constraints."""
import argparse
import hashlib
import json
import time
from collections import Counter
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    def symmetric(degree,weights):
        c=Counter()
        for mon in cwr(range(len(weights)),degree):c[tuple(sum(weights[i][j]for i in mon)for j in range(2))]+=1
        return c
    def core(degree):
        E=[(2,0),(1,1),(0,2)];out=Counter()
        for n,sign in [(degree,1),(degree-1,-1)]:
            if n<0:continue
            c=symmetric(n,E)
            for a,na in c.items():
                for b,nb in c.items():out[(degree+a[0]-b[0],degree+a[1]-b[1])]+=sign*na*nb
        assert all(n>=0 for n in out.values())
        return Counter({a:n for a,n in out.items()if n})
    def decompose(c):
        c=c.copy();out=[]
        while c:
            a,b=max(c);n=c[a,b];assert a>=b and n>0;out.append((a,b,n))
            for j in range(a-b+1):
                wt=a-j,b+j;c[wt]-=n;assert c[wt]>=0
                if not c[wt]:del c[wt]
        return out
    H=[Counter(),Counter()];totals=[0,0];table=[]
    for vertical in range(5):
        cy=core(4-vertical);cv=symmetric(vertical,[(3,1),(2,2),(1,3)]);assert sum(cy.values())==(5-vertical)**3
        product=Counter()
        for a,na in cy.items():
            for b,nb in cv.items():product[a[0]+b[0],a[1]+b[1]]+=na*nb
        for a,b,n in decompose(product):
            x,y=a+2,b+1
            if not x*y or abs(x)==abs(y):co=None
            else:
                degree=sum(v<0 for v in [x-y,y,x+y,x]);assert degree in [0,1]
                A,B=sorted([abs(x),abs(y)],reverse=True);A-=2;B-=1
                dim=(A-B+1)*(B+1)*(A+2)*(A+B+3)//6
                co=[degree,[A,B],dim];H[degree][A,B]+=n;totals[degree]+=n*dim
            table.append(dict(vertical_degree=vertical,Levi_weight=[a,b],multiplicity=n,cohomology=co))
    assert totals==[41621,3051]
    primary=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json').read_text())
    assert table==primary['Borel_Weil_Bott']
    ss=s.symbols('s0:3');hh=s.symbols('h0:5');gg=s.symbols('g0:3');fiber=ss+hh+gg
    weights=[(-2,0),(-1,-1),(0,-2),(-3,1),(-2,0),(-1,-1),(0,-2),(1,-3),(-3,-1),(-2,-2),(-1,-3)]
    a,b,d,e=s.symbols('a b d e');cv=(a,b,d,e,*gg)
    chart=[-12*(a*e+2*d),24*(-a*d-2*b*e),12*(-a*a*d-a*b*e+2*b*d),
           e,a*e-2*d,-3*a*d,-a*a*d-a*b*e-2*b*d,-a*b*d-b*b*e,*gg]
    sub=dict(zip(fiber,chart))
    levi=[0,-2*ss[0],-ss[1],0,-4*hh[0],-3*hh[1],-2*hh[2],-hh[3],0,-2*gg[0],-gg[1]]
    rad=[0]*8+[12*hh[0],6*hh[1]-ss[0]/2,2*hh[2]-ss[1]/4]
    def deriv(f,co):return s.expand(sum(s.diff(f,var)*value for var,value in zip(fiber,co)))
    out=[]
    for A,B in H[1]:
        forms=[s.prod(fiber[i]for i in mon)for mon in cwr(range(11),4)
               if tuple(sum(weights[i][j]for i in mon)for j in range(2))==(-A,-B)]
        polys=[s.Poly(f.subs(sub,simultaneous=True),*cv)for f in forms]
        mons=sorted(set().union(*(set(p.monoms())for p in polys)))
        mat=s.Matrix([[p.coeff_monomial(m)for p in polys]for m in mons]);piv=list(mat.to_DM().rref()[1]);forms=[forms[i]for i in piv]
        ops=[]
        for f in forms:
            q=f
            for _ in range(B+1):q=deriv(q,rad)
            ops.append([s.Poly(deriv(f,levi).subs(sub,simultaneous=True),*cv),s.Poly(q.subs(sub,simultaneous=True),*cv)])
        rows=[]
        for j in range(2):
            mons=sorted(set().union(*(set(p[j].monoms())for p in ops)))
            rows.extend([[p[j].coeff_monomial(m)for p in ops]for m in mons])
        mat=s.Matrix(rows);rank=mat.to_DM().rank();bound=len(forms)-rank
        expected=next(x for x in primary['P_map_constraints']if x['partition']==[A,B])
        assert [len(forms),rank,bound]==[expected['weight_dimension'],expected['constraint_rank'],expected['Hom_upper']]
        out.append(dict(partition=[A,B],weight_dimension=len(forms),constraint_rank=rank,Hom_upper=bound))

    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest for path,digest in primary['input_sha256'].items())
    paths=list(primary['input_sha256'])+[
        'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json',
        'research/notes/AB_ORIGIN_DEGREE_FOUR_20261101.md',
        'research/notes/AB_ORIGIN_DEGREE_THREE_REVIEW_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/AB_ORIGIN_DEGREE_FOUR_COHOMOLOGY_REVIEW_20261101.md']
    return dict(id='AB_ORIGIN_DEGREE_FOUR_COHOMOLOGY_REVIEW_20261101',all_passed=True,
        independent_rank_one_trace_character=True,independent_Bott_table=table,
        associated_graded_cohomology_dimensions=totals,independent_manual_shear_constraints=out,
        P_cyclicity_and_evaluation_proof_accepted=True,
        accepted_H1_vanishing=True,accepted_normalization_degree_four_dimension=38570,
        evaluation_minors_independently_checked_here=False,
        limitations=['This is an independent upper-bound/cohomology review; coordinate-product lower-bound evaluation is audited separately.',
                     'Degree four alone does not prove all-degree generation or origin normality.'],
        input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('AB FOURTH COHOMOLOGY REVIEW PASS',result['seconds'])
