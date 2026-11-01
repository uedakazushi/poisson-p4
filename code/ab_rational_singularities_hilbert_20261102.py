#!/usr/bin/env python3
"""Finite certificate accompanying the AB rational-singularity/Hilbert proof."""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from ab_origin_degree_four_backend_20261101 import Counter, character, core_weights, vertical_weights, leading, decompose, bott, s

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.time()
    old=json.loads((ROOT/'reports/AB_ORIGIN_NORMALITY_20261102.json').read_text())
    assert old['all_passed']
    fiber=[(0,0,1),(1,2,9),(2,3,16),(3,4,9),(4,6,1)]
    tor=[Counter() for _ in range(39)]
    for i,d,r in fiber:
        for k in range(35):tor[i+k][d+k]+=r*math.comb(34,k)
    assert max(d for row in tor for d in row)==40
    assert max(d-i for i,row in enumerate(tor) for d in row)==2
    shifts={}
    for p in range(39):
        allowed={d for j in range(4) if p+j<len(tor) for d in tor[p+j]}
        if allowed:
            assert max(allowed)<=p+5
            shifts[p]=max(allowed)
    assert shifts[1]==6
    dims=[];tables=[]
    for n in range(9):
        value=0;rows=[]
        for j in range(n+1):
            cy=character(n-j,core_weights,leading)
            cv=character(j,vertical_weights);prod=Counter()
            assert sum(cy.values())==(n-j+1)**3
            for a,ca in cy.items():
                for b,cb in cv.items():prod[a[0]+b[0],a[1]+b[1]]+=ca*cb
            for a,b,m in decompose(prod):
                co=bott(a,b)
                if co:
                    degree,wt,dim=co
                    assert degree in (0,1)
                    value+=(-1)**degree*m*dim
                rows.append(dict(vertical_degree=j,Levi_weight=[a,b],multiplicity=m,cohomology=co))
        dims.append(value);tables.append(rows)
    assert dims[:6]==[1,45,761,6691,38570,166818]
    numerator=[sum((-1)**j*math.comb(10,j)*dims[n-j] for j in range(min(n,10)+1)) for n in range(9)]
    assert numerator==[1,35,356,986,715,91,0,0,0]
    h=numerator[:6]
    for n in range(9):assert dims[n]==sum(h[i]*math.comb(n-i+9,9) for i in range(min(n,5)+1))
    n=s.symbols('n')
    assert s.expand(((n+2)*(n+1)/2)**2-((n+1)*n/2)**2-(n+1)**3)==0
    paths=[str(Path(__file__).relative_to(ROOT)),
           'code/ab_origin_degree_four_backend_20261101.py',
           'reports/AB_ORIGIN_NORMALITY_20261102.json',
           'reports/AB_ORIGIN_NORMALITY_REVIEW_20261102.json',
           'research/notes/AB_INCIDENCE_20261101.md',
           'research/notes/AB_ORIGIN_NORMALITY_20261102.md',
           'research/notes/AB_RATIONAL_SINGULARITIES_HILBERT_20261102.md']
    report=dict(all_passed=True,seconds=time.time()-started,
                all_local_Tor=[dict(row) for row in tor],maximum_internal_shift=40,
                maximum_internal_minus_homological=2,global_resolution_shift_upper=shifts,
                Hilbert_values=dims,Euler_decompositions=tables,dimension=10,ambient_variables=45,
                Hilbert_numerator=h,Hilbert_denominator_power=10,degree=sum(h),
                regularity=5,a_invariant=-5,projective_dimension=35,ideal_generator_degree_upper=6,
                quadratic_relations=274,geometric_claim='AB has rational singularities; see the explicit blowup proof.',
                limitations=['Reduced AB only; no T01 nonreduced fiber assertion.',
                             'Euler operations use the stated frozen backend.',
                             'Geometric resolution/cohomology reasoning is a written proof, not a program-only inference.'],
                input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    assert report['degree']==2184
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print('AB RATIONAL/HILBERT PASS',report['seconds'])

if __name__=='__main__':main()
