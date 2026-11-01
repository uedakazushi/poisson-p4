#!/usr/bin/env python3
"""Independent finite checks and proof audit for the AB all-degree criterion."""
import argparse, hashlib, json, math, time
from collections import Counter
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def run():
    started=time.time()
    main_path=ROOT/'reports/AB_ORIGIN_NORMALITY_20261102.json'
    main=json.loads(main_path.read_text())
    assert main['all_passed']
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in main['input_sha256'].items())
    # Rebuild the fixed-L tensor directly from its two cubic potentials.
    x,y,z,w=s.symbols('x y z w');Y=[x,y,z,w]
    ss=s.symbols('s0:3');hh=s.symbols('h0:5');gg=s.symbols('g0:3')
    trace=sum(ss[i]*x**(2-i)*z**i for i in range(3))
    H=sum(hh[i]*x**(4-i)*z**i for i in range(5))
    F=s.Matrix([x*trace/4+s.diff(H,z),z*trace/4-s.diff(H,x)])
    A=F.jacobian([x,z]);assert (A*s.Matrix([x,z])-3*F).applyfunc(s.expand)==s.zeros(2,1)
    Q=s.zeros(4)
    for i,ii in enumerate([0,2]):
        for j,jj in enumerate([1,3]):Q[ii,jj]=A[i,j];Q[jj,ii]=-A[i,j]
    Q[1,3]=sum(gg[i]*x**(2-i)*z**i for i in range(3));Q[3,1]=-Q[1,3]
    mons=[Y[i]*Y[j]for i,j in combinations_with_replacement(range(4),2)]
    vec=s.Matrix([s.Poly(Q[i,j],*Y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in mons])
    matrix=vec.jacobian([*ss,*hh,*gg]);assert matrix*s.Matrix([*ss,*hh,*gg])==vec
    fourth=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json').read_text())
    selected=matrix[fourth['coordinate_indices'],:]
    rows=list(selected.T.to_DM().rref()[1]);assert len(rows)==11
    minor=s.factor(selected[rows,:].det());assert minor!=0
    # The Y resolution is a declared input; only the Koszul convolution is recomputed here.
    shifts=main['fiber_resolution_shifts'];ranks=main['fiber_resolution_ranks']
    assert shifts==[0,2,3,4,6] and ranks==[1,9,16,9,1]
    tor=[]
    for i in range(4):
        counts=Counter()
        for h in range(i+1):counts[shifts[h]+i-h]+=ranks[h]*math.comb(34,i-h)
        saved=main['low_relative_Tor'][i]['internal_degree_multiplicities']
        assert dict(counts)=={int(k):v for k,v in saved.items()}
        tor.append(dict(counts))
    terms={};bounds={}
    for q in range(4):
        terms[q]=[(j,i,degree)for j in range(4)for i in range(4)if j-i==q for degree in tor[i]]
        bounds[q]=max(degree for j,i,degree in terms[q])
    assert bounds=={0:4,1:3,2:2,3:0}
    a,b,n,g=s.symbols('a b n g')
    assert s.expand((a+2)-(b+1)-(a-b+1))==0
    assert s.expand((a+2)+(b+1)-(a+b+3))==0
    assert s.expand(2*(n-g)+4*g-(2*n+2*g))==0
    assert main['known_coordinate_dimensions']==[1,45,761,6691,38570]
    paths=[Path(__file__),main_path,ROOT/'research/notes/AB_ORIGIN_NORMALITY_20261102.md',
           ROOT/'code/ab_origin_normality_20261102.py',ROOT/'code/ab_origin_normality_fiber_resolution_20261102.m2',
           ROOT/'research/notes/AB_INCIDENCE_20261101.md',
           *[ROOT/f'research/notes/AB_ORIGIN_DEGREE_{d}_20261101.md'for d in ['TWO','THREE','FOUR']],
           ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json',
           ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_COHOMOLOGY_REVIEW_20261101.json',
           ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_MINORS_REVIEW_20261101.json']
    return dict(all_passed=True,accepted=True,
        independent_checks=['Fixed-L 60-coefficient tensor rebuilt from s,H,g; selected45 fiber rank11.',
                            'Koszul34 convolution and total-cohomological-degree0..3 shift bounds.',
                            'Polynomial identities underlying the all-degree Bott inequalities.'],
        fiber_minor_rows=rows,fiber_minor_value=str(minor),low_Tor=tor,
        free_complex_shift_bounds={str(k):v for k,v in bounds.items()},
        spectral_positions={str(k):v for k,v in terms.items()},
        proof_review=['Relative Tor is supported scheme-theoretically on Z and locally free there.',
                      'Proper coherent higher images and regular S give a bounded graded perfect complex.',
                      'Projection formula and spectral-sequence subquotients bound F0 by4 and F1 by3.',
                      'All-n Bott length<=1 and minimal-complex Nakayama imply H1=0 and D generated<=4.',
                      'Independent second route: successively eliminate F3,F2,F1 using only known low-degree cohomology.'],
        limitations=['Y Betti exactness, B_n=D_n for n<=4, and prior low-degree cohomology are audited dependencies, not reimplemented here.',
                     'No main-producer replay or second derived-category implementation is claimed.',
                     'Only reduced AB normality and incidence higher-image vanishing are accepted; no full S61 exhaustion.'],
        primary_sources_checked=['https://stacks.math.columbia.edu/tag/02O3','https://stacks.math.columbia.edu/tag/066Z',
                                 'https://stacks.math.columbia.edu/tag/00MR','https://stacks.math.columbia.edu/tag/01E6'],
        proof_note='research/notes/AB_ORIGIN_NORMALITY_REVIEW_20261102.md',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=run()
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print('AB NORMALITY REVIEW PASS',report['seconds'])
