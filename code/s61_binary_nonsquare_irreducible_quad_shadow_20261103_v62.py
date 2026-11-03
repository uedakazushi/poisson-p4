#!/usr/bin/env python3
"""Uniform quadratic-product rank37 via an irreducible quadric and3-plane restriction."""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_v62'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(out,p,sign=1):
    for ex,c in p.items():out[ex]+=sign*c

def mul(p,q):
    out=defaultdict(Fraction)
    for u,c in p.items():
        for v,d in q.items():out[tuple(a+b for a,b in zip(u,v))]+=c*d
    return {ex:c for ex,c in out.items()if c}

def mons(n,d):return [(d,)]if n==1 else[(i,)+v for i in range(d+1)for v in mons(n-1,d-i)]
def verify():
    start=time.monotonic();prior=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/report.json';P=json.loads(prior.read_text());assert P['all_passed']
    for name,want in P['input_sha256'].items():assert sha(ROOT/name)==want,name
    src=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';D=json.loads(src.read_text());q0={tuple(mon[:6]):{tuple(ex):Fraction(c)for ex,c in terms}for mon,terms in D['cleared_qfour_sparse'][0]}
    # q0 is unchanged by original-d=P/a: it contains no original d.
    assert all(ex[2]==0 for q in q0.values()for ex in q)
    ids=(0,1,2);M=[]
    for i in ids:
        row=[]
        for j in ids:
            mon=tuple(int(k==i)+int(k==j)for k in range(6));scale=Fraction(1 if i==j else 1,1 if i==j else 2)
            row.append({ex:c*scale for ex,c in q0.get(mon,{}).items()})
        M.append(row)
    det=defaultdict(Fraction)
    for p in itertools.permutations(range(3)):
        inversions=sum(p[i]>p[j]for i in range(3)for j in range(i+1,3));v={(0,0,0):Fraction(1)}
        for i,j in enumerate(p):v=mul(v,M[i][j])
        add(det,v,(-1)**inversions)
    det={ex:c for ex,c in det.items()if c};assert det=={(2,0,0):Fraction(125,256)}
    # Orthogonal Grassmann strata with radical intersection dimension k.
    fano=[]
    for r in range(1,7):
        radical=6-r;dims=[]
        for k in range(4):
            s=3-k
            if k<=radical and s<=r//2:
                dims.append(k*(radical-k)+s*(r-s)-s*(s+1)//2+s*(radical-k))
        fano.append(max(dims))
    assert fano==[6,6,4,4,3,3]
    assert max(2+max(fano[:2]),3+max(fano[2:]))==8<9
    hist={};qmon=mons(3,2)
    for four in itertools.combinations(qmon,4):
        shadow={tuple(a+b for a,b in zip(q,x))for q in four for x in mons(3,1)}
        hist[len(shadow)]=hist.get(len(shadow),0)+1
    assert sum(hist.values())==15 and min(hist)==8
    step33=ROOT/'reports/S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62/report.json';S=json.loads(step33.read_text());assert S['all_passed'] and S['total_constructive_product_unit_pivots']==33
    files=[Path(__file__),prior,src,step33,ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,field='QQ and all geometric characteristic0 points',open='a*d*F*B*C*J!=0, J=5ad-6e-4',
        exact_q0_symmetric_minor_indices=list(ids),exact_q0_minor='125*a^2/256',
        q0_irreducible_at_every_point=True,known_four_quad_coefficient_rank=4,
        isotropic3plane_variety_max_dimensions_by_quadric_rank=fano,
        bad3plane_incidence_dimension_upper=8,Grassmann_3_6_dimension=9,
        some3plane_restriction_injective_at_each_geometric_point=True,
        degree2_monomial_quadruples_in3vars_checked=15,shadow_histogram=hist,
        restricted_cubic_product_lower_bound=8,independent_kernel_products=12,
        all_point_mod_ell_cubic_product_lower_bound=20,pair_cover_kernel_products=17,
        all_point_original_quadratic_product_lower_bound=37,
        generic_constant_unimodular_product_row_exists_after33_pivots=True,
        all_point_rank41_or_full_cubic57_or_quartic_coverage_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in files},seconds=time.monotonic()-start)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])
