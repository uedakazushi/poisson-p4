#!/usr/bin/env python3
"""Finite Tor-shift certificate supporting the all-degree AB normality proof."""
import argparse,hashlib,json,math,os,subprocess,time
from pathlib import Path
from ab_origin_degree_four_backend_20261101 import *

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.time()
    stem='code/ab_origin_normality_fiber_resolution_20261102.m2'
    env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    proc=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(ROOT/stem)],
                        cwd=ROOT,env=env,text=True,capture_output=True,check=True,timeout=30)
    assert 'AB_Y_BETTI_EXACT_QQ_PASS' in proc.stdout
    degrees=[0,2,3,4,6];ranks=[1,9,16,9,1];tor=[]
    for i in range(4):
        counts=Counter()
        for j in range(min(i,4)+1):counts[degrees[j]+i-j]+=ranks[j]*math.comb(34,i-j)
        assert max(counts)<= (0 if i==0 else i+1)
        tor.append(dict(homological_degree=i,internal_degree_multiplicities={str(d):v for d,v in counts.items()}))
    bounds={}
    for q in [0,1]:
        possible=[d for i in range(4)for j in range(4)if j-i==q for d in map(int,tor[i]['internal_degree_multiplicities'])]
        bounds[q]=max(possible)
    assert bounds=={0:4,1:3}
    # The fiber embedding is linear of rank11 already in the chosen 45 coordinates.
    g0,g1,g2=s.symbols('g0 g1 g2');fiber=(*pars,g0,g1,g2)
    Q=Qof(N,g0*x*x+g1*x*z+g2*z*z)
    spatial=[Y[i]*Y[j]for i,j in combinations_with_replacement(range(4),2)]
    coeff=s.Matrix([s.Poly(Q[i,j],*Y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in spatial])
    A=coeff.jacobian(fiber);assert A*s.Matrix(fiber)==coeff
    fourth=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_20261101.json').read_text())
    selected=fourth['coordinate_indices'];A=A[selected,:];rows=list(A.T.to_DM().rref()[1]);assert len(rows)==11
    determinant=A[rows,:].det();assert determinant.is_Rational and determinant
    # The low H1 vanishing not requiring the degree2/3 parabolic argument.
    low=[]
    for n in [0,1]:
        hc=Counter()
        for j in range(n+1):
            cy=character(n-j,core_weights,leading);cv=character(j,vertical_weights);prod=Counter()
            for a,ca in cy.items():
                for b,cb in cv.items():prod[a[0]+b[0],a[1]+b[1]]+=ca*cb
            for a,b,m in decompose(prod):
                out=bott(a,b)
                if out:hc[out[0]]+=m*out[2]
        assert not any(v for q,v in hc.items()if q>0)
        low.append(dict(degree=n,cohomology_dimensions=dict(hc)))
    second=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_TWO_20261101.json').read_text())
    third=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_THREE_20261101.json').read_text())
    upper=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_COHOMOLOGY_REVIEW_20261101.json').read_text())
    minor=json.loads((ROOT/'reports/AB_ORIGIN_DEGREE_FOUR_MINORS_REVIEW_20261101.json').read_text())
    assert all(r.get('all_passed')or r.get('accepted')for r in [second,third,fourth,upper,minor])
    assert second['normalization_degree_two_H1']==third['normalization_degree_three_H1']==0
    assert second['normalization_degree_two_dimension']==761 and third['normalization_degree_three_dimension']==6691
    assert fourth['coordinate_degree_four_dimension']==fourth['normalization_degree_four_dimension']==38570
    paths=[str(Path(__file__).relative_to(ROOT)),stem,'code/ab_origin_degree_four_backend_20261101.py',
           'research/notes/AB_ORIGIN_NORMALITY_20261102.md','research/notes/AB_INCIDENCE_20261101.md',
           'research/notes/AB_ORIGIN_DEGREE_TWO_20261101.md','research/notes/AB_ORIGIN_DEGREE_THREE_20261101.md',
           'research/notes/AB_ORIGIN_DEGREE_FOUR_20261101.md',
           *[f'reports/{x}_20261101.json'for x in ['AB_ORIGIN_DEGREE_TWO','AB_ORIGIN_DEGREE_THREE','AB_ORIGIN_DEGREE_FOUR','AB_ORIGIN_DEGREE_FOUR_COHOMOLOGY_REVIEW','AB_ORIGIN_DEGREE_FOUR_MINORS_REVIEW']]]
    report=dict(all_passed=True,conclusion='Reduced AB is normal everywhere, including the origin; all higher incidence direct images vanish.',
                ambient_coordinate_count=45,linear_fiber_rank=11,linear_kernel_rank=34,
                fiber_minor_rows=rows,fiber_minor_value=str(determinant),
                fiber_resolution_ranks=ranks,fiber_resolution_shifts=degrees,fiber_resolution_exactness_checked=True,
                low_relative_Tor=tor,derived_free_generator_degree_bounds={str(k):v for k,v in bounds.items()},
                low_cohomology=low,known_coordinate_dimensions=[1,45,761,6691,38570],
                all_degree_argument='See proof note: all-n Bott length<=1; minimal perfect complex; H1 generated<=3; D generated<=4.',
                macaulay2_stdout=proc.stdout,seconds=time.time()-started,
                proof_note='research/notes/AB_ORIGIN_NORMALITY_20261102.md',
                limitations=['The all-degree statement uses the geometric and derived proof, not finite dimension agreement alone.',
                             'No scheme-theoretic T01 zero-fiber description, full S61 exhaustion, or component upper bound.'],
                primary_foundations=['https://stacks.math.columbia.edu/tag/02O3','https://stacks.math.columbia.edu/tag/066Z','https://stacks.math.columbia.edu/tag/00MR','https://stacks.math.columbia.edu/tag/01E6'],
                input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    with args.output.open('x')as out:json.dump(report,out,indent=2);out.write('\n')
    print('AB ALL-DEGREE NORMALITY CERTIFICATE PASS',report['seconds'],flush=True)

if __name__=='__main__':main()
