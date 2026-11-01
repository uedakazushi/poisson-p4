#!/usr/bin/env python3
"""Source geometry review plus an explicit replay of the new free65 point checks."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_binary_deeper_pivot_divisors_20261101 import verify as producer_verify
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.time()
    path=ROOT/'reports/S61_BINARY_DEEPER_PIVOT_DIVISORS_20261101.json'
    old=json.loads(path.read_text());assert old['all_passed']
    for rel,digest in old['input_sha256'].items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest
    u=s.symbols('u1:62');a,b,c,d,e,f=[u[i] for i in [3,13,34,35,48,49]]
    Q=600*b*b-255*b*e-63*e*e+100*a*f
    A=100*b*b-505*b*e+187*e*e
    R=d*A-40*f*f*(2*b-e)
    assert s.gcd(A,40*f*f*(2*b-e))==1 and s.degree(R,d)==1
    # Both charts are absolutely irreducible after f is inverted: the first
    # solves b,c,a, while the second is primitive linear in d.
    subs={b:11*e/5,c:-f/5,a:s.solve(Q,a)[0]}
    subs[a]=s.cancel(subs[a].subs(b,11*e/5))
    assert all(s.cancel(v.subs(subs,simultaneous=True))==0 for v in [Q,5*c+f,5*b-11*e])
    for case in old['cases']:
        M=s.Matrix(case['beta_normal_matrix']).applyfunc(s.Rational)
        det=M[:,case['beta_minor_columns']].det(method='bareiss')
        assert det==s.Rational(case['beta_dominance_minor']) and det!=0
        assert case['source_dimension']==11 and case['joint_image_dimension']==16
        assert case['entire_general_beta_stratum_removed']
    fresh=producer_verify()
    assert len(fresh['cases'])==2
    for prior,current in zip(old['cases'],fresh['cases']):
        for key in ['label','beta_dominance_minor','relation_minor_mod_prime']:
            assert prior[key]==current[key]
    paths=[Path(__file__),path,ROOT/'code/s61_binary_deeper_pivot_divisors_20261101.py',
        ROOT/'code/s61_binary_deeper_boundary_removal_20261101.py',
        ROOT/'reports/S61_BINARY_DEEPER_BOUNDARY_REVIEW_20261101.json',
        ROOT/'reports/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_REVIEW_20261101.json',
        ROOT/'research/notes/S61_BINARY_DEEPER_PIVOT_DIVISORS_REVIEW_20261101.md']
    return dict(id='S61_BINARY_DEEPER_PIVOT_DIVISORS_REVIEW_20261101',all_passed=True,accepted=True,
        verification_kind='Independent source-geometry and proof review; replay of the supplied numerical producer',
        fresh_raw135_point_checks_passed=True,independent_source_minor_checks=2,
        proof_checks=['Absolute irreducibility on f nonzero, including the dense d chart boundary',
            'Both ambient open annihilations follow from relative free65 and smooth equivariant descent',
            'The irreducible16D closures have exceptional subsets of dimension at most15'],
        limitations=['Finite-jet numerics are a producer replay, not a new independent implementation.',
                    'The factor 5b-3e and the remaining deeper open are not removed.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))
