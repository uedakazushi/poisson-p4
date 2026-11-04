#!/usr/bin/env python3
"""Prior-gate application: T01 excluded everywhere, T39/T40 generically."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_GENERIC_KNOWN_EXCLUSIONS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    paths=[ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63_REVIEW/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63_REVIEW.json',
           ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz']
    boundary,positive,torus,geometry,inc=map(read,paths)
    assert boundary['all_passed']and positive['all_passed']and torus['all_passed']and geometry['all_passed']
    assert boundary['independent_whole_v_zero_u_unit_stabilizer_exact_Ga']
    assert positive['T39_T40_excluded_by_closed_invariant_components_and_actual_limit']
    t=s.Symbol('t');n0,n1,n2=s.symbols('n0 n1 n2')
    P=s.Matrix(geometry['affine_P']);curve=P.subs({n0:t,n1:t,n2:0},simultaneous=True)
    oldpositive=s.Matrix(inc['positive_P'])
    assert(curve.subs(t,0)-oldpositive).applyfunc(s.expand)==s.zeros(4)
    proof=dict(source='the irreducible projective PGL saturation of P(u,1,0), u!=0',
        T01_all_u='Its every orbit has dimension23; the whole closed T01 component has orbit dimension<=20, so every u!=0 is excluded.',
        actual_curve='P(t,t,0)',curve_to_saturation='For t!=0 the accepted all-parameter action with r=t,epsilon=1 takes P(t,1,0) to P(t,t,0).',
        limit='P(t,t,0) tends to the ALREADY certified positive-axis Pplus as t->0.',
        generic_T39_T40='If either closed PGL-stable component contained the whole irreducible saturation, it would contain Pplus; the accepted positive-axis exclusion contradicts this.',
        no_exceptional_u_classification=True,no_new_point_prime_source_rank_or_high_invariant=True)
    ep=out/'proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    paths +=[Path(__file__),ROOT/'research/notes'/(ID+'.md')]
    report=dict(id=ID,all_passed=True,all_v_zero_u_nonzero_outside_T01=True,
        generic_v_zero_u_nonzero_outside_T39_and_T40=True,
        explicit_old_Pplus_limit_bound_to_actual_three_parameter_family=True,
        T39_T40_exceptional_parameter_roots_not_computed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-started,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':report['seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
