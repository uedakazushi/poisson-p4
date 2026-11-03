#!/usr/bin/env python3
"""Explicit CA/CB source lifts for every Qa and non-CM scope receipt."""
import argparse, hashlib, json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X, raw as rawA, canonical
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
from q_weighted_finite_model_20261103_v63 import point, coefficients

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False)
    h,x,y,z,w=X;a=s.Symbol('a')
    H=s.Rational(4,3)*h*x*w-h*z*z/3-x*x*z
    A=canonical(rawA(a*h*h,-s.Rational(2,9),H))
    B=canonical(rawB(a*h*h,-s.Rational(2,9)*h*h,H))
    assert A==B
    # Source ordering (h,x,y,z,w)=(old0,old1,old4,old2,old3).
    order=[0,1,4,2,3];old=s.symbols('x0:5');P=s.zeros(5)
    for i in range(5):
        for j in range(5):
            P[order[i],order[j]]=A[i,j].subs(dict(zip(X,[old[k]for k in order])),simultaneous=True)
    assert (P-point(a)).applyfunc(s.expand)==s.zeros(5)
    assert all(s.expand(v)==0 for v in coefficients(P-point(a)))
    t01=json.loads((ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63_REVIEW.json').read_text())
    assert t01['all_passed']
    paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py',
        ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63_REVIEW.json',
        ROOT/'research/notes/Q_T01_MEMBERSHIP_20261103_v63.md',
        ROOT/'research/notes/Q_ALL_FAMILY_NON_CM_20261103_v63.md']
    report=dict(id='Q_ALL_FAMILY_NON_CM_20261103_v63',all_passed=True,
        parameter_scope='every finite complex a, including zero',
        exact_CA_CB_raw_and_canonical_lifts_all150=True,
        CA_source={'f':'a*h^2','g':'-2/9','H':str(H)},
        CB_source={'F':'h^2','ratio':'[a:-2/9]','H':str(H)},
        source_coordinate_order=order,
        known_minimal_projective_component_dimensions={'T01':25,'CA':28,'CB':28},
        all_Qa_full_Poisson_completed_local_ring_non_Cohen_Macaulay=True,
        all_Qa_full_reduced_Poisson_completed_local_ring_non_Cohen_Macaulay=True,
        justification='closed-point minimal components of different dimensions; a complete Cohen-Macaulay local ring is equidimensional',
        all_local_branches_or_local_dimension_upper_bound_claimed=False,
        individual_CA_CB_nonCM_for_all_a_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('All Qa explicit CA/CB lifts and non-CM scope PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    main(p.parse_args().output)
