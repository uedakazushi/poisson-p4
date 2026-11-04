#!/usr/bin/env python3
"""Bind the pure length3 application to the independent actual joint105 audit."""
import hashlib,json,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    start=time.monotonic()
    base=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63'
    primary=json.loads((base/'report.json').read_text());cert=json.loads((base/'certificate.json').read_text())
    peer=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json'
    own=json.loads(peer.read_text());assert primary['all_passed']and own['all_passed']
    for rel,h in primary['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in primary['evidence_sha256'].items():assert sha(base/rel)==h,rel
    independent=json.loads((peer.parent/'independent-certificate.json').read_text())
    assert cert['exact_quadratic_coefficient_matrix']==independent['quadratic_matrix']
    assert F(cert['determinant'])==F(independent['quadratic_minor'])==F(625,768)
    assert own['exact_joint_zero_ideal_quotient_length']==3
    assert own['symmetric_intrinsic_T3_annihilation_and_reverse_weight3_stability_accepted']
    assert cert['quotient_length']==3 and cert['quotient_Hilbert_function']==[1,2]
    assert primary['T3_length_upper_bound']==36 and not primary['T3_nonzero_depth19_colon_or_actual_length_claimed']
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',peer,peer.parent/'independent-certificate.json',
           ROOT/'research/notes/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63.md']
    out=ROOT/'reports'/(ID+'.json');assert not out.exists()
    result=dict(id=ID,all_passed=True,actual_joint_ideal_equals_z7_plus_r_s_squared=True,
        exact_joint_quotient_length=3,T3_finite_length=True,T3_length_upper=36,maximal_ideal_squared_annihilates_T3=True,
        proof='The independent actual105 audit verifies F subset m^2 and F+m^3=m^2, so Nakayama gives F=m^2. Both actual zero ideals annihilate the intrinsic T3; its12 generators give length<=36.',
        T3_nonzero_actual_length_socle_depth19_or_full_nil_not_claimed=True,
        no_new_jets_point_GB_or_author_rerun=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start,frozen=True)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':result['seconds']}))
if __name__=='__main__':main()
