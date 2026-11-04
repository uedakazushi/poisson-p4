#!/usr/bin/env python3
"""Serialization-only completion of the exact already-passed contact review.
Does not evaluate polynomials, matrices or sources.
"""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_CRITICAL_V1_CONTACT_20261104_v63_REVIEW'
out=ROOT/'reports'/STEM
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    ev=out/'independent-N-adaptation-contact-and-local-proof.json'
    assert sha(ev)=='756abfd67c8f63ff77135039f9c6080948ca8fc885fdf342bfaaf064d34de8db'
    E=read(ev);assert E['all3300_actual_source_slots_verified'] and E['all150_ambient_PiDelta_zero'] and E['generic_CA_critical_kappa']==1
    bad=out.with_name(out.name+'_METADATA_ATTEMPT');failure=read(bad/'failure.json')
    assert failure['source_arithmetic_completed'] and failure['completed_evidence_sha256']==sha(ev) and not failure['arithmetic_rerun_for_metadata_repair']
    files=[Path(__file__),ROOT/'code/ambient_CA_critical_V1_contact_review_20261104_v63.py',ROOT/'research/notes'/f'{STEM}.md',
      ROOT/'code/ambient_CA_critical_V1_relative_tangent_review_20261104_v63.py',ROOT/'code/ambient_B18_fullV56_comparison_unit_review_20261104_v63.py']
    hashes=['b30b7e43b30e503fb8ec1ce0ac541288ab723178ed495a0d620e22762613fb8c','06cafbe66b9022bb37cc82e505e1a711f16ea8084f28b93f68e4786317cf614d','74d00c4095a71da7f42c86e5d9146823621f0c7491a5432b0964489f6a9e4cf8']
    for name,h in zip(['CUSP_DIRECTION_PREFLIGHT','CONTACT_OPERATOR_PREFLIGHT','CONTACT_VECTOR'],hashes):
        folder=ROOT/'reports'/('AMBIENT_CA_CRITICAL_V1_'+name+'_20261104_v63');p=folder/'report.json';r=read(p);assert sha(p)==h and r['all_passed'];files.append(p)
        for rel,hh in r['input_sha256'].items():assert sha(ROOT/rel)==hh;files.append(ROOT/rel)
        for rel,hh in r['evidence_sha256'].items():assert sha(folder/rel)==hh;files.append(folder/rel)
    for d in [out.with_name(out.name+'_FIRST_ASSERTION_ATTEMPT'),bad]:files+=list(d.iterdir())
    (out/'phase.json').write_text(json.dumps({'phase':'FINAL_all_exact_arithmetic_already_passed_metadata_only_completed','completed_evidence_sha256':sha(ev),'arithmetic_rerun':False})+'\n')
    report={'id':STEM,'all_passed':True,'independent_actual_N_frame_adjoint_source_replay_and_contact_proof':True,
      'first_reviewer_Wc_equality_mistake_preserved_and_exact_tangent_correction_verified':True,
      'all3300_source_slots_and_all150_graph_bits_checked':True,'c_v':'0','nonzero_graph_count':0,
      'actual_generic_CA_critical_kappa':1,'transverse_normalizer_kernel_length':1,'transverse_comparison_C_length':1,
      'no_authorproducer_newpoint_jetorder_source_solve_GB':True,'origin_D0_depth_socle_fullnil_claimed':False,
      'metadata_only_repair_of_completed_math':True,'completed_math_seconds':0.246442502,'max_rss_kib_not_recorded':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{ev.name:sha(ev),'phase.json':sha(out/'phase.json')}}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'id':STEM,'all_passed':True,'report_sha256':sha(out/'report.json'),'arithmetic_rerun':False}))
if __name__=='__main__':main()
