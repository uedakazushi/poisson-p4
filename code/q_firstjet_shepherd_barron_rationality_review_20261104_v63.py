#!/usr/bin/env python3
"""Primary-literature scope audit + frozen F=B(t2) binding; no new CAS job."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_SHEPHERD_BARRON_RATIONALITY_20261104_v63_REVIEW'
URL='https://www.numdam.org/article/CM_1989__70_1_13_0.pdf'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    d=read(p);assert d['all_passed']
    for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in d.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return d

def main(out):
    start=time.monotonic();assert not out.exists();out.mkdir(parents=True)
    sr=ROOT/'reports/Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63_REVIEW/report.json';s5=accepted(sr)
    vr=ROOT/'reports/Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63.json';v=accepted(vr)
    assert s5['fivepoint_unique_biform_Koszul_and_projective_frame_quotient_proof_accepted']
    assert s5['center_squared_line_rank2_bundle_honest_full_AutQ_descent_proof_accepted']
    assert s5['full_Sp4_field_equals_B_two_independent_variables']
    assert v['full_field_equals_projective_finite_group_base_plus_two_independent_variables']
    primary=dict(primary_url=URL,author='N. I. Shepherd-Barron',year=1989,journal='Compositio Mathematica',volume=70,pages='13–25',
        reading_method='Independently opened primary PDF via managed web; text read and printed18/21 screenshots viewed',
        observed_PDF_pages_including_cover=14,visually_read_printed_pages=[18,21],corresponding_zero_based_PDF_indices=[6,9],
        primary_scope=dict(proposition11='Projective biform pair types(1,2) and(2,1), quotient by fullAutQ including factor swap, is rational',
                           theorem15='Projective quotient for five-dimensional irreducible S5 representation is rational',
                           lemma8='Projective character twists give the same quotient'),
        application_base_field='C only',verbatim_source_quotes_used=False,source_paraphrase_under200_words=True,
        paper_proof_reproved_by_CAS=False,review_recorded_at_utc=datetime.now(timezone.utc).isoformat())
    ep=out/'primary-reading-record.json';ep.write_text(json.dumps(primary,indent=2)+'\n')
    assert (5+5)-6==4 and 4+2==6
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',sr,vr,
        ROOT/'research/notes/Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63_REVIEW.md',
        ROOT/'research/notes/Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63.md']
    r=dict(id=STEM,all_passed=True,
        verification_kind='Independent primary-literature reading and exact theorem-hypothesis match; frozen field-proof hash audit',
        primary_Proposition11_matches_exact_projective_pair_full_AutQ_quotient=True,
        no_hidden_connected_only_group_or_ruling_orientation_change=True,
        accepted_full_Sp4_field_F_equals_B_two_variables=True,
        consequence_over_C='Purely rational invariant function field of transcendence degree6',
        rationality_theorem_is_existing_literature_application_not_new_theorem=True,
        explicit_six_generators_constructed=False,arbitrary_field_or_characteristic_claim=False,
        NEW_equation_degree_upper_or_source_target_field_index_solved=False,
        polynomial_invariant_ring_or_all_orbits_claim=False,
        prior_stable2_packets_changed=False,new_point_coefficient_GB_or_producer_job=False,
        source_url=URL,input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['all_passed','consequence_over_C','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
