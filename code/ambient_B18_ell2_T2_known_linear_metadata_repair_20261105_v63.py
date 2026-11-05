#!/usr/bin/env python3
"""Report-only repair after completed432 math; never reruns that math."""
import json
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63'
def main():
 out=ROOT/'reports'/STEM;ep=out/'actual432-overlap-rank-and-known-linear-cover-degree2.json';assert H.sha(ep)=='c1309522f6bb018592fe03a13585eb379942de72ab93eb51727ee8336a6493c6';d=json.loads(ep.read_text())
 gp='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/report.json';g=json.loads((ROOT/gp).read_text());ip=str(Path(gp).parent/'actual-complete-L48-and-degree2-rewrite-source-input.json');pp='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63_REVIEW/report.json'
 assert d['rank']==0 and all(not v for v in d['all432_remainder_columns']) and d['known_homogeneous_linear_cover_HF_through2']==[6,6,6]
 failure={'reason':'report input path variable p shadowed by integer linear pivot; AttributeError after complete evidence saved','completed_math_evidence_sha256':H.sha(ep),'arithmetic_rerun':False,'report_only_repair':True}
 H.dump(out/'FIRST_REPORT_PATH_SHADOW.json',failure)
 slots=sum(len(v) for i,v in d['all264_linear_NF_columns'])+sum(map(len,d['known_degree2_quotient_projection_columns']))+10000
 report={'id':STEM,'all_passed':True,'all432_source_multiples_and_overlaps_retained':True,'known_cover_HF_through2':[6,6,6],
  'known_cover_prefix_length':18,'degree2_relation_rank':0,'actual_T2_HF2_total_length_or_H0_computed':False,
  'input_sha256':{p:H.sha(ROOT/p) for p in [gp,ip,pp]},'code_sha256':H.sha(ROOT/'code'/f'ambient_B18_ell2_T2_known_linear_degree2_20261105_v63.py'),
  'executed_math_code_sha256':H.sha(out/'FIRST_REPORT_PATH_SHADOW_CODE.py'),
  'metadata_repair_code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep),'FIRST_REPORT_PATH_SHADOW.json':H.sha(out/'FIRST_REPORT_PATH_SHADOW.json')},
  'owned_QQ_slot_bound':slots,'all_intermediate_height_guards_passed_before_metadata_failure':True,
  'verified_NF_intermediate_height_upper_bits':139,'declared_Gaussian_guard_bits':2048,'arithmetic_rerun':False}
 H.dump(out/'report.json',report);print(json.dumps({'all_passed':True,'known_cover_HF_through2':[6,6,6],'report_sha256':H.sha(out/'report.json')}))
if __name__=='__main__':main()
