#!/usr/bin/env python3
"""Small exact dimensional corollary, no new coefficient computation or point."""
import hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H_LINEAR2_RANK126_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 rb=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';r=read(rb/'report.json');assert r['total_proved_H_linear_unit_lower']>=2 and r['actual16_raw_H_inputs_source_bound']
 hb=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';h=read(hb/'certificate.json');hr=read(hb/'report.json');assert sha(hb/'certificate.json')==hr['certificate_sha256'];assert hr['global_H2_knownunit_sourcepair_and_deletion_verified']
 review=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW/report.json';assert read(review)['all_passed'];assert h['target117_projection_preserves_Hdegree_exact_structural_support']
 mon=h['parent_joint_target239_monomials'];assert sum(m[7]==1 for m in mon)==28 and len(h['target122_Hdegree_classes']['1'])==23
 assert len(h['all112_genuine_cubic_H1_scalar_roots'])==112
 span=28-15;quotient_loss=28-23;assert span==13 and quotient_loss==5;lower=117+1+(span-quotient_loss);assert lower==126
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rb/'report.json',rb/'preflight.json',rb/'checkpoint.json',hb/'report.json',hb/'certificate.json',review]
 report={'id':ID,'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','raw_Hlinear_unit_lower_primary':2,'two_linearforms_times7light_dimension':13,'NF117_H1_projection_dimension_loss':5,'pointwise_H1_cubic_image_lower':8,'full319x239_pointwise_rank_lower_primary':126,
  'same_single_global126_minor_or_eighth_constructiveH_pivot_NOT_claimed':True,'prior_constructive_unit_lower125_unchanged':True,'no_source_free_basis_count_change_or_full239_claim':True,'independent_rawH2_review_required_before_ledger_promotion':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'seconds':time.monotonic()-t,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
