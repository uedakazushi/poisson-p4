#!/usr/bin/env python3
"""Record independent conceptual review and bind all three actual proof gates."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63_REVIEW'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63'
 ap=base/'report.json';author=json.loads(ap.read_text())
 assert author['all_three_final_independent_gates_passed']
 cp=base/'application-certificate.json';assert sha(cp)==author['application_certificate_sha256']
 checked={}
 for field in ['input_sha256','checked_gate_dependency_sha256']:
  for p,want in author[field].items():
   if p not in checked:checked[p]=sha(ROOT/p)
   assert checked[p]==want,p
 cert=json.loads(cp.read_text());assert cert['parameter_order']==['xi','r','s','q','B']
 assert cert['group_action'].startswith('Actual rho is finite pullback')
 refs=[ROOT/'reports'/n/'report.json' for n in [
  'Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW',
  'Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW',
  'Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW']]
 for p in refs:assert json.loads(p.read_text())['all_passed']
 note=ROOT/'research/notes'/(ID+'.md')
 report=dict(id=ID,all_passed=True,method='Independent mathematical review in the note; current gate hashes only, no author rerun',
  same_actual_quotient_pullback_and_section_identity_reviewed=True,
  dense_saturation_excludes_pole_and_nonfree_locus_reviewed=True,
  full_geometric_stabilizer_not_only_Lie_algebra_reviewed=True,
  dominant_equal_dimension_separable_degree1_over_QQ_reviewed=True,
  torsor_trivialized_by_actual_section_not_quotient_rationality_assumption=True,
  pullback_anti_action_and_right_source_multiplication_reviewed=True,
  independent_Cayley_rational10_argument_in_note=True,
  B_NEW15_pure_rational_over_QQ_and_C=True,full25_component_rationality_claim=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),note,ap,cp,*refs]},
  checked_dependency_sha256=checked,seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(dict(all_passed=True,current_hash_checks=len(checked),seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
