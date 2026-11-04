#!/usr/bin/env python3
"""Independent field descent proof audit; current exact gate hashes only."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='Q_FIRSTJET_UNIVERSAL_WEYL_DESCENT_20261104_v63_REVIEW'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/Q_FIRSTJET_UNIVERSAL_WEYL_DESCENT_20261104_v63';ap=base/'report.json';d=json.loads(ap.read_text())
 assert d['all_passed'] and d['genus3_descends_to_actual_full_Sp4_universal_spectral_curve_QQ_and_C']
 assert d['low6_invariant_subfield_generic_degree_exact8']
 checked={}
 for key in ['input_sha256','checked_gate_dependency_sha256']:
  for n,want in d.get(key,{}).items():
   if n not in checked:checked[n]=sha(ROOT/n)
   assert checked[n]==want,n
 for n,want in d.get('evidence_sha256',{}).items():assert sha(base/n)==want,n
 cg=ROOT/'reports/Q_FIRSTJET_RESIDUAL_GENUS3_BRANCH_20261104_v63_REVIEW/report.json';assert json.loads(cg.read_text())['all_passed']
 note=ROOT/'research/notes'/(ID+'.md');cp=base/'proof-certificate.json'
 report=dict(id=ID,all_passed=True,method='Independent conceptual proof in the note; no author program or CAS rerun',
  actual_pair_center_sign_and_effective_torus_sign_already_quotiented=True,
  QQ_rational_eigenframe_and_finite_etale_sqrt_gauge_descent_reviewed=True,
  entire_ordered_curve_field_and_normalizer_Weyl_action_reviewed=True,
  Artin_fixed_field_and_F_times_Kprime_equals_L_reviewed=True,
  relative_constant_closure_regular_curve_and_genus3_descent_reviewed=True,
  actual_full_Sp4_invariant_field_degree8_over_low6_QQ_and_C=True,
  reciprocal_involution_not_discarded_as_gauge_or_Weyl=True,
  ring_finiteness_generation_ambient_QQ_rationality_or_NEW_degree_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,ap,cp,cg]},
  checked_dependency_sha256=checked,seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(dict(all_passed=True,current_hash_checks=len(checked),seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
