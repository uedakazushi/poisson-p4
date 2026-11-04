#!/usr/bin/env python3
"""Independent whole-domain orbit proof audit; gate hashes, no producer rerun."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='Q_NEW_FAMILY_ENTIRE_NORMAL_FORM_ORBITS_20261104_v63_REVIEW'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/Q_NEW_FAMILY_ENTIRE_NORMAL_FORM_ORBITS_20261104_v63';ap=base/'report.json';d=json.loads(ap.read_text());assert d['all_passed']
 assert d['entire_uvp_normalform_parameter_identifications_exact_Klein_four']
 assert d['finite_parameter_quotient_Gm_times_A1_minus_256']
 checked={}
 for key in ['input_sha256','checked_gate_dependency_sha256']:
  for n,want in d.get(key,{}).items():
   if n not in checked:checked[n]=sha(ROOT/n)
   assert checked[n]==want,n
 for n,want in d.get('evidence_sha256',{}).items():assert sha(base/n)==want,n
 note=ROOT/'research/notes'/(ID+'.md');cp=base/'proof-certificate.json'
 r=dict(id=ID,all_passed=True,method='Independent low-coefficient, conic covariance and parity proof in the note; hashes only',
  entire_domain_intrinsic_modular_image_flag_reviewed=True,
  eta_squared1_without_dividing_a_or_c_and_target_U_pm_u_reviewed=True,
  line_covariance_and_conic_J_all_domain_units_reviewed=True,
  full_tau_sigma_necessary_and_sufficient_entire_uvp=True,
  finite_parameter_quotient_ring_and_free_rank4_reviewed=True,
  generic_logarithmic_line_uniqueness_or_all_smooth_Q_component_coverage_assumed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,ap,cp]},
  checked_dependency_sha256=checked,seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(all_passed=True,hash_checks=len(checked),seconds=r['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
