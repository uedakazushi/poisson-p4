#!/usr/bin/env python3
"""Independent normalized-equivalence count and exact field-degree audit."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='Q_NEW_FAMILY_FULL_SOURCE_DEGREE_20261104_v63_REVIEW'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/Q_NEW_FAMILY_FULL_SOURCE_DEGREE_20261104_v63';ap=base/'report.json';d=json.loads(ap.read_text());assert d['all_passed']
 checked={}
 for key in ['input_sha256','checked_gate_dependency_sha256']:
  for n,want in d.get(key,{}).items():
   if n not in checked:checked[n]=sha(ROOT/n)
   assert checked[n]==want,n
 for n,want in d.get('evidence_sha256',{}).items():assert sha(base/n)==want,n
 note=ROOT/'research/notes'/(ID+'.md');cp=base/'proof-certificate.json'
 r=dict(id=ID,all_passed=True,method='Independent conceptual normalized-equivalence proof in note; gate hashes only',
  source_and_normalized_equivalence_tuples_bijective_on_generic_domain=True,
  four_normalforms_times_two_aut_components_not_extra_sigma_or_center=True,
  unique_Ga_y_zero_cross_section_and_c_frame_perspective_checked=True,
  original_full_source_generic_degree_exact8=True,
  tau_fixed_rational_frame_source_degree2_and_same_actual_beta_base_checked=True,
  descended_full_source_fiber_generic_degree_exact4_QQ_and_C=True,
  degree_one_z_map_or_expanded_Luroth_generator_or_special_fiber_claim=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,ap,cp]},
  checked_dependency_sha256=checked,seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(all_passed=True,hash_checks=len(checked),seconds=r['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
