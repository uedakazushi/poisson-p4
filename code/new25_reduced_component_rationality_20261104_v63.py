#!/usr/bin/env python3
"""Combine the full S61 rational-fibre proof and the exact plus-nine product."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='NEW25_REDUCED_COMPONENT_RATIONALITY_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 stems=['Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW',
        'S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63_REVIEW']
 gates=[ROOT/'reports'/n/'report.json'for n in stems]
 records=[json.loads(p.read_text())for p in gates]
 assert all(d['all_passed']for d in records)
 # The note, rather than these existence checks alone, supplies the application.
 checked={}
 for d in records:
  for key in ['input_sha256','checked_dependency_sha256']:
   for name,want in d.get(key,{}).items():
    if name not in checked:checked[name]=sha(ROOT/name)
    assert checked[name]==want,name
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,
  target='Reduced projective Poisson component R4_POSITIVE_WEIGHT6_25, not only its firstjet image',
  full_S61_Y16_rational_over_QQ_and_C_by_same_base_Luroth=True,
  reduced_projective_X25_birational_Y16_times_A9=True,
  reduced_projective_new_component_pure_rational_dimension25=True,
  coefficient_fields=['QQ','C'],
  one_variable_source_degree1_or_expanded_Luroth_generator_assumed=False,
  Jacobi_scheme_reducedness_or_generic_length5_changed=False,
  NEW_firstjet_equation_or_global_component_exhaustion_claimed=False,
  no_new_CAS_point_prime_invariant_or_source26_computation=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,*gates]},
  checked_dependency_sha256=checked,seconds=time.monotonic()-start,
  independent_combined_application_review_pending=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(all_passed=True,checked_dependencies=len(checked),seconds=result['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
