#!/usr/bin/env python3
"""Bind the exact hypotheses of the multiplicity-one reduced-smoothness proof."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 paths=[Path(__file__),ROOT/'research/notes/R4_NEW25_POSITIVE_AXIS_SMOOTH_REDUCED_COMPONENT_20261103_v63.md']
 receipts=['Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63/report.json','Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW/report.json','Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63/report.json','Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/report.json','Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json','Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63.json','Q_NEW_FAMILY_FULL_COMPLETION_ROOT_REVIEW_20261103_v63.json','Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json']
 data={}
 for name in receipts:
  p=ROOT/'reports'/name;d=json.loads(p.read_text());assert d['all_passed'];data[name]=d;paths.append(p)
  for source,want in d.get('input_sha256',{}).items():assert sha(ROOT/source)==want,source
 fifth=data[receipts[0]];assert fifth['through5_ideal_dimension']==3 and fifth['through5_ideal_degree']==5
 assert fifth['radical_unique3plane_free_coordinates']==['z0','z3','z5']
 result=dict(id='R4_NEW25_POSITIVE_AXIS_SMOOTH_REDUCED_COMPONENT_20261103_v63',all_passed=True,proof='Hilbert-Samuel associativity plus generic length5 and actual multiplicity<=5; excellent reduced equidimensional multiplicity-one ring is regular',actual_transverse_dimension=3,full_transverse_multiplicity_exact=5,reduced_new_component_transverse_multiplicity_exact=1,reduced_new_component_smooth_on_entire_nonzero_positive_axis=True,unique_dimension25_minimal_component_at_positive_axis=True,formal_new_component_transverse_tangent_space='L3=(z1,z2,z4,z6,z7)',lower_dimensional_minimal_or_embedded_components_excluded=False,full_positive_Poisson_ring_CM=False,global_component_normality_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
 (ROOT/'reports/R4_NEW25_POSITIVE_AXIS_SMOOTH_REDUCED_COMPONENT_20261103_v63.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
