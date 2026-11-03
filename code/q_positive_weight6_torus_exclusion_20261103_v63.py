#!/usr/bin/env python3
"""A universal stabilizer-rank exclusion of T01,T02,T03 from the positive ray."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 cp=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
 catalog=json.loads(cp.read_text());evidence={}
 for name in ['T01','T02','T03']:
  record=next(r for r in catalog['records']if r['component_id']==name);diagram=record['source_diagram'];chars=[]
  for edge in diagram['edges']:
   char=[0]*5
   for i in edge['edge']:char[i]+=1
   for i,power in edge['orders'].items():char[int(i)]-=power
   assert sum(char)==0;chars.append(char)
  rank=s.Matrix(chars).rank()if chars else 0
  assert rank==diagram['smoothing_weight_rank']==(0 if name=='T01'else 1)
  torus=4-rank;assert torus==diagram['generic_torus_stabilizer_dimension']
  evidence[name]=dict(diagonal_characters=chars,character_rank=rank,projective_torus_stabilizer_dimension_at_least=torus,whole_component_projective_orbit_dimension_at_most=24-torus)
 orbit=ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63/report.json'
 od=json.loads(orbit.read_text());assert od['all_passed']and od['projective_PGL5_orbit_dimension']==22
 for name,digest in od['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
 assert all(v['whole_component_projective_orbit_dimension_at_most']<22 for v in evidence.values())
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),cp,orbit,ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63/certificate.json',ROOT/'research/notes/MPS_S61_CATALOGUE_20261103_v62.md',ROOT/'research/notes/Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63.md']
 result=dict(id='Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63',all_passed=True,
  whole_closed_components_excluded=['T01','T02','T03'],
  positive_nonzero_ray_projective_orbit_dimension=22,
  closed_action_rank_upper_bounds={k:v['whole_component_projective_orbit_dimension_at_most']for k,v in evidence.items()},
  proof='Dense fullMPS source has indicated diagonal torus; rank-upper-bound locus is Zariski closed, hence includes full component boundary',
  no_generic_dimension_or_sample_rank_substitution=True,
  other_MPS_membership_or_new_component_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print('POSITIVE RAY OUTSIDE WHOLE T01/T02/T03 BY UNIVERSAL ACTION RANK PASS')

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
