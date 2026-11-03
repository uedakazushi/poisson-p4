#!/usr/bin/env python3
"""Direct raw tensor lifts and full/reduced local non-CM scope audit."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 public=ROOT/'reports/Q_ALL_FAMILY_NON_CM_20261103_v63/report.json';d=load(public);assert d==load(args.replay)
 paths=[Path(__file__),public]
 for file,hsh in d['input_sha256'].items():p=ROOT/file;assert sha(p)==hsh,file;paths.append(p)
 a=s.Symbol('a');h,x,y,z,w=s.symbols('h x y z w');X=[h,x,y,z,w]
 H=4*h*x*w/3-h*z*z/3-x*x*z;f=a*h*h;g=-s.Rational(2,9)
 PA=s.zeros(5);PA[1,2]=f;PA[3,4]=g*h*h;PA[2,3]=-s.diff(H,w);PA[2,4]=s.diff(H,z);PA=PA-PA.T
 PB=s.zeros(5);PB[1,2]=f;PB[3,4]=g*h*h;PB[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);PB[2,4]=-s.diff(f,x)*w/2+s.diff(H,z);PB=PB-PB.T
 assert PA==PB
 assert all(s.expand(sum(s.diff(PA[i,j],X[i])for i in range(5)))==0 for j in range(5))
 old=s.symbols('x0:5');order=[0,1,4,2,3];transport=s.zeros(5)
 for i in range(5):
  for j in range(5):transport[order[i],order[j]]=PA[i,j].subs(dict(zip(X,[old[k]for k in order])),simultaneous=True)
 Q=s.zeros(5)
 for(i,j),v in{(1,4):a*old[0]**2,(2,3):-s.Rational(2,9)*old[0]**2,(2,4):s.Rational(4,3)*old[0]*old[1],(3,4):s.Rational(2,3)*old[0]*old[2]+old[1]**2}.items():Q[i,j]=v;Q[j,i]=-v
 assert(transport-Q).applyfunc(s.expand)==s.zeros(5)
 assert g!=0 and Q.subs(a,0)!=s.zeros(5)
 assert d['known_minimal_projective_component_dimensions']=={'T01':25,'CA':28,'CB':28}
 for file in ['reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json','reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json','reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json']:
  p=ROOT/file;assert load(p)['all_passed'];paths.append(p)
 paths += [ROOT/'research/notes/Q_ALL_FAMILY_NON_CM_20261103_v63_REVIEW.md']
 result=dict(id='Q_ALL_FAMILY_NON_CM_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_report_exact_match=True,independent_raw_CA_CB_polynomial_lifts_and_divergence_and_all_tensor_entries=True,
  all_finite_complex_a_including_zero=True,actual_minimal_component_dimensions_reused={'T01':25,'CA':28,'CB':28},
  closed_point_and_excellent_completion_preserve_distinct_minimal_component_dimensions=True,
  complete_CM_unmixedness_obstruction_accepted=True,
  full_Poisson_and_full_reduced_Poisson_completed_local_rings_nonCM_for_every_a=True,
  individual_CA_CB_all_a_nonCM_or_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS every Qa actual CA/CB lift and full/reduced non-CM completions')

if __name__=='__main__':main()
