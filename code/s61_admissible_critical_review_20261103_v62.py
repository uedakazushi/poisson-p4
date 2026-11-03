#!/usr/bin/env python3
"""Direct original134 slice/Jacobian and critical-locus scope audit."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_quintic_open_review_20261103_v62 import invert_det
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

P=32003
def mod(v):
 q=s.Rational(v);assert int(q.q)%P;return int(q.p)*pow(int(q.q),-1,P)%P

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 paths=[Path(__file__),public,ROOT/'code/s61_quintic_open_review_20261103_v62.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
 data=read(ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json');record=next(v for v in read(ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json')['points']if v['component_id']=='T01');point=list(map(s.Rational,record['S61_coordinates_u1_through_u61']))
 first=[i for i,w in enumerate(data['parameter_weights'])if w==1];second=[i for i,w in enumerate(data['parameter_weights'])if w==2];fi={v:i for i,v in enumerate(first)};qi={v:i for i,v in enumerate(second)}
 chart=read(folder/'admissible-slice.json');chosen=[first.index(i-1)for i in chart['varying_old_u_one_based']];center=s.Matrix([point[i]for i in first]);assert list(map(str,center))==chart['center']
 z,u,beta=primitive_beta(S61());pairs=list(itertools.combinations(range(4),2));coords=[(ij,k)for ij in pairs for k in range(4)]
 basis=s.Matrix([[s.expand(beta[i,j]).coeff(z[k]).coeff(v)for v in u]for(i,j),k in coords]);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();B=beta.subs(dict(zip(u,center)),simultaneous=True)
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]);directions=[]
 for i in range(4):
  for j in range(i,4):
   S=s.zeros(4);S[i,j]=S[j,i]=1;A=J*S
   # Differentiate actual pullback g^t B(gz)g, using its polynomial chain rule.
   D=A.T*B+B*A+sum(((A*s.Matrix(z))[k]*B.diff(z[k])for k in range(4)),s.zeros(4))
   vector=s.Matrix([s.expand(D[i,j]).coeff(z[k])for(i,j),k in coords]);directions.append(left*vector[rows,:])
 dominance=s.Matrix.hstack(*directions,center).row_join(s.eye(16)[:,chosen]);assert dominance==s.Matrix([[s.Rational(v)for v in row]for row in chart['dominance_matrix']])
 assert dominance.det()==s.Rational(3564223982,1029)and mod(dominance.det())==9752
 R,*variables=ring(','.join([f'ss{i}'for i in range(1,6)]+[f'q{i}'for i in range(45)]),s.QQ);forms=[R(point[i])for i in range(61)]
 for j,col in enumerate(chosen):forms[first[col]]+=variables[j]
 for j,i in enumerate(second):forms[i]=variables[5+j]
 expected=[]
 for equation in data['equations']:
  value=R.zero
  for mon,n,den in equation:
   term=R(s.Rational(n,den))
   for i in mon:term*=forms[i]
   value+=term
  expected.append(value)
 encoded=read(folder/'restricted-original134.json.gz');actual=[]
 for row in encoded:
  value=R.zero
  for qmon,coefficient in row:
   part=R.zero
   for exponents,c in coefficient:
    term=R(s.Rational(c))
    for j,power in enumerate(exponents):
     if power:term*=variables[j]**power
    part+=term
   for i in qmon:part*=variables[5+i]
   value+=part
  actual.append(value)
 assert actual==expected and len(actual)==134
 critical=read(folder/'critical-determinant.json');C=np.array(critical['constant_row_combination'],dtype=np.int64);five=read(ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json');vs=s.symbols('v0:3');slope=s.Matrix([[s.sympify(p).coeff(v)for v in vs]for p in five['affine_plane_substitution']]);known=read(ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json')['known_five_mod_32003'];dets=[]
 for row in known:
  delta=[sum(mod(slope[i,j])*int(row[j+1])for j in range(3))%P for i in range(45)];values=[mod(v)for v in point]
  for i,v in zip(second,delta):values[i]=(values[i]+v)%P
  Jac=np.zeros((134,45),dtype=np.int64)
  for r,equation in enumerate(data['equations']):
   value=0
   for mon,n,den in equation:
    c=int(n)*pow(int(den),-1,P)%P;v=c
    for i in mon:v=v*values[i]%P
    value=(value+v)%P
    for pos,i in enumerate(mon):
     if i in qi:
      v=c
      for k,j in enumerate(mon):
       if k!=pos:v=v*values[j]%P
      Jac[r,qi[i]]=(int(Jac[r,qi[i]])+v)%P
   assert not value
  inverse,det=invert_det((C@Jac%P).tolist());dets.append(det)
 assert dets==critical['known_five_determinants_mod32003']==[18384,6387,6706,13504,28620]
 paths += [ROOT/'reports/S61_AFFINE_QUINTIC_OPEN_20261103_v62_REVIEW.json',ROOT/'research/notes/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62_REVIEW.md']
 report=dict(id='S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62_REVIEW',all_passed=True,
  fresh_producer_and_all_evidence_match=True,
  independent_primitive_polynomial_pullback_dominance_det='3564223982/1029',
  independent_QQ_50variable_restriction_of_ALL_original134=True,
  independent_original61coordinate_evaluation_and_raw_q45_Jacobians_at_knownfive=True,
  global_E_definition='det(C*Jac_q(original134)); degree_q<=45; no cyclic/frame denominators',
  all_positive_dimensional_fixed_beta_components_contained_in_E_zero=True,
  unit_original134_plus_E_sufficient_not_equivalent_to_dimension_zero=True,
  finite_model_frames_regular_on_nonempty_open_of_NEW_through_beta_slice=True,
  fixed_pair_center_finiteness_full_unit_or_wholefive_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS admissible5 dominance,all134 exact restriction and global raw-q critical determinant scope')

if __name__=='__main__':main()
