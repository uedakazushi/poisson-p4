#!/usr/bin/env python3
"""Independent bounded checks for the two exact global annihilators."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def trim(a):
 a=[int(x)%P for x in a]
 while a and a[-1]==0:a.pop()
 return a


def resultant(f,g):
 """Plain coefficient-list Euclidean resultant; no CAS/native engine."""
 f,g=trim(f),trim(g);answer=1
 while g:
  n,m=len(f)-1,len(g)-1
  if m==0:return answer*pow(g[0],n,P)%P
  if n<m:
   if n*m%2:answer=-answer%P
   f,g=g,f;continue
  remainder=f.copy();inverse=pow(g[-1],-1,P)
  for j in range(n-m,-1,-1):
   c=remainder[m+j]*inverse%P
   if c:
    for k,v in enumerate(g):remainder[k+j]=(remainder[k+j]-c*v)%P
  remainder=trim(remainder)
  if not remainder:return 0
  if n*m%2:answer=-answer%P
  answer=answer*pow(g[-1],n-(len(remainder)-1),P)%P
  f,g=g,remainder
 return 0


def rank(matrix):
 a=np.asarray(matrix,dtype=np.int64).copy()%P;r=0
 for col in range(a.shape[1]):
  choices=np.flatnonzero(a[r:,col])
  if not len(choices):continue
  pivot=r+int(choices[0]);a[[r,pivot]]=a[[pivot,r]]
  a[r]=a[r]*pow(int(a[r,col]),-1,P)%P
  a[r+1:]=(a[r+1:]-a[r+1:,col,None]*a[r,None,:])%P
  r+=1
  if r==a.shape[0]:break
 return r


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--resultant-replay',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 rdir=ROOT/'reports/S61_AFFINE_GLOBAL_RESULTANT_ANNIHILATOR_20261103_v62'
 gdir=ROOT/'reports/S61_AFFINE_T_DETERMINANT_ANNIHILATOR_20261103_v62'
 paths=[rdir/'report.json',gdir/'report.json'];audited=[]
 for path in paths:
  d=json.loads(path.read_text());assert d['all_passed']
  for name,expected in d['input_sha256'].items():
   p=ROOT/name;assert sha(p)==expected,name;audited.append(p)
  for name,expected in d['evidence_sha256'].items():
   p=path.parent/name;assert sha(p)==expected,str(p);audited.append(p)
 # Native/resultant replay was fresh; binary logs/times may differ.
 replay=args.resultant_replay
 assert (rdir/'native.txt').read_bytes()==(replay/'native.txt').read_bytes()
 assert (rdir/'annihilator-circuit.json').read_bytes()==(replay/'annihilator-circuit.json').read_bytes()
 inp=list(map(int,(rdir/'polynomials.in').read_text().split()))
 assert inp[0]==1031
 chi,pt=inp[1:1033],inp[1033:]
 assert len(chi)==1032 and len(pt)==1031
 second=[2*c%P for c in pt];second[0]=(second[0]-1)%P
 independent=resultant(chi,second);assert independent==30737
 circuit=json.loads((rdir/'annihilator-circuit.json').read_text())
 u=circuit['independent_lambda2_bezout_coefficients']['chi']
 w=circuit['independent_lambda2_bezout_coefficients']['two_pt_minus_one']
 # Integer coefficient products, independently of FLINT and NumPy convolution.
 product=[0]*max(len(u)+len(chi)-1,len(w)+len(second)-1)
 for left,right in [(u,chi),(w,second)]:
  for i,a in enumerate(left):
   for j,b in enumerate(right):product[i+j]+=a*b
 assert product[0]%P==1 and all(c%P==0 for c in product[1:])
 with np.load(gdir/'central-action-blocks.npz')as archive:
  L=archive['linear_to_quadratic'].astype(np.int64)
  C=archive['quadratic_to_cubic'].astype(np.int64)
  CL=archive['linear_to_cubic'].astype(np.int64)
 assert L.shape==(301,29) and C.shape==(700,301)
 assert np.array_equal(C@L%P,CL)
 ranks=[1031,1+rank(L)+rank(C),int(np.any(L[:,0]))+rank(CL),int(np.any(CL[:,0])),0]
 assert ranks==[1031,330,30,1,0]
 delta=sum(range(1031))-(29+2*301+3*700)
 value=pow(pow(2,delta,P)*17765%P,2061,P)*independent%P
 assert delta==528234 and value==27775
 paths+=audited+[Path(__file__),
   ROOT/'research/notes/S61_AFFINE_GLOBAL_ANNIHILATORS_20261103_v62.md',
   ROOT/'research/notes/S61_GLOBAL_ANNIHILATORS_20261103_v62_REVIEW.md',
   ROOT/'reports/S61_GLOBAL_ACTION_CYCLIC_OPEN_20261103_v62_REVIEW.json']
 report=dict(id='S61_GLOBAL_ANNIHILATORS_20261103_v62_REVIEW',all_passed=True,
             strict_original_ideal_membership_by_adjugate_accepted=True,
             G_degree_upper=1031,G_constant_QQ=-1,cyclic_localization_needed=False,
             independent_saved_central_block_rank_sequence=ranks,
             central_blocks_independently_reconstructed_from_original_quadrics=False,
             fresh_native_full1032_radial_resultant_replay_matches=True,
             independent_plain_Python_Euclidean_resultant_at_two=independent,
             independent_plain_integer_Bezout_coefficient_identity=True,
             independent_radial_scaling_and_global_R_at_two_v=value,
             generic_spread_uses_corrected_full_quartic_inverse=True,
             generic_huge_coefficients_serialized=False,
             q_equals_radial_G_assumed=False,annihilator_zero_locus_exhausted=False,
             affine_surface_exclusion_claimed=False,full_affine_exhaustion_claimed=False,
             input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(paths)})
 with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
 print('PASS exact global ideal-membership scope, central ranks and independent resultant/Bezout')


if __name__=='__main__':main()
