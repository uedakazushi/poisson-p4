#!/usr/bin/env python3
"""Independently normalize the degree42 determinant from its QQ linear entries."""
import argparse,gzip,hashlib,json,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];P=32003

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/S61_AFFINE_IDEMPOTENT42_20261103_v63';report=json.loads((folder/'report.json').read_text());assert report['all_passed']
 for name,want in report['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want
 cert=json.loads((folder/'idempotent42-certificate.json').read_text())
 old=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63'
 qq=json.loads(gzip.decompress((old/'linear129x42-QQ.json.gz').read_bytes()))
 direction=json.loads((old/'residual42-certificate.json').read_text())
 def residue(v):
  v=Fraction(v);assert v.denominator%P
  return v.numerator*pow(v.denominator,-1,P)%P
 factors=[[[ (tuple(m),residue(c))for m,c in entry]for entry in row]for row in qq['all129_linear_factor_rows']]
 def evaluate(point,constant=1):return [[sum(c*(point[m[0]]if m else constant)for m,c in entry)%P for entry in row]for row in factors]
 def det(a):
  a=[r[:]for r in a];d=1
  for j in range(len(a)):
   k=next((k for k in range(j,len(a))if a[k][j]%P),None)
   if k is None:return 0
   if k!=j:a[j],a[k]=a[k],a[j];d=-d
   pivot=a[j][j]%P;d=d*pivot%P;inv=pow(pivot,-1,P)
   for i in range(j+1,len(a)):
    ratio=a[i][j]*inv%P
    for k in range(j+1,len(a)):a[i][k]=(a[i][k]-ratio*a[j][k])%P
  return d%P
 knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json';known=json.loads(knownpath.read_text())['known_five_mod_32003']
 selected=cert['selected_normal_rows'];slot=cert['replaced_slot'];replacements=cert['replacement_rows'];C=[];matrices=[]
 for point in known:
  A=evaluate(point[1:]+[0]*42);matrices.append(A);row=[]
  for k in replacements:
   square=[A[j][:]for j in selected];square[slot]=A[k][:];row.append(det(square))
  C.append(row)
 assert C==cert['coefficient_matrix_at_knownfive']and det(C)==91
 aug=[row[:]+[1]for row in C]
 for j in range(5):
  k=next(k for k in range(j,5)if aug[k][j]);aug[j],aug[k]=aug[k],aug[j]
  inv=pow(aug[j][j],-1,P);aug[j]=[a*inv%P for a in aug[j]]
  for k in range(5):
   if k!=j:
    c=aug[k][j];aug[k]=[(a-c*b)%P for a,b in zip(aug[k],aug[j])]
 weights=[row[-1]for row in aug];assert weights==cert['replacement_weights_at_frozen_prime']
 def normalized(A):
  square=[A[j][:]for j in selected];square[slot]=[sum(w*A[r][j]for r,w in zip(replacements,weights))%P for j in range(42)]
  return det(square)
 assert [normalized(A)for A in matrices]==[1]*5
 point=direction['leading_degree42_direction_z'][1:]+direction['leading_degree42_direction_n']
 leading=normalized(evaluate(point,constant=0));assert leading==17048
 paths=[Path(__file__),folder/'report.json',folder/'idempotent42-certificate.json',old/'linear129x42-QQ.json.gz',old/'residual42-certificate.json',knownpath,ROOT/'research/notes/S61_AFFINE_IDEMPOTENT42_20261103_v63_REVIEW.md']
 result=dict(id='S61_AFFINE_IDEMPOTENT42_20261103_v63_REVIEW',all_passed=True,
  independent_source='Frozen exact QQ129x42 coefficients; no original projection engine or producer determinant imports',
  all25_replacement_determinants_compared=True,normalization_5x5_determinant=91,
  independently_solved_normalization_weights=weights,all_known_five_determinants_one=True,
  exact_degree42_leading_coefficient_mod32003=leading,
  full_generic_determinant_equals_global_idempotent_scope_reviewed=True,
  generic_beta_weights_are_finite_recipe_not_specialized_constants=True,
  actual_residual_unit_or_exhaustion_proved=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT DEGREE42 IDEMPOTENT REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
