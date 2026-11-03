#!/usr/bin/env python3
"""Independent original-S61 evaluation and critical rank-five certificate audit."""
import argparse, hashlib, json, time
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 stem='S61_AFFINE_CRITICAL_QUINTIC_PRESENTATION_20261103_v63'
 folder=ROOT/'reports'/stem
 report=json.loads((folder/'report.json').read_text())
 cert=json.loads((folder/'center-certificate.json').read_text());p=cert['prime']
 assert p==32003 and report['all_passed']
 checked=0
 for name,want in report['input_sha256'].items():
  assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,name
  checked+=1
 rawpath=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
 raw=json.loads(rawpath.read_text())
 pointpath=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
 base=next(v for v in json.loads(pointpath.read_text())['points'] if v['component_id']=='T01')['S61_coordinates_u1_through_u61']
 def residue(v):
  v=Fraction(v);assert v.denominator%p
  return v.numerator*pow(v.denominator,-1,p)%p
 base=list(map(residue,base));qi=[i for i,w in enumerate(raw['parameter_weights'])if w==2]
 critpath=ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62/critical-determinant.json'
 critical=json.loads(critpath.read_text());assert critical['original_q_indices_zero_based']==qi
 C=critical['constant_row_combination']
 def det(a):
  a=[[x%p for x in r]for r in a];d=1
  for j in range(len(a)):
   k=next((k for k in range(j,len(a))if a[k][j]),None)
   if k is None:return 0
   if j!=k:a[j],a[k]=a[k],a[j];d=-d
   pivot=a[j][j];d=d*pivot%p;inv=pow(pivot,-1,p)
   for i in range(j+1,len(a)):
    ratio=a[i][j]*inv%p
    for k in range(j+1,len(a)):a[i][k]=(a[i][k]-ratio*a[j][k])%p
  return d%p
 values=[];rootchecks=0
 for centered in cert['original_q_points']:
  point=base[:]
  for j,i in enumerate(qi):point[i]=(point[i]+centered[j])%p
  J=[]
  for equation in raw['equations']:
   value=0;gradient=[0]*61
   for mon,n,d in equation:
    co=n*pow(d,-1,p)%p;term=co
    for i in mon:term=term*point[i]%p
    value=(value+term)%p
    for pos,i in enumerate(mon):
     derivative=co
     for pos2,k in enumerate(mon):
      if pos2!=pos:derivative=derivative*point[k]%p
     gradient[i]=(gradient[i]+derivative)%p
   assert value==0
   rootchecks+=1;J.append([gradient[i]for i in qi])
  M=[[sum(c*J[k][j]for k,c in enumerate(row))%p for j in range(45)]for row in C]
  values.append(det(M))
 assert values==cert['original_E_values']==critical['known_five_determinants_mod32003']
 def trim(f):
  f=[v%p for v in f]
  while f and not f[-1]:f.pop()
  return f
 def add(f,g):return trim([(f[i]if i<len(f)else 0)+(g[i]if i<len(g)else 0)for i in range(max(len(f),len(g)))])
 def mul(f,g):
  h=[0]*max(0,len(f)+len(g)-1)
  for i,a in enumerate(f):
   for j,b in enumerate(g):h[i+j]=(h[i+j]+a*b)%p
  return trim(h)
 def rem(f,g):
  f=trim(f)
  while len(f)>=len(g):
   shift=len(f)-len(g);c=f[-1]*pow(g[-1],-1,p)%p
   for i,b in enumerate(g):f[shift+i]=(f[shift+i]-c*b)%p
   f=trim(f)
  return f
 def ev(f,x):
  value=0
  for a in reversed(f):value=(value*x+a)%p
  return value
 g=cert['g5'];e=cert['critical_remainder_e5'];roots=cert['primitive_roots']
 product=[1]
 for h in roots:product=mul(product,[-h,1])
 assert product==g and len(set(roots))==5 and g[-1]==1
 assert [ev(e,h)for h in roots]==values
 assert add(mul(e,cert['Bezout_e']),mul(g,cert['Bezout_g']))==[1]
 cols=[]
 for j in range(5):
  f=rem([0]*j+e,g);cols.append(f+[0]*(5-len(f)))
 M=[list(row)for row in zip(*cols)];assert M==cert['multiplication_by_e']
 assert det(M)==cert['norm']==9485
 inverse=cert['multiplication_inverse']
 assert [[sum(M[i][k]*inverse[k][j]for k in range(5))%p for j in range(5)]for i in range(5)]==[[int(i==j)for j in range(5)]for i in range(5)]
 recipe=json.loads((folder/'finite-recipe.json').read_text())
 assert recipe['exact_critical_module_presentation_shape']==[5,680]
 assert len(recipe['residuals'])==135 and recipe['etale_discriminant_not_needed_for_module_presentation']
 inputs=[Path(__file__),rawpath,pointpath,critpath,folder/'report.json',folder/'finite-recipe.json',folder/'center-certificate.json',ROOT/'research/notes'/f'{stem}.md']
 result=dict(id=stem+'_REVIEW',all_passed=True,producer_input_hash_checks=checked,
  original_raw_S61_equation_evaluations=rootchecks,original_q_jacobian_determinants=values,
  exact_center_Bezout_and_full5x5_inverse=True,center_norm=9485,
  full_critical5x680_presentation_scope_reviewed=True,
  general_coefficients_are_finite_inverse_recipe_not_expanded=True,
  no_quintic_discriminant_required=True,cyclic_highblock_tnorm_complements_resolved=False,
  generic_affine_exhaustion_proved=False,
  input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in inputs},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print('INDEPENDENT CRITICAL QUINTIC REVIEW PASS',result['seconds'])

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);main(parser.parse_args().output)
