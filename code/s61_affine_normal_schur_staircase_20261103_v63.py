#!/usr/bin/env python3
"""Companion staircase gives degree7 delta/vector,14square,30scalar135."""
import argparse
import gzip
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_affine_normal_global_quintic_20261103_v63 import ROOT,uadd,umul,sha
from fractions import Fraction

def det(M):
 D=DomainMatrix.from_Matrix(M);return s.expand(D.domain.to_sympy(D.det()))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 deps=[ROOT/'reports/S61_AFFINE_NORMAL_SCALAR_FITTING_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63/report.json']
 for dep in deps:
  old=json.loads(dep.read_text());assert old['all_passed']
  for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
  for rel,want in old.get('evidence_sha256',{}).items():assert sha(dep.parent/rel)==want,rel
 source=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz';data=json.loads(gzip.decompress(source.read_bytes()));decode=lambda poly:{d:{tuple(m):Fraction(c)for m,c in coeff}for d,coeff in poly};G=decode(data['monic_quintic_G']);R=list(map(decode,data['all134_remainders']))
 def reduce(poly):
  while poly and max(poly)>=5:
   k=max(poly);poly=uadd(poly,umul({k-5:poly[k]},G),Fraction(-1))
  return poly
 q_witness=None;l4_witness=None
 for i,r in enumerate(R):
  Q={m:c for m,c in r.get(0,{}).items()if len(m)==2};assert all(44 not in m for m in Q)
  if Q and q_witness is None:q_witness=dict(row=i,coefficient=str(next(iter(Q.values()))),monomial=list(next(iter(Q))))
  if r.get(4)and l4_witness is None:l4_witness=dict(row=i,coefficients=[[list(m),str(c)]for m,c in r[4].items()])
  assert all(len(m)<=1 for k,coeff in r.items()if k>=1 for m in coeff)
  for j in range(5):
   row=reduce({k+j:coeff for k,coeff in r.items()})
   for k,coeff in row.items():
    quad={m:c for m,c in coeff.items()if len(m)==2};assert all(len(m)<=2 for m in coeff)
    if k>j:assert not quad
    if k==j:assert quad==Q
 assert q_witness is not None and l4_witness is not None
 # Universalweightedcompanioncheck: Qhasweight2, allL/Hweight1,gcoeffs0.
 H,Q,t=s.symbols('H Q t');L=s.symbols('L0:5');gs=s.symbols('g0:5');T=s.zeros(5)
 for i in range(4):T[i+1,i]=1
 T[:,4]=s.Matrix([H-gs[0],-gs[1],-gs[2],-gs[3],-gs[4]]);r=s.Matrix([L[0]+Q,*L[1:]]);V=s.Matrix.hstack(*(T**j*r for j in range(4)));scaling={H:t*H,Q:t*t*Q,**{v:t*v for v in L}};scaled=V.subs(scaling,simultaneous=True);D=scaled[1:,:];delta=det(D);degree=s.Poly(delta,t).degree();top=s.expand(s.Poly(delta,t).coeff_monomial(t**7));assert degree==7 and s.expand(top+L[4]*Q**3)==0
 replacement=s.Matrix.vstack(-scaled[0,:],scaled[2:,:]);U=det(replacement);assert s.Poly(U,t).degree()<=7
 certificate=dict(all134_initial_coefficient_shapes_verified=True,all670_all3350_companion_staircase_entries_QQ_checked=True,
  universal_rule='quadraticparts(T^j*r)_k=0 fork>j;at k=j exactlyQ_i(L);initialr_h1..h4linear',
  delta_weighted_degree=7,universal_delta_top_degree='-L4*Q^3',source_nonzero_Q_witness=q_witness,source_nonzero_r4_witness=l4_witness,
  Cramer_U_weighted_degree_atmost=7,updated_bounds=dict(delta_vector=7,delta_square_scalar=14,delta_cubic_scalar=21,cleared134_scalar_generators=30,norm_generator=10),
  generic_degree7_exact=True,actual_source_delta_not_expanded=True,after_nonlinear_normal_basechange_degree_not_claimed=True,
  all_earlier_bounds8_16_34_remain_valid=True,no_actual_membership_or_unit=True)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),*deps,source,ROOT/'research/notes/S61_AFFINE_NORMAL_SCHUR_STAIRCASE_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_SCHUR_STAIRCASE_20261103_v63',all_passed=True,all3350_QQ_coefficient_pattern_verified=True,universal_companion_determinant_degree7_top_checked=True,vector7_square14_scalar135_degree30_bounds=True,memberships_and_units_uncomputed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('STAIRCASEDELTA7/VECTOR7/SQUARE14/SCALAR30 PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
