#!/usr/bin/env python3
"""Independent c-monomial and finite24-coefficient cocharacter audit."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
 r=read(p);assert r['all_passed']
 for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
 for n,h in r.get('evidence_sha256',{}).items():
  q=ROOT/n if (ROOT/n).exists() else p.parent/n
  assert sha(q)==h,n
 return r
def matmul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]

def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 assert not out.exists();out.mkdir(parents=True)
 base=ROOT/'reports/Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63';ap=base/'report.json';accepted(ap);cp=base/'certificate.json';c=read(cp)
 src=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json';d=read(src)
 inv=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';rep=read(inv)
 six=ROOT/'reports/Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW/report.json';accepted(six)
 const=ROOT/'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json';accepted(const)
 lowerp=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63_REVIEW/report.json';lower=accepted(lowerp)
 assert lower['NEW25_firstjet_defining_equation_degree_lower_bound']==26
 bir=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW/report.json';accepted(bir)
 R,*_=ring('u,v,c,X,W,B',s.QQ);F=R.to_field()
 E=[-1,-1,0,0,1,0,1,0,1,2,1,2,0,1,0,1]
 assert E==d['exact_c_exponents_beta16']==c['source_c_scaling_exponents']
 records=[]
 for i,text in enumerate(d['primitive_beta16']):
  f=F.from_expr(s.sympify(text))
  if not f:records.append(dict(index=i,identically_zero=True,exponent_fixed_by_representation=E[i]));continue
  n={m[2] for m in f.numer};de={m[2] for m in f.denom}
  assert len(n)==len(de)==1
  ne=next(iter(n));den=next(iter(de));assert ne-den==E[i]
  records.append(dict(index=i,identically_zero=False,numerator_c_degree=ne,denominator_c_degree=den,c_power=E[i]))
 B=[[Q(v) for v in row] for row in rep['primitive24_by16_basis']]
 pairs=list(itertools.combinations(range(4),2));coords=[(i,j,k) for i,j in pairs for k in range(4)]
 h=[-1,1,-1,1];weights=[2*e-1 for e in E]
 for row,(i,j,k) in enumerate(coords):
  actual=h[i]+h[j]+h[k]
  for col in range(16):assert actual*B[row][col]==weights[col]*B[row][col]
 for col,w in enumerate(rep['Cartan_weights']):assert -sum(w)==weights[col]
 J=[[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]
 assert all(not J[i][j] or h[i]+h[j]==0 for i in range(4) for j in range(4))
 generators=[]
 for i in range(4):
  for j in range(i,4):
   sym=[[0]*4 for k in range(4)];sym[i][j]=sym[j][i]=1;generators.append(matmul(J,sym))
 assert [[-generators[1][i][j]-generators[8][i][j] for j in range(4)] for i in range(4)]==[[h[i] if i==j else 0 for j in range(4)] for i in range(4)]
 rho=[[[Q(v) for v in row] for row in M] for M in rep['all_ten_representation_matrices']]
 actual=[[-rho[1][i][j]-rho[8][i][j] for j in range(16)] for i in range(16)]
 assert actual==[[weights[i] if i==j else 0 for j in range(16)] for i in range(16)]
 sectionp=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json';section=read(sectionp)
 assert s.Symbol('c') not in s.sympify(section['section_v']).free_symbols
 assert not c['ZN_primitivity_or_target_index_one_proved']
 charp=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/independent-certificate.json'
 chars=read(charp)['characters_through24'];assert next(x['invariant_dimension'] for x in chars if x['degree']==14)==3
 certificate=dict(all16_c_checks=records,standard_cocharacter_exponents=h,primitive24_coordinate_order=coords,
  all384_finite_torus_basis_coefficients_exactly_checked=True,cocharacter_weights16=weights,
  independent_standard_h_equals_minus_generator1_minus_generator8=True,original_rho_h_is_exact_diagonal=True,
  c_independent_section=section['section_v'],general_invariant_character='degree2n -> c^n',
  conditional_shape_field='C(r,s,q,B)',conditional_index_one='If full shape field is in actual L, nonzeroI14/F4^3 recoversxi, henceL=K',
  shape_generation_ZN_primitivity_or_index_one_proved=False)
 ep=out/'independent-certificate.json';ep.write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,cp,src,inv,six,const,lowerp,bir,sectionp,charp,
        ROOT/'code/q_new_family_general_source_scale_20261104_v63.py',ROOT/'research/notes/Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63.md']
 report=dict(id=STEM,all_passed=True,independent_rational_numerator_denominator_c_monomial_checks16=True,
  independent_finite_actual24_tensor_cocharacter_not_only_stored_Lie_weight_metadata=True,
  original_rho_h_by_minus_generators1_and8_full256_entries_checked=True,
  all_degrees_invariant_source_homogeneity_c_power_D_over2_proved=True,c_independent_section_preserves_covariance=True,
  degree14_form_xi7_times_shape_and_nonvanishing_by_lower26_accepted=True,
  target_index_one_criterion_is_conditional_not_premises_proved=True,
  actual_shape_generation_ZN_primitivity_target_index_one_or_NEW_equation_claim=False,
  no_invariant_source26_expansion_point_GB_or_author_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k] for k in ['all_passed','all_degrees_invariant_source_homogeneity_c_power_D_over2_proved','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
