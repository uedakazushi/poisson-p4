#!/usr/bin/env python3
"""ONE fixed row2-row4 difference, named-unit stripping only; no search/GB/Gaussian."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_ROW_DIFFERENCE_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';cp=base/'checkpoint.json';assert sha(cp)==read(base/'report.json')['evidence_sha256'][cp.name];x=read(cp)
 sp=ROOT/'reports/S61_NONSQUARE_RAW_H_SPARSE_BEZOUT_20261104_v63/certificate.json';sr=read(sp.parent/'report.json');assert sha(sp)==sr['certificate_sha256'];sparse=read(sp)
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();dec=lambda z:R.from_dict({tuple(m):s.QQ(c)for m,c in z});enc=lambda q:[[list(m),str(c)]for m,c in q.items()];encK=lambda q:[enc(q.numer),enc(q.denom)]
 F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4;known=[a,d,F,B,C,J]
 assert x['remaining_columns']==[0,2,3,4,5]and 2 in x['remaining_rows']and 4 in x['remaining_rows'];rows=[[dec(z)for z in x['matrix15x6'][i]]for i in[2,4]];diff=[p-q for p,q in zip(*rows)];assert all(not diff[j]for j in range(5));tail=diff[5];assert tail and len(tail)<=129
 c=tail.LC;primitive=tail/c;powers=[]
 for f in known:
  n=0
  while not primitive.rem(f):primitive=primitive.exquo(f);n+=1
  powers.append(n)
 factor=R.ground_new(c)
 for f,n in zip(known,powers):factor*=f**n
 assert factor*primitive==tail;isunit=primitive.is_ground
 # Previously computed rank3 proves uniqueness of this constant column0 nullspace direction.
 attempt=sparse['constant_then_affine_attempts'][0];assert attempt['cofactor_degree']==0 and attempt['coefficient_unknowns']==4 and attempt['coefficient_rank']==3
 assert sparse['primitive_scalar_column'][0]==sparse['primitive_scalar_column'][1]
 result={'coefficient_domain':'FULLU','fixed_source_rows':[2,4],'difference_coefficients':[1,-1],'remaining_target_columns':x['remaining_columns'],
  'all_six_native_target_coefficients_difference':[enc(q)for q in diff],'all_but_column5_literal_zero':True,'tail':enc(tail),'tail_terms':len(tail),
  'exact_named_factorization':{'QQ_constant':str(c),'a_d_F_B_C_J_powers':powers,'remaining_polynomial':enc(primitive)},'literal_factor_reconstruction_verified':True,
  'tail_is_known_named_unit':bool(isunit),'only_constant_column0_nullspace_direction_up_to_scale':True,'no_geometric_rankdefect_or_nonunit_in_FULLU_claim_if_residual_nonconstant':True,
  'no_further_Gaussian_targetcompletion_unknownfactor_inversion_or_source_growth':True}
 if isunit:
  pre=read(base/'preflight.json');delta=dec(pre['source_delta_change']['delta']);eta=dec(pre['pair_values']['eta']);clear=[dec(z)for z in pre['other15_full_source_clearings_numerators_over_eta']];genuine=pre['remaining15_genuine_indices'];cd=[K.zero]*16;cd[6]=K(delta);cd[1]=-K(delta);S=[]
  for i,g in enumerate(genuine):
   v=[-K(clear[i])/K(eta)*q for q in cd];v[g]+=K.one;S.append(v)
  def u(z):
   q=K.from_expr(s.sympify(z['constant']))
   for f,n in zip(known,z['a_d_F_B_C_J_powers']):q*=K(f)**n
   return q
  for i,z in enumerate(x['initial_row_divisors']):S[i]=[q/u(z)for q in S[i]]
  for step in x['steps']:
   pr=step['pivot_row'];p=u(step['pivot_unit']);saved=S[pr][:]
   for el in step['eliminations']:
    i=el['row'];z=K(dec(el['coefficient_before']));h=u(el['row_divisor_after']);S[i]=[(p*q-z*t)/h for q,t in zip(S[i],saved)]
  full=[p-q for p,q in zip(S[2],S[4])];original=[[K(dec(q))for q in row]for row in pre['raw16x7']];old=[sum((full[i]*original[i][j]for i in range(16)),K.zero)for j in range(7)];H=pre['constant_coefficient_column_frame7x7'];X=[sum((H[j][k]*old[k]for k in range(7)),K.zero)for j in range(7)]
  zz=K(dec(pre['pair_values']['delta_z']));ww=K(dec(pre['pair_values']['delta_w']));cZ=K(dec(pre['pair_values']['cZ']));cW=K(dec(pre['pair_values']['cW']));gamma=K(dec(pre['pair_values']['gamma']));new=X[:5]+[(cZ*X[5]+cW*X[6])/K(eta),(-ww*X[5]+zz*X[6])/gamma];assert new==[K.zero]*6+[K(tail)]
  result.update({'full_original_genuine16_source_coefficients':[encK(q)for q in full],'all7_transformed_actualH_coefficients_verified':[encK(q)for q in new],'pure_source_parts_retained_lazy_with_saved16_coefficients':True})
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),cp,base/'report.json',sp,sp.parent/'report.json']
 report={'id':ID,'fixed_difference_all_other_coordinates_zero':True,'tail_is_known_named_unit':bool(isunit),'tail_terms':len(tail),'named_unit_powers':powers,'remaining_polynomial_terms':len(primitive),'no_general_ranknegative_GB_point_or_Gaussian_claim':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
