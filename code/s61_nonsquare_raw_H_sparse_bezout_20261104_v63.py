#!/usr/bin/env python3
"""ONE fixed4-entry scalar-column constant/affine QQ span; no Gaussian/GB/point."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_SPARSE_BEZOUT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';rp=base/'report.json';p=base/'checkpoint.json';r=read(rp);assert sha(p)==r['evidence_sha256'][p.name];cp=read(p);pre=read(base/'preflight.json');assert cp['completed_pivots']==1
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();dec=lambda ts:R.from_dict({tuple(m):s.QQ(c)for m,c in ts});enc=lambda q:[[list(m),str(c)]for m,c in q.items()];encK=lambda q:[enc(q.numer),enc(q.denom)]
 F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4;known=[a,d,F,B,C,J]
 chosen=[2,4,5,8];assert all(i in cp['remaining_rows']for i in chosen)and cp['remaining_columns']==[0,2,3,4,5]
 matrix=[[dec(q)for q in row]for row in cp['matrix15x6']];actual=[[matrix[i][j]for j in[0,5]]for i in chosen]
 assert all(not matrix[i][j]for i in chosen for j in[1,2,3,4]);assert sum(len(q)for row in actual for q in row)<10000
 primitive=[];divisors=[];unit_records=[]
 for row in actual:
  q=row[0];assert q;c=q.LC;q=q/c;ex=[]
  for f in known:
   n=0
   while not q.rem(f):q=q.exquo(f);n+=1
   ex.append(n)
  u=R.ground_new(c)
  for f,n in zip(known,ex):u*=f**n
  assert q*u==row[0];primitive.append(q);divisors.append(u);unit_records.append({'constant':str(c),'a_d_F_B_C_J_powers':ex,'literal_divisor':enc(u)})
 targets=[R.one,a,d,F,B,C,J,a*d];names=['1','a','d','F','B','C','J','ad'];attempts=[];chosen_solution=None
 for degree,affine in[(0,[R.one]),(1,[R.one,a,e,d])]:
  basis=[q*m for q in primitive for m in affine];mons=sorted(set().union(*(set(q)for q in basis+targets)),reverse=True);n=len(basis)
  data=[[s.Rational(str(q.get(mon,0)))for q in basis+targets]for mon in mons]
  dm=DomainMatrix.from_Matrix(s.Matrix(data)).convert_to(s.QQ);rr,piv=dm.rref();M=rr.to_Matrix();rank=sum(j<n for j in piv);solutions=[]
  for ti,T in enumerate(targets):
   bad=any(all(M[i,j]==0 for j in range(n))and M[i,n+ti]!=0 for i in range(M.rows))
   if bad:continue
   v=[s.QQ.zero]*n
   for i,j in enumerate(piv):
    if j<n:v[j]=s.QQ(M[i,n+ti])
   cof=[sum((R.ground_new(v[i*len(affine)+k])*mon for k,mon in enumerate(affine)),R.zero)for i in range(4)]
   assert sum((u*q for u,q in zip(cof,primitive)),R.zero)==T
   solutions.append({'target':names[ti],'target_polynomial':enc(T),'primitive_cofactors':[enc(q)for q in cof]})
  attempts.append({'cofactor_degree':degree,'monomial_basis':['1']if degree==0 else['1','a','e','d'],'coefficient_unknowns':n,'coefficient_equations':len(mons),'coefficient_rank':rank,'feasible_target_solutions':solutions})
  if solutions:chosen_solution=min(solutions,key=lambda z:sum(len(q)for q in z['primitive_cofactors']));break
 result={'source_matrix_rows':chosen,'source_target_column':0,'retained_tail_column':5,'variable_order':['a','e','d'],'exact_target_class':names,'raw4x2':[[enc(q)for q in row]for row in actual],
  'known_unit_entry_divisors':unit_records,'primitive_scalar_column':[enc(q)for q in primitive],'constant_then_affine_attempts':attempts,'success':chosen_solution is not None,
  'no_GB_point_prime_minorsearch_Gaussian_or_unknown_inverse':True,'no_full_H7_full57_full239_or_source_basis_count_claim':True}
 if chosen_solution:
  cofs=[dec(q)for q in chosen_solution['primitive_cofactors']];T=dec(chosen_solution['target_polynomial']);source_cofs=[K(q)/K(u)for q,u in zip(cofs,divisors)]
  combined=[sum((c*K(row[j])for c,row in zip(source_cofs,actual)),K.zero)for j in range(2)];assert combined[0]==K(T)
  # Transport the actual full-source combination back to the original genuine16.
  delta=dec(pre['source_delta_change']['delta']);eta=dec(pre['pair_values']['eta']);clear=[dec(q)for q in pre['other15_full_source_clearings_numerators_over_eta']];genuine=pre['remaining15_genuine_indices'];Cdelta=[K.zero]*16;Cdelta[6]=K(delta);Cdelta[1]=-K(delta)
  S=[]
  for idx,g in enumerate(genuine):
   row=[-(K(clear[idx])/K(eta))*q for q in Cdelta];row[g]+=K.one;S.append(row)
  def unit(q):
   z=K.from_expr(s.sympify(q['constant']))
   for f,n in zip(known,q['a_d_F_B_C_J_powers']):z*=K(f)**n
   return z
  for i,u in enumerate(cp['initial_row_divisors']):S[i]=[q/unit(u)for q in S[i]]
  for step in cp['steps']:
   row=step['pivot_row'];pivot=unit(step['pivot_unit']);old=[q for q in S[row]]
   for elim in step['eliminations']:
    i=elim['row'];c=K(dec(elim['coefficient_before']));div=unit(elim['row_divisor_after']);S[i]=[(pivot*x-c*y)/div for x,y in zip(S[i],old)]
  full=[sum((c*S[i][j]for c,i in zip(source_cofs,chosen)),K.zero)for j in range(16)]
  original=[[K(dec(q))for q in row]for row in pre['raw16x7']];oldH=[sum((full[i]*original[i][j]for i in range(16)),K.zero)for j in range(7)];H=pre['constant_coefficient_column_frame7x7'];X=[sum((H[j][k]*oldH[k]for k in range(7)),K.zero)for j in range(7)]
  cZ=K(dec(pre['pair_values']['cZ']));cW=K(dec(pre['pair_values']['cW']));zz=K(dec(pre['pair_values']['delta_z']));ww=K(dec(pre['pair_values']['delta_w']));gamma=K(dec(pre['pair_values']['gamma']));frame=X[:5]+[(cZ*X[5]+cW*X[6])/K(eta),(-ww*X[5]+zz*X[6])/gamma]
  expected=[K(T),K.zero,K.zero,K.zero,K.zero,K.zero,combined[1]];assert frame==expected
  # Every denominator in the saved source/tail values must remain a named FULLU unit.
  for value in source_cofs+full+combined:
   q=value.denom
   for f in known:
    while not q.rem(f):q=q.exquo(f)
   assert q.is_ground and q
  result.update({'selected_solution':chosen_solution,'actual4_source_coefficients':[encK(q)for q in source_cofs],'literal_full4x2_row_combination':[encK(q)for q in combined],
   'all7_actual_transformed_H_coefficients_verified':[encK(q)for q in frame],'full_original_genuine16_source_coefficients':[encK(q)for q in full],
   'full_source_pure_parts':'Use these16 actual scalar coefficients on the original genuine cubic polynomials; do not discard their pure coefficients. No pure84 expansion here.',
   'target_tail_may_be_nonzero':bool(combined[1]),'no_source_completion_or_target_deletion_performed':True})
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,p,base/'preflight.json']
 report={'id':ID,'fixed4x1_exact_span_success':chosen_solution is not None,'selected_target':None if chosen_solution is None else chosen_solution['target'],'last_cofactor_degree':attempts[-1]['cofactor_degree'],'literal_full4x2_and_genuine16_transport_verified':chosen_solution is not None,
  'no_unknown_inverse_GB_point_Gaussian_or_automatic_source_completion':True,'input_sha256':{str(q.relative_to(ROOT)):sha(q)for q in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(q.stat().st_size for q in out.iterdir()if q.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
