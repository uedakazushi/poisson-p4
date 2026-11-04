#!/usr/bin/env python3
"""ONE point-free generic rank certificate via a sparse2x2 block and7x3 Bareiss."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63'
class Stop(Exception):pass
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('57s_internal_boundary')));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';cp=base/'checkpoint.json';assert sha(cp)==read(base/'report.json')['evidence_sha256'][cp.name];c=read(cp);assert c['completed_pivots']==1
 review=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW/report.json';assert read(review)['all_passed']
 R,a,e,d=ring('a,e,d',s.QQ);dec=lambda q:R.from_dict({tuple(m):s.QQ(v)for m,v in q});enc=lambda q:[[list(m),str(v)]for m,v in q.items()]
 M=[[dec(q)for q in row]for row in c['matrix15x6']];Bcore=(50*a+18*e**2-9*e-2)*d+60*e-40;p=M[2][0];assert p and M[4][0]==p
 assert M[2][5]-M[4][5]==s.QQ(6,5)*Bcore*p;pairdet=M[2][0]*M[4][5]-M[2][5]*M[4][0];assert pairdet==-s.QQ(6,5)*Bcore*p**2 and pairdet
 brows=[7,9,10,11,12,13,14];bcols=[2,3,4];assert all(i in c['remaining_rows']for i in brows)
 actual=[i for i in c['remaining_rows']if any(M[i][j]for j in bcols)];assert actual==brows;assert all(not M[i][j]for i in[2,4]for j in bcols)
 B=[[M[i][j]for j in bcols]for i in brows];A=[row[:]for row in B];rows=list(range(7));cols=list(range(3));steps=[];previous=R.one;stop=None;complete=False;certificate={}
 def guard():
  if sum(len(q)for row in A for q in row)>150000:raise Stop('150000_TOTAL_stage_support')
  if time.monotonic()-start>50:raise Stop('50s_arithmetic_boundary')
 def mul(x,y):
  if not x or not y:return R.zero
  dx=max(sum(m)for m in x);dy=max(sum(m)for m in y);bound=min(len(x)*len(y),(dx+dy+1)*(dx+dy+2)*(dx+dy+3)//6)
  if bound>150000:raise Stop('150000_preproduct_support_bound')
  z=x*y
  if len(z)>150000:raise Stop('150000_product_terms')
  if time.monotonic()-start>50:raise Stop('50s_arithmetic_boundary')
  return z
 def det(matrix):
  n=len(matrix)
  if not n:return R.one
  return sum(((-1)**j*mul(matrix[0][j],det([r[:j]+r[j+1:]for r in matrix[1:]]))for j in range(n)),R.zero)
 def checkpoint():
  q={'completed_generic_pivots':len(steps),'steps':steps,'remaining_rows7':rows,'remaining_cols3':cols,'matrix7x3':[[enc(x)for x in row]for row in A],
   'previous_Bareiss_divisor':enc(previous),'current_TOTAL_support':sum(len(x)for row in A for x in row),'pivot_divisors_nonzero_QQ_polynomials_NOT_global_units':True,'seconds':time.monotonic()-start}
  (out/'checkpoint.json').write_text(json.dumps(q,indent=2)+'\n')
 checkpoint()
 try:
  while rows and cols:
   candidates=[(len(A[i][j]),i,j)for i in rows for j in cols if A[i][j]]
   if not candidates:break
   _,i,j=min(candidates);pivot=A[i][j];nextA=[row[:]for row in A]
   for u in rows:
    if u==i:continue
    for v in cols:
     if v==j:nextA[u][v]=R.zero;continue
     numerator=mul(pivot,A[u][v])-mul(A[u][j],A[i][v]);nextA[u][v]=numerator.exquo(previous);assert nextA[u][v]*previous==numerator
   if sum(len(x)for row in nextA for x in row)>150000:raise Stop('150000_TOTAL_next_stage_support')
   newrows=[s['pivot_row7']for s in steps]+[i];newcols=[s['pivot_col3']for s in steps]+[j];minor=det([[B[u][v]for v in newcols]for u in newrows]);assert minor==pivot and minor
   steps.append({'pivot_row7':i,'pivot_col3':j,'pivot':enc(pivot),'previous_exact_Bareiss_divisor':enc(previous),'original_minor_rows7':newrows,'original_minor_cols3':newcols,'original_minor_nonzero':enc(minor)});A=nextA;previous=pivot;rows.remove(i);cols.remove(j);guard();checkpoint()
  r=len(steps);I=[q['pivot_row7']for q in steps];J=[q['pivot_col3']for q in steps];P=[[B[i][j]for j in J]for i in I];Delta=det(P);assert Delta
  adj=[[(-1)**(i+j)*det([row[:i]+row[i+1:]for row in P[:j]+P[j+1:]])for j in range(r)]for i in range(r)]
  L=[]
  for row in B:
   lam=[sum((mul(row[J[k]],adj[k][j])for k in range(r)),R.zero)for j in range(r)]
   assert all(mul(Delta,row[k])==sum((mul(lam[j],B[I[j]][k])for j in range(r)),R.zero)for k in range(3));L.append(lam)
  syz=[]
  for free in range(3):
   if free in J:continue
   v=[R.zero]*3;v[free]=Delta
   for k in range(r):v[J[k]]=-sum((mul(adj[k][j],B[I[j]][free])for j in range(r)),R.zero)
   assert all(sum((mul(row[k],v[k])for k in range(3)),R.zero)==0 for row in B);syz.append(v)
  fullminor=mul(pairdet,Delta);assert fullminor
  certificate={'input14_remaining_rows':c['remaining_rows'],'input5_columns':c['remaining_columns'],'sparse_pair_rows':[2,4],'sparse_pair_columns':[0,5],
   'p':enc(p),'literal_pair_determinant':enc(pairdet),'pairdet_formula':'-(6/5)*B*p^2','all_other3_column_support_exact_rows':brows,'other3_original_columns':bcols,
   'actual_block7x3':[[enc(x)for x in row]for row in B],'block_generic_rank':r,'selected_rows7':I,'selected_cols3':J,'nonzero_selected_minor_Delta':enc(Delta),
   'cleared_upper_row_coordinates7xr':[[enc(x)for x in row]for row in L],'all21_cleared_upper_identities_verified':True,
   'polynomial_right_kernel_basis3':[[enc(x)for x in v]for v in syz],'all7_right_kernel_identities_verified':True,
   'nonzero_minor_of14x5_rows':[2,4]+[brows[i]for i in I],'nonzero_minor_of14x5_cols':[0,5]+[bcols[j]for j in J],'nonzero_minor_of14x5':enc(fullminor),
   'generic_rank14x5':2+r,'generic_rank_actual16_H7':4+r,'rank_decomposition_proof':'Two previously certified knownunit H directions + invertible 2x2 sparse pair overQQ(a,e,d) + the7x3 projected block. The pair rows vanish on the other3 columns, so those coefficients are unchanged by clearing the pair columns.',
   'unknown_p_and_Delta_inverted_ONLY_in_QQ_coefficient_field':True,'no_all_FULLU_rank_lower_or_new_globalunit_claim':True,'full_H7_uniform_route_impossible_if_generic_rank_below7':4+r<7,
   'no_newpoint_prime_GB_full239_or_newkernel_count_claim':True}
  complete=True
 except(Stop,MemoryError)as exc:stop=str(exc)or type(exc).__name__
 signal.alarm(0);(out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',cp,review]
 report={'id':ID,'complete_exact_generic_rank_certificate':complete,'block7x3_input_terms':sum(len(x)for row in B for x in row),'completed_generic_pivots':len(steps),'generic_actual_raw_H_rank':certificate.get('generic_rank_actual16_H7'),
  'full_H7_uniform_route_impossible':certificate.get('full_H7_uniform_route_impossible_if_generic_rank_below7'),'stop_reason':stop,'no_allpoint_rank_lower_or_global_unknown_inverse':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'checkpoint_sha256':sha(out/'checkpoint.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
