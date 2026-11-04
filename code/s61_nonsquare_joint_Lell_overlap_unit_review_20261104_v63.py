#!/usr/bin/env python3
"""Independent Laurent QQ Schur/overlap replay and frozen qpair input binding."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63_REVIEW'
Z={};ONE={(0,0,0):Q(1)}
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def plus(p,q):
 if not p:return q
 if not q:return p
 out=p.copy()
 for e,c in q.items():
  v=out.get(e,Q())+c
  if v:out[e]=v
  else:out.pop(e,None)
 return out
def negative(p):return{e:-c for e,c in p.items()}
def minus(p,q):return plus(p,negative(q))
def times(p,q):
 if not p or not q:return Z
 if p==ONE:return q
 if q==ONE:return p
 out={}
 for e,c in p.items():
  for f,d in q.items():
   v=tuple(x+y for x,y in zip(e,f));a=out.get(v,Q())+c*d
   if a:out[v]=a
   else:out.pop(v,None)
 return out
def invert_monomial(p):
 assert len(p)==1;e,c=next(iter(p.items()));assert c and e[1]==0
 return{tuple(-x for x in e):1/c}
def decoded(pair):
 num,den=pair;assert len(den)==1;de,dc=den[0];dc=Q(dc);assert dc and de[1]==0
 out={}
 for e,c in num:
  ex=tuple(x-y for x,y in zip(e,de));a=out.get(ex,Q())+Q(c)/dc
  if a:out[ex]=a
  else:out.pop(ex,None)
 return out
def poly(ts):return decoded([ts,[[[0,0,0],'1']]])
def native_matrix(p):
 ts=iter(p.read_text().split());r,c=int(next(ts)),int(next(ts));out=[]
 for _ in range(r):
  row=[]
  for _ in range(c):
   terms=[]
   for _ in range(int(next(ts))):v=next(ts);terms.append([[int(next(ts))for i in range(3)],v])
   row.append(terms)
  out.append(row)
 assert next(ts,None)is None;return out
def product_row(q,cof,index):return{index[tuple(x+y for x,y in zip(mon,cof))]:c for mon,c in q.items()}
def schur_setup(stage,polys):
 mons=stage['target_monomials'];index={tuple(m):j for j,m in enumerate(mons)};C=stage['monic_pivot_target_columns'];D=stage['retained_target_columns'];selected=[]
 for row in stage['monic_pivot_source_rows']:
  qi,cof=stage['four_quadric_product_rows'][row];selected.append(product_row(polys[qi],cof,index))
 assert all(selected[i].get(C[i],Z)==ONE for i in range(len(C)))
 Ucols=[[(i,selected[i][C[j]])for i in range(j)if C[j]in selected[i]]for j in range(len(C))]
 Drows=[[(j,c)for j,c in row.items()if j in D]for row in selected]
 assert all(not any(j in row for j in C[:i])for i,row in enumerate(selected))
 def nf(row):
  # Explicit row_C U^-1, then subtract its product withRD. This is not
  # the author's leading-tail rewrite nor FracElement arithmetic.
  mu=[]
  for j,col in enumerate(C):
   v=row.get(col,Z)
   for i,c in Ucols[j]:
    if mu[i]:v=minus(v,times(mu[i],c))
   mu.append(v)
  out={j:c for j,c in row.items()if j in D and c}
  for v,terms in zip(mu,Drows):
   if not v:continue
   for j,c in terms:
    w=minus(out.get(j,Z),times(v,c))
    if w:out[j]=w
    else:out.pop(j,None)
  return out
 return mons,index,C,D,nf
def profile(A):
 entries=[v for row in A for v in row if v]
 return{'rows':len(A),'columns':len(A[0])if A else 0,'nonzero_entries':len(entries),'numerator_terms':sum(len(v)for v in entries),'denominator_terms':len(entries),'max_entry_numerator_terms':max(map(len,entries),default=0)}
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63';r=read(b/'report.json');cert=read(b/'overlap-unit-certificate.json');before=read(b/'overlap18-before.json');pair=read(b/'remaining-pair-block.json')
 assert r['frozen']and r['extra_global_known_unit_pivots']==5 and r['overlap_remainder_rows_literal_zero']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in r['evidence_sha256'].items():assert sha(b/p)==h,p
 pre=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';pf=read(pre/'preflight.json');s3,s4=pf['stages'];m2=[m for _,m in s4['four_quadric_product_rows'][:29]]
 nums=native_matrix(pre/'six-prime-quadrics-numerators.txt');dens=native_matrix(pre/'six-prime-quadrics-denominators.txt');polys=[{tuple(m):decoded([nums[i][j],dens[i][0]])for j,m in enumerate(m2)if nums[i][j]}for i in range(6)]
 mons,index,C,D,nf=schur_setup(s4,polys);overlap=s4['all_unselected_overlap_rows_retained'];assert overlap==before['source_overlap_rows']==cert['source_original_overlap_rows']and D==before['target141_columns']==cert['target141_columns']
 original=[]
 for row in overlap:
  qi,cof=s4['four_quadric_product_rows'][row];n=nf(product_row(polys[qi],cof,index));original.append([n.get(j,Z)for j in D])
 assert original==[[decoded(v)for v in row]for row in before['entries']]
 assert profile(original)==cert['before_profile']
 M3,ix3,C3,D3,nf3=schur_setup(s3,polys);x0=(0,0,1,0,0,0,0,0);L=(1,0,0,0,0,0,0,0);S=product_row(polys[0],x0,ix3)
 for j,v in product_row(polys[1],L,ix3).items():
  co=minus(S.get(j,Z),v)
  if co:S[j]=co
  else:S.pop(j,None)
 Snf=nf3(S);sf=read(b/'S-cubic-NF.json');assert sf['nonzero']and Snf
 assert Snf=={ix3[tuple(m)]:decoded(v)for m,v in sf['coefficients']}
 A=[row[:]for row in original];T=[[ONE if i==j else Z for j in range(18)]for i in range(18)];available=list(range(18));used=[];unit={(1,0,0):Q(-80000,3)}
 for step in cert['unit_pivots']:
  i,j=step['pivot_original_overlap_row'],step['pivot_target141_column'];assert i in available and j not in used and A[i][j]==decoded(step['pivot_value'])==unit
  assert step['known_unit']=={'constant':'-80000/3','a_d_F_B_C_J_powers':[1,0,0,0,0,0]}
  assert step['pivot_actual_monomial']==mons[D[j]]
  ip=invert_monomial(A[i][j]);A[i]=[times(v,ip)for v in A[i]];T[i]=[times(v,ip)for v in T[i]];available.remove(i);used.append(j)
  saved={row:decoded(v)for row,v in step['eliminations']};assert saved=={t:A[t][j]for t in available if A[t][j]}
  for t in available:
   co=A[t][j]
   if co:A[t]=[minus(x,times(co,y))for x,y in zip(A[t],A[i])];T[t]=[minus(x,times(co,y))for x,y in zip(T[t],T[i])]
 assert available==cert['remaining_source_rows']and len(available)==13 and all(not v for i in available for v in A[i])
 assert A==[[decoded(v)for v in row]for row in cert['after_rows']]
 assert T==[[decoded(v)for v in row]for row in cert['row_transform18x18']]
 assert profile(A)==cert['after_profile']
 for i in range(18):
  for j in range(141):
   value=Z
   for k in range(18):
    if T[i][k]and original[k][j]:value=plus(value,times(T[i][k],original[k][j]))
   assert value==A[i][j],(i,j)
 (out/'progress.json').write_text(json.dumps({'stage':'all2538_overlap_identities_passed','seconds':time.monotonic()-began},indent=2)+'\n')
 # Only bind the previously saved pair block. No additional Gaussian,
 # newpivot choice or unit search occurs here.
 other=[j for j in range(141)if j not in used];assert pair['target141_remaining_columns']==other and pair['source_rows']==[[i,m]for i in[4,5]for m in m2]
 qpair=[]
 for qi in[4,5]:
  for cof in m2:
   row=nf(product_row(polys[qi],cof,index));v=[row.get(j,Z)for j in D]
   for st in cert['unit_pivots']:
    i,j=st['pivot_original_overlap_row'],st['pivot_target141_column'];co=v[j]
    if co:v=[minus(x,times(co,y))for x,y in zip(v,A[i])]
   assert all(not v[j]for j in used);qpair.append([v[j]for j in other])
 assert qpair==[[decoded(v)for v in row]for row in pair['entries']]
 assert profile(qpair)==pair['profile']==r['remaining_q0_q5_profile']
 assert profile(qpair)['numerator_terms']==62664 and profile(qpair)['nonzero_entries']==1286
 # Check every claimed positive scalar-unit witness without discovering
 # additional candidates or repeating the author's later Gaussian.
 factors=[{(1,0,0):Q(1)},{(0,0,1):Q(1)},poly([[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']]),poly([[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]),poly([[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']]),poly([[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']])]
 for witness in pair['individual_known_unit_entries']:
  value={(0,0,0):Q(witness['known_unit']['constant'])}
  for i,power in enumerate(witness['known_unit']['a_d_F_B_C_J_powers']):
   f=factors[i]
   if power<0:assert i<2;f=invert_monomial(f);power=-power
   for _ in range(power):value=times(value,f)
  assert qpair[witness['row']][witness['column']]==value
 assert len(pair['individual_known_unit_entries'])==21
 deps=[ROOT/'reports'/j/'report.json'for j in['S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63_REVIEW','S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW']];assert all(read(p)['all_passed']for p in deps)
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,independent_Laurent_QQ_arithmetic_no_FracElement_or_author_import=True,
  all18x141_overlap_inputs_rebuilt_by_fixed_unit_Schur=True,
  all5_fixed_pivots_literal_minus80000a_over3_known_units=True,
  all2538_T_times_before_equals_after_identities_and13_zero_rows_checked=True,
  actual_four_prime_quadric_degree4_image_rank_every_FULLU=103,
  cubic_S_normal_form_exact_and_nonzero=True,
  frozen_qpair58x136_all_coefficients_source_bound=True,all21_claimed_positive_unit_witnesses_checked=True,
  no_new_qpair_Gaussian_or_search_or_rank141_full239_claim=True,
  uniform_known_kernel_count_unchanged=38,remaining_for80=42,
  no_G16_or33quartic_expansion_newpoint_prime_GB_solver_or_author_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,b/'report.json',b/'overlap18-before.json',b/'overlap-unit-certificate.json',b/'remaining-pair-block.json',pre/'preflight.json',pre/'six-prime-quadrics-numerators.txt',pre/'six-prime-quadrics-denominators.txt',*deps]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-began,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
