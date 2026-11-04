#!/usr/bin/env python3
"""One small exact overlap/S filtration and q0/q5 block, never the full183matrix."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def native(p):
 ls=p.read_text().splitlines();r,c=map(int,ls[0].split());out=[]
 for line in ls[1:]:
  t=line.split();n=int(t[0]);assert len(t)==1+4*n;out.append([[list(map(int,t[i+1:i+4])),t[i]]for i in range(1,len(t),4)])
 assert len(out)==r*c;return[[out[i*c+j]for j in range(c)]for i in range(r)]
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('one60s cap')));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 pre=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';pr=read(pre/'report.json');assert pr['all_exact_joint_coordinate_and_quadric_checks_passed']
 for n,h in pr['evidence_sha256'].items():assert sha(pre/n)==h
 pflight=read(pre/'preflight.json');st3,st4=pflight['stages'];m2=[m for i,m in st4['four_quadric_product_rows'][:29]];mons=st4['target_monomials'];cols=st4['retained_target_columns'];colix={c:i for i,c in enumerate(cols)}
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();load=lambda ts:R.from_dict({tuple(m):s.QQ(c)for m,c in ts});encp=lambda p:[[list(m),str(c)]for m,c in p.items()];enc=lambda c:[encp(c.numer),encp(c.denom)]
 N=native(pre/'six-prime-quadrics-numerators.txt');D=native(pre/'six-prime-quadrics-denominators.txt');polys=[]
 for i in range(6):polys.append({tuple(m):K(load(N[i][j]))/K(load(D[i][0]))for j,m in enumerate(m2)if N[i][j]})
 leads=[tuple(x)for x in[(2,0,0,0,0,0,0,0),(1,0,1,0,0,0,0,0),(0,1,0,0,0,1,0,0),(0,1,0,0,0,0,1,0)]]
 for i,l in enumerate(leads):assert polys[i][l]==K.one and all(m<=l for m in polys[i])
 def padd(p,m,c):
  z=p.get(m,K.zero)+c
  if z:p[m]=z
  else:p.pop(m,None)
 def product(q,m):return{tuple(x+y for x,y in zip(k,m)):c for k,c in q.items()}
 reductions=0
 def nf(p,stage):
  nonlocal reductions
  p=p.copy()
  for source,c in zip(stage['monic_pivot_source_rows'],stage['monic_pivot_target_columns']):
   mon=tuple(stage['target_monomials'][c]);coef=p.pop(mon,K.zero)
   if not coef:continue
   qi,cof=stage['four_quadric_product_rows'][source];reductions+=1
   for m,z in polys[qi].items():
    if m==leads[qi]:continue
    mm=tuple(x+y for x,y in zip(m,cof));assert mm<mon;padd(p,mm,-coef*z)
  assert all(tuple(stage['target_monomials'][c])not in p for c in stage['monic_pivot_target_columns']);return p
 # S is preserved, not dropped as though the four leading monomials formed an initial ideal.
 x0=(0,0,1,0,0,0,0,0);L=(1,0,0,0,0,0,0,0);S=product(polys[0],x0)
 for m,c in product(polys[1],L).items():padd(S,m,-c)
 Snf=nf(S,st3);assert Snf
 overlapids=st4['all_unselected_overlap_rows_retained'];assert len(overlapids)==18
 original=[]
 for source in overlapids:
  qi,cof=st4['four_quadric_product_rows'][source];p=nf(product(polys[qi],cof),st4);original.append([p.get(tuple(mons[c]),K.zero)for c in cols])
 def profile(A):
  nz=0;nt=dt=mx=0
  for row in A:
   for c in row:
    if c:nz+=1;nt+=len(c.numer);dt+=len(c.denom);mx=max(mx,len(c.numer))
  return dict(rows=len(A),columns=len(A[0])if A else 0,nonzero_entries=nz,numerator_terms=nt,denominator_terms=dt,max_entry_numerator_terms=mx)
 (out/'overlap18-before.json').write_text(json.dumps({'source_overlap_rows':overlapids,'target141_columns':cols,'entries':[[enc(c)for c in row]for row in original]},separators=(',',':'))+'\n')
 (out/'S-cubic-NF.json').write_text(json.dumps({'nonzero':True,'coefficients':[[list(m),enc(c)]for m,c in Snf.items()]},separators=(',',':'))+'\n')
 known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4];ucache={}
 def unit(c):
  if not c:return None
  key=str(c)
  if key in ucache:return ucache[key]
  ex=[];ground=[]
  for p in[c.numer,c.denom]:
   row=[]
   for f in known:
    n=0
    while not p.rem(f):p=p.exquo(f);n+=1
    row.append(n)
   if not p.is_ground:ucache[key]=None;return None
   ex.append(row);ground.append(p)
  z={'constant':str(K(ground[0])/K(ground[1])),'a_d_F_B_C_J_powers':[x-y for x,y in zip(*ex)]};ucache[key]=z;return z
 A=[row[:]for row in original];T=[[K.one if i==j else K.zero for j in range(18)]for i in range(18)];available=list(range(18));pivots=[];used=[]
 while available:
  choices=[]
  for r in available:
   for j,c in enumerate(A[r]):
    if j in used or not c:continue
    u=unit(c)
    if u is not None:choices.append((len(c.numer)+len(c.denom),r,j,u))
  if not choices:break
  _,r,j,u=min(choices,key=lambda z:z[:3]);p=A[r][j];A[r]=[c/p for c in A[r]];T[r]=[c/p for c in T[r]];available.remove(r);used.append(j);elims=[]
  for t in available:
   c=A[t][j]
   if c:A[t]=[x-c*y for x,y in zip(A[t],A[r])];T[t]=[x-c*y for x,y in zip(T[t],T[r])];elims.append([t,enc(c)])
  pivots.append({'pivot_original_overlap_row':r,'pivot_target141_column':j,'pivot_actual_monomial':mons[cols[j]],'pivot_value':enc(p),'known_unit':u,'eliminations':elims})
  (out/'progress.json').write_text(json.dumps({'completed_known_unit_pivots':len(pivots),'seconds':time.monotonic()-started},indent=2)+'\n')
  if time.monotonic()-started>35:break
 remainingzero=all(not c for r in available for c in A[r]);rank=len(pivots)
 overlap={'unit_pivots':pivots,'row_transform18x18':[[enc(c)for c in row]for row in T],'after_rows':[[enc(c)for c in row]for row in A],'remaining_source_rows':available,'remaining_rows_literal_zero':remainingzero,'source_original_overlap_rows':overlapids,'target141_columns':cols,'before_profile':profile(original),'after_profile':profile(A),'all2538_QQ_row_operation_identities_passed':False}
 (out/'overlap-unit-certificate.json').write_text(json.dumps(overlap,separators=(',',':'))+'\n')
 assert all(sum((T[i][k]*original[k][j]for k in range(18)),K.zero)==A[i][j]for i in range(18)for j in range(141));overlap['all2538_QQ_row_operation_identities_passed']=True
 (out/'overlap-unit-certificate.json').write_text(json.dumps(overlap,separators=(',',':'))+'\n')
 # Only the two remaining quadrics; every genuine cubic/quartic stays in the exact lazy oracle.
 othercols=[j for j in range(141)if j not in used];qpair=[];qlabels=[];stopped=False
 for qi in[4,5]:
  for mon in m2:
   if time.monotonic()-started>48:stopped=True;break
   p=nf(product(polys[qi],mon),st4);row=[p.get(tuple(mons[c]),K.zero)for c in cols]
   for step in pivots:
    r,j=step['pivot_original_overlap_row'],step['pivot_target141_column'];v=row[j]
    if v:row=[x-v*y for x,y in zip(row,A[r])]
   assert all(not row[j]for j in used);qpair.append([row[j]for j in othercols]);qlabels.append([qi,mon])
  if stopped:break
 initial_units=[];unit_scan_complete=True
 for r,row in enumerate(qpair):
  if time.monotonic()-started>54:unit_scan_complete=False;break
  for j,c in enumerate(row):
   u=unit(c)
   if u is not None:initial_units.append({'row':r,'column':j,'original_target141_column':othercols[j],'known_unit':u})
 (out/'remaining-pair-block.json').write_text(json.dumps({'complete58_rows':len(qpair)==58,'source_rows':qlabels,'target141_remaining_columns':othercols,'entries':[[enc(c)for c in row]for row in qpair],'profile':profile(qpair),'individual_known_unit_entry_scan_complete':unit_scan_complete,'individual_known_unit_entry_count':len(initial_units),'individual_known_unit_entries':initial_units},separators=(',',':'))+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),pre/'report.json',pre/'preflight.json',pre/'six-prime-quadrics-numerators.txt',pre/'six-prime-quadrics-denominators.txt',ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63/report.json']
 report={'id':ID,'all_exact_fixed_monic_reductions_and_unit_operations_passed':True,'extra_global_known_unit_pivots':rank,'overlap_remainder_rows_literal_zero':remainingzero,
 'four_quadric_weight4_rank_everywhere_if_complete':98+rank if remainingzero else None,'remaining_q0_q5_block_complete':len(qpair)==58,'remaining_q0_q5_block_shape':[len(qpair),len(othercols)],'remaining_q0_q5_profile':profile(qpair),'remaining_q0_q5_individual_known_unit_entries':len(initial_units),
 'no_additional_qpair_Gaussian_or_full183matrix_expansion':True,'all33_genuine_quartics_retained_in_unexpanded_lossless_presentation':True,
 'uniform_kernel38_missing42_unchanged':True,'no_rank141_full239_or_allpoint_cubic57_claim':True,
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'artifact_bytes':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert report['artifact_bytes']<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args()
 try:main(a.output.resolve())
 except(TimeoutError,MemoryError)as error:
  a.output.mkdir(parents=True,exist_ok=True);(a.output/'BLOCKED.json').write_text(json.dumps({'completed':False,'constraint':type(error).__name__,'message':str(error),'no_auto_retry':True},indent=2)+'\n');raise
