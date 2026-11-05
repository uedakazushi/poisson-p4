#!/usr/bin/env python3
"""Independent source/fixed-receipt replay of raw33 and raw45 actions."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','value':f'reports/{A}/actual-ell2-product33-and-raw45-actions.json',
 'input':f'reports/{A}/source-and-size-preflight.json','flatten':f'reports/{A}/actual-unit-flattening-through3.json',
 'source':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'lemma_peer':'reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def sp(p):return {int(i):Q(c) for i,c in p if Q(c)}
def enc(p):return [[i,str(c)] for i,c in sorted(p.items())]
def add(*ps):
 z={}
 for p in ps:
  for k,c in p.items():
   a=z.get(k,Q(0))+c
   if a:z[k]=a
   else:z.pop(k,None)
 return z
def scale(p,c):return {i:c*a for i,a in p.items() if c*a}
def mul(a,b):
 z={}
 for (r,e),c in a.items():
  for (s,f),d in b.items():
   if e+f<2 and r+s+e+f<=5:z=add(z,{(r+s,e+f):c*d})
 return z
def sub(p,x,y):
 out={}
 for mon,c in p:
  v={(0,0):Q(c)}
  for _ in range(mon[0]):v=mul(v,x)
  for _ in range(mon[1]):v=mul(v,y)
  out=add(out,v)
 return out
def mv(A,v):return add(*(scale(A[i],c) for i,c in v.items()))
def mp(A,B):return [mv(A,v) for v in B]
def eye(n):return [{i:Q(1)} for i in range(n)]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'SOURCE_AND_FIXED_RANK_RECEIPTS_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['lemma_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 S=D['source'];E=D['value'];ambient=[tuple(m) for m in E['basis74']];assert len(ambient)==74 and len(set(ambient))==74
 expected=sorted([(sector,a,e) for sector in range(8) for e in range(2) for a in range(6-e-(sector!=0))],key=lambda m:(sum(m[1:])+(m[0]!=0),-m[0],m[2],m[1]))
 assert ambient==expected;ix={m:i for i,m in enumerate(ambient)}
 plane2=[[[m[7:],c] for m,c in p if not any(m[:7])] for p in S['beta2_before_ell']]
 plane3=[[[m[7:],c] for m,c in p if not any(m[:7])] for p in S['beta3_before_ell']]
 r={(1,0):Q(1)};ell={(0,1):Q(1)};s=add(r,scale(ell,Q(1,5)))
 beta=[add(sub(a,r,s),sub(b,r,s)) for a,b in zip(plane2,plane3)]
 gA=sub(S['gA_plane3'],r,s);bot=[[Q(c) for c in row[7:]] for row in S['YB_in_YA_linear'][7:]]
 rb=add(scale(r,bot[0][0]),scale(s,bot[0][1]));sb=add(scale(r,bot[1][0]),scale(s,bot[1][1]));gB=sub(S['gB_own_plane3'],rb,sb)
 assert [[list(m),str(c)] for m,c in sorted(gA.items())]==D['input']['gA'] and [[list(m),str(c)] for m,c in sorted(gB.items())]==D['input']['gB']
 assert [[[list(m),str(c)] for m,c in sorted(b.items())] for b in beta]==D['input']['beta7_after_s_equals_r_plus_ell_over5']
 def sector(i,p):return {ix[i,a,e]:c for (a,e),c in p.items() if (i,a,e) in ix}
 base=[]
 for i in range(7):
  for j in range(i+1,7):base.append((['pair',i,j],add(sector(i+1,beta[j]),scale(sector(j+1,beta[i]),-1))))
 base.extend((['gA',i],add(sector(i+1,gA),sector(0,mul(gA,beta[i])))) for i in range(7))
 base.extend((['gB',i],sector(i+1,gB)) for i in range(7))
 def plane_times(v,p):
  out={}
  for k,c in v.items():
   sec,a,e=ambient[k]
   for (b,f),d in p.items():
    if (sec,a+b,e+f) in ix:out=add(out,{ix[sec,a+b,e+f]:c*d})
  return out
 def ztimes(v,j):
  out={}
  for k,c in v.items():
   sec,a,e=ambient[k]
   out=add(out,sector(j+1,{(a,e):c}) if sec==0 else sector(j+1,scale(mul(beta[sec-1],{(a,e):Q(1)}),-c)))
  return out
 rows=[];labels=[]
 for label,p in base:
  if not p:continue
  deg=min(ambient[k][1]+ambient[k][2]+(ambient[k][0]!=0) for k in p)
  for e in range(2):
   for a in range(6-deg-e):
    q=plane_times(p,{(a,e):Q(1)})
    if q:rows.append(q);labels.append({'relation':label,'plane_multiplier':[a,e]})
 assert len(rows)==126 and list(map(enc,rows))==E['relation_rows_sparse'] and labels==E['relation_row_sources']
 # Fixed receipts reconstruct41 independent normalized rows, without any
 # new pivot selection or Gaussian quotient computation.
 basis={}
 for receipt in E['rank41_receipts']:
  q=rows[receipt['source_row']]
  for p,c in receipt['subtract']:q=add(q,scale(basis[p],-Q(c)))
  p=receipt['pivot'];assert p not in basis and q.get(p)==Q(receipt['divide_by'])
  assert all(k>=p for k in q)
  q=scale(q,1/Q(receipt['divide_by']));assert enc(q)==receipt['normalized'];basis[p]=q
 assert len(basis)==41
 free=E['free_column_indices'];assert free==[i for i in range(74) if i not in basis] and len(free)==33
 normal=list(map(sp,E['all74_normal_forms']));fidx={k:i for i,k in enumerate(free)}
 def nf(v):return add(*(scale(normal[k],c) for k,c in v.items()))
 for k in range(74):
  if k in fidx:assert normal[k]=={fidx[k]:Q(1)}
  else:assert nf(basis[k])=={}
 assert all(not nf(p) for p in rows)
 # This supplies both rank bounds:41 independent source rows plus33
 # well-defined independent quotient images of free monomials.
 actions=[[nf(ztimes({k:Q(1)},i)) for k in free] for i in range(7)]
 actions += [[nf(plane_times({k:Q(1)},p)) for k in free] for p in [r,ell]]
 assert [[enc(p) for p in M] for M in actions]==E['coordinate_actions9_sparse_columns']
 assert all(mp(a,b)==mp(b,a) for a in actions for b in actions) and all(not p for p in mp(actions[8],actions[8]))
 assert all(not nf(ztimes(p,i)) for label,p in base for i in range(7))
 def polynomial_action(p):
  out=[{} for _ in range(33)]
  for (a,e),c in p.items():
   T=eye(33)
   for _ in range(a):T=mp(actions[7],T)
   for _ in range(e):T=mp(actions[8],T)
   out=[add(x,scale(y,c)) for x,y in zip(out,T)]
  return out
 IA=actions[:7]+[polynomial_action(gA)]
 IB=[[add(a,b) for a,b in zip(actions[i],polynomial_action(beta[i]))] for i in range(7)]+[polynomial_action(gB)]
 assert all(all(not p for p in mp(a,b)) for a in IA for b in IB)
 degree=lambda k:ambient[k][1]+ambient[k][2]+(ambient[k][0]!=0)
 assert all(degree(free[j])>=degree(k) for k,p in enumerate(normal) for j in p)
 assert all(degree(free[k])>=degree(free[j])+1 for A in actions for j,p in enumerate(A) for k in p)
 assert max(degree(k) for k in free)==4
 branches=[]
 for side,g,claimed in zip(['CA','CB'],[gA,gB],E['canonical_branch_quotients6']):
  b6=[(a,e) for e in range(2) for a in range(3)];bi={m:i for i,m in enumerate(b6)};c=-g.get((2,1),Q(0))/g[3,0]
  def red(p):
   z={}
   for (a,e),v in p.items():
    if e>=2 or a>=4 or a>=3 and e:continue
    if a==3:a,e,v=2,1,v*c
    if v:z=add(z,{bi[a,e]:v})
   return z
  def bmul(p,q):return red(mul({b6[i]:v for i,v in p.items()},{b6[j]:v for j,v in q.items()}))
  zz=[{} for _ in range(7)] if side=='CA' else [red(scale(p,-1)) for p in beta]
  BA=[[bmul(p,{i:Q(1)}) for i in range(6)] for p in zz+[red(r),red(ell)]]
  Qmap=[]
  for k in free:
   sec,a,e=ambient[k];v=red({(a,e):Q(1)})
   if sec:v=bmul(zz[sec-1],v)
   Qmap.append(v)
  assert [enc(p) for p in Qmap]==claimed['quotient_map6_by33_columns']
  assert [[enc(p) for p in A] for A in BA]==claimed['coordinate_actions9_sparse_columns']
  assert all(mp(Qmap,A)==mp(B,Qmap) for A,B in zip(actions,BA))
  piv={}
  for receipt in claimed['surjectivity_receipts']:
   v=Qmap[receipt['source_column']]
   for k,a in receipt['subtract']:v=add(v,scale(piv[k],-Q(a)))
   k=receipt['pivot'];assert k not in piv and v.get(k)==Q(receipt['divide_by'])
   piv[k]=scale(v,1/Q(receipt['divide_by']));assert enc(piv[k])==receipt['normalized']
  assert len(piv)==6;branches.append((BA,Qmap))
 raw=[]
 for i in range(9):
  M=[];offset=0
  for block in [branches[0][0][i],branches[1][0][i],actions[i]]:
   M.extend([{k+offset:c for k,c in p.items()} for p in block]);offset+=len(block)
  raw.append(M)
 assert [[enc(p) for p in M] for M in raw]==E['raw45_coordinate_actions9_sparse_columns']
 assert [enc(add(a,scale(b,Q(1,5)))) for a,b in zip(actions[7],actions[8])]==E['sA_action_sparse_columns']
 hf=[sum(degree(k)==i for k in free) for i in range(5)]
 ep=o/'independent-source126-rows-fixed41-quotient33-and-raw45-actions.json'
 dump(o,ep.name,{'inputs':inputs,'all126_source_rows74_reconstructed':True,'all41_fixed_rank_receipts_checked_no_pivot_search':True,
  'all74_normal_forms_independent_and_complete':True,'actual_product_dimension':33,'all9_actions_and64products_checked':True,
  'both_branch6_maps_and_surjectivity_receipts_intertwiners_checked':True,'all9_block_raw45_actions_verified':True,
  'source_sA_rA_ell_frame_exact':True,'degree_filtration_HF':hf,'m5_product_zero_from_strict_degree_actions':True,
  'raw_cokernel_dimension':45,'physical12_actions_C2_T2_H0_not_computed':True,'new_source_or_implicit_coefficients':0})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'actual_product33_andraw45_source_actions_independently_accepted':True,'fixed_rank_receipts_only_no_new_gaussian':True,
  'actual_product_HF':hf,'actual_product_m5_zero':True,'physical12_C2_T2_H0_computed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
