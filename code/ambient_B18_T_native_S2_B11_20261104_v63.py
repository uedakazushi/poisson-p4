#!/usr/bin/env python3
"""ONE B11 scalar-gauge correction and native-S2 representative.

No N2/R12/T2/source105 solve.  E=(L0+L1)S_actual is valid in the
zero-origin source combinations, not an assertion about individual actual X2.
"""
from fractions import Fraction as F
from pathlib import Path
import collections,hashlib,json,resource,signal,time
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS, Algebra as OldAlgebra
from ambient_B18_T_native_B11_source_preflight_20261104_v63 import BASIS,bm,rec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63'
BPATH='reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json'
EXTRA={
 'B11':BPATH,'B11_report':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/report.json',
 'gauge1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase3-gauge30-firsts.json',
 'J30':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
 'native_source1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
 'native_adapter1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'native_array1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'L2_cancel':'reports/AMBIENT_B18_T_NATIVE_ADAPTER2_CANCELLATION_PREFLIGHT_20261104_v63/report.json'}
LINEAR=[(0,0),(1,8),(3,1),(4,2),(5,5),(6,6)]
QUAD=[2,7,8,9,10]
def dense(a,n=11):
 z=[F(0)]*n
 for i,c in a:z[int(i)]=F(c)
 return tuple(z)
def first_to_B(p):
 v=[F(0)]*11
 for bi,di in LINEAR:v[bi]=F(p[di])
 return tuple(v)
def encode(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [encode(v) for v in x]
 return x
def counts(x):
 if isinstance(x,F):return 1
 if isinstance(x,dict):return sum(counts(v) for v in x.values())
 if isinstance(x,(list,tuple)):return sum(counts(v) for v in x)
 return 0

class Finite:
 def __init__(self,m=0,nus=()):
  self.m=m;self.b=1<<m;self.n=11*self.b;self.Z=(F(0),)*self.n;self.ONE=self.c(1);self.nus=nus;self.nilcap=9;self.table={};self.inverse_checks=0
  for i in range(self.n):
   ai,ti=divmod(i,self.b)
   for j in range(self.n):
    aj,tj=divmod(j,self.b);v=bm({ai:F(1)},{aj:F(1)});carry=ti&tj;mask=ti^tj
    for k in range(m):
     if carry&(1<<k):v=bm(v,{a:c for a,c in enumerate(nus[k]) if c})
    if v:self.table[i,j]=[(self.b*a+mask,c) for a,c in v.items()]
 def c(self,x):return (F(x),)+(F(0),)*(self.n-1)
 def add(self,*ps):
  v=[F(0)]*self.n
  for p in ps:
   for i,c in enumerate(p):
    if c:v[i]+=c
  return tuple(v)
 def neg(self,p):return tuple(-c for c in p)
 def scale(self,p,c):return tuple(c*a for a in p)
 def mul2(self,p,q):
  if p==self.Z or q==self.Z:return self.Z
  if p==self.ONE:return q
  if q==self.ONE:return p
  v=[F(0)]*self.n;qq=[(j,c) for j,c in enumerate(q) if c]
  for i,c in enumerate(p):
   if c:
    for j,d in qq:
     for k,a in self.table.get((i,j),()):v[k]+=c*d*a
  return tuple(v)
 def mul(self,*ps):
  v=self.ONE
  for p in ps:v=self.mul2(v,p)
  return v
 def power(self,p,n):
  if n<0:return self.power(self.inv(p),-n)
  if n>=self.nilcap and not p[0]:return self.Z
  v=self.ONE
  while n:
   if n&1:v=self.mul2(v,p)
   n//=2
   if n:p=self.mul2(p,p)
  return v
 def inv(self,p):
  assert p[0],('nonunit',p);q=list(p);q[0]=F(0);q=self.scale(tuple(q),-1/p[0]);v=self.ONE;t=self.ONE
  for _ in range(8):t=self.mul(t,q);v=self.add(v,t)
  v=self.scale(v,1/p[0]);assert self.mul(p,v)==self.ONE;self.inverse_checks+=1;return v
 def dual(self,p):
  v=[F(0)]*self.n
  for i,c in enumerate(p):v[self.b*i]=F(c)
  return tuple(v)
 def tau(self,k):
  v=[F(0)]*self.n;v[1<<k]=F(1);return tuple(v)
 def install(self):
  for k in ['c','add','neg','scale','mul','power','inv']:setattr(H,k,getattr(self,k))
  H.Z=self.Z;H.ONE=self.ONE;H.TA=self.tau(0) if self.m else self.Z;H.TB=self.tau(1) if self.m>1 else self.Z

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0;bytes_saved=0
 paths={**OLD_PATHS,**EXTRA}
 def save(name,d,live=()):
  nonlocal saved,bytes_saved
  n=counts(d);assert saved+n+counts(live)<150000,('QQguard',name,saved,n,counts(live))
  text=json.dumps(encode(d),separators=(',',':'))+'\n';assert bytes_saved+len(text.encode())<32<<20
  p=out/name;p.write_text(text);saved+=n;bytes_saved+=len(text.encode())
  (out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ_slots':saved,'bytes':bytes_saved,'seconds':time.monotonic()-start})+'\n')
  print(json.dumps({'phase':name,'seconds':time.monotonic()-start}),flush=True)
  return {'path':str(p.relative_to(ROOT)),'sha256':H.sha(p)}
 (out/'phase.json').write_text(json.dumps({'phase':'authorized_worker_persisted_before_math','coefficient_cap':150000,'no_N2_or_new105_solve':True})+'\n')
 D={k:H.read(ROOT/p) for k,p in paths.items() if not p.endswith('.gz')}
 assert D['B11_report']['all_passed'] and D['L2_cancel']['all_passed']
 S,I,P,K,prep,Op,E=(D[k] for k in ['source','interface','proper','second5','preparation','Q_operator','normalizer'])
 B=D['B11'];first={int(i):dense(v) for i,v in B['all105_actual_B11_root_values']};nu=[dense(v) for v in B['nu_A_B_actual_B11']]
 pf=D['preflight'];nodes=S['DAG_scalar_nodes'];order=pf['dual9_source_evaluation']['closure_order']
 roots=[i for f in S['families'].values() for pair in f['complete_normalized_physical_w_coefficients'][0] for i in pair]
 uses=collections.Counter()
 for i in order:
  if nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
 uses.update(roots);base=Finite();vals={}
 for i in order:
  n=nodes[i];o=n[0]
  if o=='v':v=first[i]
  elif o=='c':v=base.c(F(n[1],n[2]))
  elif o=='a':v=base.add(vals[n[1]],vals[n[2]])
  elif o=='m':v=base.mul(vals[n[1]],vals[n[2]])
  elif o=='n':v=base.neg(vals[n[1]])
  elif o=='i':v=base.inv(vals[n[1]])
  else:raise ValueError(n)
  vals[i]=v
  if o in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j]
 pairs={f:[[vals[a],vals[b]] for a,b in q['complete_normalized_physical_w_coefficients'][0]] for f,q in S['families'].items()}
 oldpairs=H.read(ROOT/'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json')['physical150_tracepairs_first9']
 for family,ps in pairs.items():
  for newpair,oldpair in zip(ps,oldpairs[family]):
   for v,old in zip(newpair,oldpair):assert all(v[bi]==F(old[di]) for bi,di in LINEAR)
 p1=save('phase1-source-B11-values.json',{'physical150_tracepairs_B11':pairs,'nu_A_B_B11':nu,'all_source_firsts_match_native128':True},live=first)
 del vals,oldpairs
 positions=I['actual_group15_positions'];GL=pf['GL18_source_parameter_positions'];rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices'];centers=list(map(F,P['actual_CA_source19_center']))
 def frame_family(a,family,tauindex):
  a.install();f=S['families'][family];tau=a.tau(tauindex);ss=a.dual(first[f['p_quadratic_relation']['s']]);p=H.add(tau,H.scale(ss,F(1,2)));G=H.ident(5)
  for name,(aa,bb) in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    index=int(name[5:]);i,j=divmod(GL[index-1],5);G[i][j]=H.add(G[i][j],a.dual(first[aa]),H.mul(a.dual(first[bb]),p))
  core=H.ident(5)
  if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
  else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
  return H.mm(G,core)
 def physical(a,family,tauindex):
  a.install();return [H.add(a.dual(x),H.mul(a.dual(y),a.tau(tauindex))) for x,y in pairs[family]]
 def gauge(a,ev,od):
  a.install();G=H.ident(5)
  for (i,j),x,y in zip(positions,ev,od):G[i][j]=H.add(G[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
  return G
 def constraints(a,G,HA,PA):
  a.install();F0=H.frame(H.mm(G,HA));tr=H.mm(F0,G);moved=H.physical_push(PA,tr,H.mi(tr));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in ar]);own=[v[0] for v in own]
  return [F0[1][0],F0[1][3],F0[1][4],F0[0][3],F0[0][4]]+[H.add(own[j],H.c(-centers[j])) for j in fixed]
 ev=[first_to_B(v) for v in D['gauge1']['even15_first9']];od=[first_to_B(v) for v in D['gauge1']['odd15_first9']]
 a22=Finite(1,[nu[0]]);HA=frame_family(a22,'CA',0);PA=physical(a22,'CA',0);res=constraints(a22,gauge(a22,ev,od),HA,PA)
 assert all(not v[2*b+t] for v in res for b in range(11) if b not in QUAD for t in range(2))
 rhs=[[res[i][2*b+t] for b in QUAD] for t in range(2) for i in range(15)];Ji=[[F(c) for c in row] for row in D['J30']['J30_inverse']]
 corr=[[-sum(Ji[i][j]*rhs[j][k] for j in range(30)) for k in range(5)] for i in range(30)]
 for i in range(15):
  x=list(ev[i]);y=list(od[i])
  for k,b in enumerate(QUAD):x[b]+=corr[i][k];y[b]+=corr[15+i][k]
  ev[i]=tuple(x);od[i]=tuple(y)
 assert all(v==a22.Z for v in constraints(a22,gauge(a22,ev,od),HA,PA))
 p2=save('phase2-gauge30-B11-correction.json',{'quadratic_base_indices5':QUAD,'RHS30_by5':rhs,'correction30_by5':corr,'even15_B11':ev,'odd15_B11':od,'all30_scalar_B_constraints_exact':True,'fixed_origin_J30_only':True},live=pairs)
 del HA,PA,res,rhs,Ji,corr
 a=Finite(2,nu);G=gauge(a,ev,od);HA=frame_family(a,'CA',0);HB=frame_family(a,'CB',1);FA=H.frame(H.mm(G,HA));assert FA==H.ident(5)
 FB=H.frame(H.mm(G,HB));PA=physical(a,'CA',0);PB=physical(a,'CB',1);moved=H.physical_push(PA,G,H.mi(G));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in ar]);own=[v[0] for v in own]
 assert all(own[j]==H.c(centers[j]) for j in fixed)
 actual=[H.add(own[j],H.c(-centers[j])) for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
 kappa=[actual[j] for j in I['retained_source9_indices']]
 def match1(new,old):return all(new[4*bi+t]==F(old[4*di+t]) for bi,di in LINEAR for t in range(4))
 assert all(match1(v,q) for v,q in zip(actual,D['native_source1']['actual14']))
 assert all(match1(v,q) for v,q in zip(kappa,D['native_source1']['kappa9']))
 F12=[H.add(PA[j],H.neg(PB[j])) for j in [18,19,28,29,60,67,72,63,64,68,73,74]]
 assert all(match1(v,q) for v,q in zip(F12,D['native_source1']['physicalF12']))
 core=[kappa[j] for j in [0,1,2,6]];du=H.add(a.ONE,H.neg(H.mul(kappa[1],kappa[1])));pu=H.add(H.c(5),H.mul(kappa[1],kappa[1]));inds=[K['source9_variable_order'].index(name) for name in K['selected5_pivot_order']];rr=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
 gs=[H.add(kappa[j],H.neg(v)) for j,v in zip(inds,rr)]+[H.evalpoly(p,core) for p in prep['generators3']]
 assert all(match1(v,q) for v,q in zip(gs,D['native_adapter1']['G8']))
 p3=save('phase3-actual-proper-G8-B11.json',{'actual14_C_B44':actual,'kappa9_C_B44':kappa,'physicalF12_C_B44':F12,'actualG8_C_B44':gs,'all_constants_and_firsts_match_frozen_native128':True,'actual_CA_gauge_frame_identity':True,'no_N2_R12_T_second_coefficients_evaluated':True},live=pairs)
 del PA,PB,HA,HB,FA,FB,G,moved,own,actual,kappa,core,du,pu,rr,F12,ev,od
 Sn=[[a.Z]*32 for _ in range(12)]
 for i,j,g,s in Op['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=H.scale(gs[g],F(s))
 for i,j,c in Op['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=H.c(c)
 # The frozen L0+L1 is a SCALAR48 matrix.  Truncating nativeN instead would
 # accidentally retain some monic-induced L2 terms and is not this gauge.
 oldnu=[[F(c) for c in v] for v in pf['nu_A_B_first9']];ad=OldAlgebra(8,2,oldnu);oldN=D['native_adapter1']['N'];L={}
 for j in range(48):
  jn,b=divmod(j,4);unit=[F(0)]*36;unit[b]=F(1)
  for i in range(12):
   v=ad.mul(tuple(map(F,oldN[jn][i])),tuple(unit))
   for t in range(4):
    p=first_to_B([v[4*k+t] for k in range(9)])
    if any(p):L[4*i+t,j]=p
 a.install();scols=[];ecols=[];oldcols=D['native_array1']['source_columns128_rows48_sparse_dual9'];nnzS=nnzE=0
 for j in range(128):
  jn,b=divmod(j,4);unit=[F(0)]*44;unit[b]=F(1);sv={}
  for i in range(12):
   v=a.mul(Sn[i][jn],tuple(unit))
   for t in range(4):
    p=tuple(v[4*k+t] for k in range(11))
    if any(p):sv[4*i+t]=p
  evv={}
  for (i,k),p in L.items():
   if k in sv:
    value=base.mul(p,sv[k]);evv[i]=base.add(evv.get(i,base.Z),value)
  evv={i:p for i,p in evv.items() if any(p)}
  old={i:{int(k):F(c) for k,c in p} for i,p in oldcols[j]}
  for i in range(48):
   p=evv.get(i,base.Z);q=old.get(i,{})
   assert all(p[bi]==q.get(di,F(0)) for bi,di in LINEAR)
  scols.append([[i,[[k,c] for k,c in enumerate(p) if c]] for i,p in sorted(sv.items())]);ecols.append([[i,[[k,c] for k,c in enumerate(p) if c]] for i,p in sorted(evv.items())]);nnzS+=sum(sum(bool(c) for c in p) for p in sv.values());nnzE+=sum(sum(bool(c) for c in p) for p in evv.values())
 p4=save('actual-native128-effective-S2-B11.json',{'source_columns128_rows48_sparse_B11_effective_second':ecols,'actual_S_source_columns128_rows48_sparse_B11':scols,
  'B11_basis_order':B['B11_basis_names'],'quadratic_B_indices':QUAD,'quadratic_pairs_in_original8':[[7,7],[0,7],[1,7],[4,7],[5,7]],
  'effective_factor':'(L0+L1)*S_actual_B11, scalar L firstjet frozen','individual_actual_native_X2_claimed':False,
  'valid_for_actual_zero_origin_source_combinations_only_by_L2_cancellation':True,'all_constants_and_firsts_match_actual_native128':True,
  'actual_G8_and_monic_values_pinned':p3,'no_source105_or_physical_jet_solve':True},live=[L,Sn,pairs])
 report={'id':STEM,'all_passed':True,'actual_B11_gauge_quadratic_unknowns':150,'all30_gauge_constraints_exact':True,
  'actual_proper_G8_and_native_S2_evaluated':True,'effective_native128_second_representative_saved':True,
  'individual_actual_native_X2_computed':False,'entire_scalar_L2_omitted_by_proved_zero_origin_combination_identity':True,
  'all_available_constant_first_source_F12_G8_and_native_output_checks':True,'N2_R12_T_second_arrays_not_computed':True,
  'output_S_nonzero_QQ':nnzS,'output_effective_X_nonzero_QQ':nnzE,'data_evidence':p4,'phase_evidence':[p1,p2,p3],
  'input_hashes':{p:H.sha(ROOT/p) for p in paths.values()},'code_sha256':H.sha(Path(__file__)),
  'new_physical_source_jet_105solve_or_rank_search':False,'T_length_or_depth_decided':False,
  'total_saved_QQ_slots':saved,'total_bytes':bytes_saved,'seconds':time.monotonic()-start,'rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
