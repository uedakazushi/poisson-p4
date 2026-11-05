#!/usr/bin/env python3
"""ONE authorized120 canonical gauge2 extension. ZERO new physical jets."""
import collections,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_honest_nil_T_gauge120_source_preflight_20261105_v63 import PATHS as PF_PATHS
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63'
MON=[];IDX={}
def counts(x):
 if isinstance(x,F):return 1
 if isinstance(x,dict):return sum(counts(v)for v in x.values())
 if isinstance(x,(list,tuple)):return sum(counts(v)for v in x)
 return 0
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def smul(p,q):
 z={}
 for i,c in p.items():
  for j,d in q.items():
   k=IDX.get(tuple(sorted(MON[i]+MON[j])))
   if k is not None:z[k]=z.get(k,F(0))+c*d
 return {i:c for i,c in z.items()if c}
class Finite:
 def __init__(self,m=0,nus=()):
  self.m=m;self.b=1<<m;self.n=len(MON)*self.b;self.Z=(F(0),)*self.n;self.ONE=self.c(1);self.table={}
  for i in range(self.n):
   ai,ti=divmod(i,self.b)
   for j in range(self.n):
    aj,tj=divmod(j,self.b);v=smul({ai:F(1)},{aj:F(1)});carry=ti&tj;mask=ti^tj
    for k in range(m):
     if carry&(1<<k):v=smul(v,{a:F(c)for a,c in enumerate(nus[k])if c})
    if v:self.table[i,j]=[(self.b*a+mask,c)for a,c in v.items()]
 def c(self,c):return (F(c),)+(F(0),)*(self.n-1)
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
  v=[F(0)]*self.n;qq=[(j,c)for j,c in enumerate(q)if c]
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
  if n>=10 and not p[0]:return self.Z
  v=self.ONE
  while n:
   if n&1:v=self.mul2(v,p)
   n//=2
   if n:p=self.mul2(p,p)
  return v
 def inv(self,p):
  assert p[0];h=list(p);h[0]=F(0);h=self.scale(tuple(h),-1/p[0]);v=self.ONE;t=self.ONE
  for _ in range(10):t=self.mul(t,h);v=self.add(v,t)
  v=self.scale(v,1/p[0]);assert self.mul(p,v)==self.ONE;return v
 def dual(self,p):
  v=[F(0)]*self.n
  for i,c in enumerate(p):v[self.b*i]=F(c)
  return tuple(v)
 def tau(self,k):
  v=[F(0)]*self.n;v[1<<k]=F(1);return tuple(v)
 def install(self):
  for k in ['c','add','neg','scale','mul','power','inv']:setattr(H,k,getattr(self,k))
  H.Z=self.Z;H.ONE=self.ONE;H.TA=self.tau(0)if self.m else self.Z;H.TB=self.tau(1)if self.m>1 else self.Z
def main():
 global MON,IDX
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0;peak=0
 paths=dict(PF_PATHS);paths['preflight']='reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_SOURCE_PREFLIGHT_20261105_v63/report.json';paths['preflight_input']='reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_SOURCE_PREFLIGHT_20261105_v63/four-pair-base25-and-last-m1-contraction-preflight.json';paths['firstpreflight']='reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json'
 D={k:json.loads((ROOT/p).read_text())for k,p in paths.items()};assert D['preflight']['all_passed']and D['preflight']['aggregate_QQ_upper']==86820
 for p,h in D['preflight']['input_sha256'].items():assert sha(ROOT/p)==h
 MON=[tuple(m)for m in D['preflight_input']['base25_original8_monomials']];IDX={m:i for i,m in enumerate(MON)};assert len(MON)==25
 def save(name,d,live=()):
  nonlocal saved,peak
  q=counts(d);peak=max(peak,saved+q+counts(live));assert peak<150000,('WHOLE_QQ_GUARD',name,peak)
  b=json.dumps(enc(d),separators=(',',':'))+'\n';assert len(b.encode())+sum(p.stat().st_size for p in out.iterdir())<32<<20
  (out/name).write_text(b);saved+=q;(out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ':saved,'saved_plus_live_peak':peak,'seconds':time.monotonic()-t})+'\n')
 save('source-manifest.json',{'input_sha256':{p:sha(ROOT/p)for p in paths.values()},'authorized_new_gauge_coefficients':120,'new_physical_source_coefficients':0})
 b=Finite();mons=[tuple(m)for m in D['source2']['source2_monomials']];vars8=[tuple(F(int(j==IDX[(i,)]))for j in range(25))if(i,)in IDX else b.Z for i in range(8)];YA=vars8[:7]+[vars8[7],vars8[7]]
 def compose(p,args):
  v=b.Z
  for mon,c in zip(mons,p):
   if F(c):v=b.add(v,b.scale(b.mul(*(args[i]for i in mon)),F(c)))
  return v
 YB=[]
 for lin,quad in zip(D['adapter2']['YB_in_YA_linear'],D['adapter2']['YB_in_YA_quadratic']):
  v=b.Z
  for i,c in enumerate(lin):
   if F(c):v=b.add(v,b.scale(YA[i],F(c)))
  for mon,c in quad:
   if F(c):v=b.add(v,b.scale(b.mul(*(YA[i]for i,e in enumerate(mon)for _ in range(e))),F(c)))
  YB.append(v)
 own={f:{i:compose(p,YA if f=='CA'else YB)for i,p in D['source2']['both57_source2'][f]['all57_source2jets']}for f in ['CA','CB']}
 assert len(set(own['CA'])&set(own['CB']))==9 and all(own['CA'][i]==own['CB'][i]for i in set(own['CA'])&set(own['CB']))
 first=own['CA']|own['CB'];assert len(first)==105
 old20={i:{j:F(c)for j,c in v}for i,v in D['source20']['all105_actual_source2_values']}
 assert all([v[i]for i in range(20)]==[old20[j].get(i,F(0))for i in range(20)]for j,v in first.items())
 S=D['source'];nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for p in f['complete_normalized_physical_w_coefficients'][0]for i in p]+[f['actual_tracefree_nu']for f in S['families'].values()]+sum((D['source2']['both57_source2'][f]['equations57']for f in ['CA','CB']),[])
 needed=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if i in first:continue
  n=nodes[i]
  if n[0]in ['a','m']:todo.extend(n[1:3])
  elif n[0]in ['n','i']:todo.append(n[1])
  else:assert n[0]=='c',(i,n)
 uses=collections.Counter()
 for i in needed:
  if i not in first and nodes[i][0]in ['a','m','n','i']:uses.update(nodes[i][1:])
 uses.update(roots);vals={};live=0;sourcepeak=0
 for i in sorted(needed):
  n=nodes[i];o=n[0]
  v=first[i]if i in first else b.c(F(n[1],n[2]))if o=='c'else b.add(vals[n[1]],vals[n[2]])if o=='a'else b.mul(vals[n[1]],vals[n[2]])if o=='m'else b.neg(vals[n[1]])if o=='n'else b.inv(vals[n[1]])
  vals[i]=v;live+=25;sourcepeak=max(sourcepeak,live);assert sourcepeak+20000<150000
  if i not in first and o in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j];live-=25
 for f,args in [('CA',YA),('CB',YB)]:
  eq=[vals[i]for i in D['source2']['both57_source2'][f]['equations57']];assert eq[:48]==[b.Z]*48 and eq[48:]==args
 pairs={f:[[vals[x],vals[y]]for x,y in data['complete_normalized_physical_w_coefficients'][0]]for f,data in S['families'].items()};nu=[vals[f['actual_tracefree_nu']]for f in S['families'].values()]
 save('phase1-common105-and-physical-source-base25.json',{'base25_monomials':MON,'all105_actual_source2_values':[[i,v]for i,v in sorted(first.items())],'physical150_tracepairs_base25':pairs,'nu_A_B_base25':nu,'all114_source_equations_exact':True,'all105_base20_projections_exact':True,'physical_source_stream_peak_QQ':sourcepeak},live=first);del vals
 ev=[tuple(map(F,v))+((F(0),)*5)for v in D['gauge20']['actual_even15_base20']];od=[tuple(map(F,v))+((F(0),)*5)for v in D['gauge20']['actual_odd15_base20']]
 for i in range(15):
  x=list(ev[i]);y=list(od[i]);x[20]=F(D['gauge1']['even15_first9'][i][7]);y[20]=F(D['gauge1']['odd15_first9'][i][7]);ev[i]=tuple(x);od[i]=tuple(y)
 I,P,K,prep=(D[k]for k in ['interface','proper','second5','preparation']);GL=D['firstpreflight']['GL18_source_parameter_positions'];positions=I['actual_group15_positions'];rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices'];centers=list(map(F,P['actual_CA_source19_center']))
 def frame_family(a,family,tauindex):
  a.install();f=S['families'][family];p=H.add(a.tau(tauindex),H.scale(a.dual(first[f['p_quadratic_relation']['s']]),F(1,2)));G=H.ident(5)
  for name,(aa,bb)in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);G[i][j]=H.add(G[i][j],a.dual(first[aa]),H.mul(a.dual(first[bb]),p))
  core=H.ident(5)
  if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
  else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
  return H.mm(G,core)
 def physical(a,family,tauindex):
  a.install();return [H.add(a.dual(x),H.mul(a.dual(y),a.tau(tauindex)))for x,y in pairs[family]]
 def gauge(a):
  a.install();G=H.ident(5)
  for(i,j),x,y in zip(positions,ev,od):G[i][j]=H.add(G[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
  return G
 def constraints(a,G,HA,PA):
  a.install();FA=H.frame(H.mm(G,HA));tr=H.mm(FA,G);moved=H.physical_push(PA,tr,H.mi(tr));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in ar]);own=[v[0]for v in own]
  return [FA[1][0],FA[1][3],FA[1][4],FA[0][3],FA[0][4]]+[H.add(own[j],H.c(-centers[j]))for j in fixed]
 a50=Finite(1,[nu[0]]);HA=frame_family(a50,'CA',0);PA=physical(a50,'CA',0);res=constraints(a50,gauge(a50),HA,PA);new=[21,22,23,24]
 assert all(not v[2*k+j]for v in res for k in range(25)if k not in new for j in range(2))
 rhs=[[res[i][2*k+j]for k in new]for j in range(2)for i in range(15)];Ji=[[F(c)for c in row]for row in D['J30']['J30_inverse']];corr=[[-sum((Ji[i][j]*rhs[j][k]for j in range(30)),F(0))for k in range(4)]for i in range(30)]
 for i in range(15):
  x=list(ev[i]);y=list(od[i])
  for k,q in enumerate(new):x[q]+=corr[i][k];y[q]+=corr[15+i][k]
  ev[i]=tuple(x);od[i]=tuple(y)
 assert all(v==a50.Z for v in constraints(a50,gauge(a50),HA,PA))
 save('phase2-canonical120-four-pair-gauge.json',{'new_pairs4':[MON[i]for i in new],'RHS30_by4':rhs,'correction30_by4':corr,'actual_even15_base25':ev,'actual_odd15_base25':od,'all30_canonical_constraints_exact':True,'old_base20_and_all_firsts_preserved':True,'correct_actual_nu_carry_and_auxiliary_etaA_frame':True},live=[pairs,first]);del HA,PA,res,rhs,Ji,corr
 a=Finite(2,nu);G=gauge(a);HA=frame_family(a,'CA',0);HB=frame_family(a,'CB',1);assert H.frame(H.mm(G,HA))==H.ident(5);FB=H.frame(H.mm(G,HB));PA=physical(a,'CA',0);PB=physical(a,'CB',1);moved=H.physical_push(PA,G,H.mi(G));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in ar]);own=[v[0]for v in own];assert all(own[j]==H.c(centers[j])for j in fixed)
 actual=[H.add(own[j],H.c(-centers[j]))for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in I['retained_source9_indices']];F12=[H.add(PA[j],H.neg(PB[j]))for j in [18,19,28,29,60,67,72,63,64,68,73,74]];du=H.add(a.ONE,H.neg(H.mul(kappa[1],kappa[1])));pu=H.add(H.c(5),H.mul(kappa[1],kappa[1]));core=[kappa[j]for j in [0,1,2,6]];inds=[K['source9_variable_order'].index(n)for n in K['selected5_pivot_order']];rr=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
 gs=[H.add(kappa[j],H.neg(v))for j,v in zip(inds,rr)]+[H.evalpoly(p,core)for p in prep['generators3']]
 for key,v in [('actual14_C80',actual),('kappa9_C80',kappa),('physicalF12_C80',F12),('actualG8_C80',gs)]:
  assert all([p[4*i+j]for i in range(20)for j in range(4)]==list(map(F,q))for p,q in zip(v,D['proper20'][key])),key
 save('phase3-actual-proper-G8-base25.json',{'actual14_C100':actual,'kappa9_C100':kappa,'physicalF12_C100':F12,'actualG8_C100':gs,'all_base20_B11_and_first_projections_exact':True,'all_canonical_CA_frame_and_fiber_conditions_exact':True},live=[pairs,first])
 r={'id':STEM,'all_passed':True,'status':'PROVED_actual120_canonical_gauge_second_extension','new_gauge_second_QQ_coefficients':120,'new_physical_source_jet_coefficients':0,'all114_source_equations_exact':True,'all30_canonical_constraints_exact':True,'all_base20_B11_and_first_overlaps_exact':True,'no_source3_or_higher_proper_or_unknown_quotient_value':True,'saved_QQ_slots':saved,'peak_saved_plus_live_QQ':peak,'input_sha256':{p:sha(ROOT/p)for p in paths.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[123]-*.json')},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
