#!/usr/bin/env python3
"""ONE fixed dual-nil derivative of the actual normalizer source adapter."""
import collections,json,resource,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_native128_first_source_preflight_20261104_v63 import PATHS as PREFPATHS
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE128_FIRST_20261104_v63'
PATHS={**PREFPATHS,
 'preflight':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
 'normalizer':'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json',
 'preparation':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',
 'X0_array':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json'}

class Algebra:
 def __init__(self,e=8,m=2,nus=None):
  self.e=e;self.m=m;self.b=1<<m;self.n=(e+1)*self.b
  self.Z=(F(0),)*self.n;self.ONE=self.c(1);self.nus=nus or [[F(0)]*(e+1) for _ in range(m)]
  self.nilcap=(m+1 if not e else [2,4,5][m]);self.table={};self.inverse_checks=0
  for i in range(self.n):
   a,u=divmod(i,self.b)
   for j in range(self.n):
    b,v=divmod(j,self.b)
    if a and b:continue
    k=a+b;carry=u&v;mask=u^v
    if not carry:self.table[i,j]=[(k*self.b+mask,F(1))]
    elif k==0 and carry&(carry-1)==0:
     q=carry.bit_length()-1;self.table[i,j]=[(t*self.b+mask,c) for t,c in enumerate(self.nus[q]) if c]
  assert all(not nu[0] for nu in self.nus)
 def c(self,x):return (F(x),)+(F(0),)*(self.n-1)
 def add(self,*ps):
  z=[F(0)]*self.n
  for p in ps:
   for i,c in enumerate(p):
    if c:z[i]+=c
  return tuple(z)
 def neg(self,p):return tuple(-x for x in p)
 def scale(self,p,c):return tuple(c*x for x in p)
 def mul2(self,p,q):
  if p==self.Z or q==self.Z:return self.Z
  if p==self.ONE:return q
  if q==self.ONE:return p
  z=[F(0)]*self.n;qs=[(j,y) for j,y in enumerate(q) if y]
  for i,x in enumerate(p):
   if x:
    for j,y in qs:
     for k,c in self.table.get((i,j),()):z[k]+=x*y*c
  return tuple(z)
 def mul(self,*ps):
  z=self.ONE
  for p in ps:z=self.mul2(z,p)
  return z
 def power(self,p,n):
  if n<0:return self.power(self.inv(p),-n)
  if n>=self.nilcap and not p[0]:return self.Z
  z=self.ONE
  while n:
   if n&1:z=self.mul(z,p)
   n//=2
   if n:p=self.mul(p,p)
  return z
 def inv(self,p):
  assert p[0],('nonunit',p)
  q=list(p);q[0]=F(0);q=self.scale(tuple(q),-1/p[0]);v=self.ONE;term=self.ONE
  for _ in range(self.nilcap-1):term=self.mul(term,q);v=self.add(v,term)
  v=self.scale(v,1/p[0]);assert self.mul(p,v)==self.ONE;self.inverse_checks+=1;return v
 def dual(self,p):
  z=[F(0)]*self.n
  for i,x in enumerate(p):z[i*self.b]=F(x)
  return tuple(z)
 def tau(self,i):
  z=[F(0)]*self.n;z[1<<i]=F(1);return tuple(z)
 def install(self):
  for k in ['c','add','neg','scale','mul','power','inv']:setattr(H,k,getattr(self,k))
  H.Z=self.Z;H.ONE=self.ONE
  H.TA=self.tau(0) if self.m else self.Z
  H.TB=self.tau(1) if self.m>1 else self.Z

def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [enc(v) for v in x]
 return x
def counts(x):
 if isinstance(x,F):return 1
 if isinstance(x,dict):return sum(counts(v) for v in x.values())
 if isinstance(x,(tuple,list)):return sum(counts(v) for v in x)
 return 0
def inverseQ(a):
 n=len(a);M=[r[:]+[F(i==j) for j in range(n)] for i,r in enumerate(a)];det=F(1)
 for j in range(n):
  k=next(k for k in range(j,n) if M[k][j]);
  if k!=j:M[k],M[j]=M[j],M[k];det=-det
  c=M[j][j];det*=c;M[j]=[v/c for v in M[j]]
  for i in range(n):
   if i!=j and M[i][j]:
    c=M[i][j];M[i]=[v-c*w for v,w in zip(M[i],M[j])]
 return [r[n:] for r in M],det
def mmQ(a,b):return [[sum((x*y for x,y in zip(r,c)),F(0)) for c in zip(*b)] for r in a]

def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True);saved=0;bytes_saved=0
 def save(name,data,live=()):
  nonlocal saved,bytes_saved
  q=counts(data);assert saved+q+counts(live)<150000,('totalQQguard',saved,q,counts(live))
  text=json.dumps(enc(data),separators=(',',':'))+'\n';assert bytes_saved+len(text.encode())<32<<20
  p=out/name;p.write_text(text);saved+=q;bytes_saved+=len(text.encode())
  (out/'progress.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-start,'saved_QQ_slots':saved,'bytes':bytes_saved})+'\n')
  print(json.dumps({'phase':name,'seconds':time.monotonic()-start}),flush=True)
  return {'path':str(p.relative_to(ROOT)),'sha256':H.sha(p)}
 D={k:H.read(ROOT/p) for k,p in PATHS.items() if not p.endswith('.gz')}
 S,I,P,Q,R,K,C,prep,Op,E,A0,origin=(D[k] for k in ['source','interface','proper','equations','first5','second5','membership','preparation','Q_operator','normalizer','X0','X0_array'])
 pf=D['preflight'];first={i:tuple(map(F,p)) for i,p in pf['all105_existing_unknown_first9']}
 nu=[list(map(F,p)) for p in pf['nu_A_B_first9']]
 nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values() for pair in f['complete_normalized_physical_w_coefficients'][0] for i in pair]
 order=pf['dual9_source_evaluation']['closure_order'];uses=collections.Counter()
 for i in order:
  if nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
 uses.update(roots);a9=Algebra(8,0);vals={}
 for i in order:
  n=nodes[i];o=n[0]
  if o=='v':v=first[i]
  elif o=='c':v=a9.c(F(n[1],n[2]))
  elif o=='a':v=a9.add(vals[n[1]],vals[n[2]])
  elif o=='m':v=a9.mul(vals[n[1]],vals[n[2]])
  elif o=='n':v=a9.neg(vals[n[1]])
  elif o=='i':v=a9.inv(vals[n[1]])
  else:raise ValueError(n)
  vals[i]=v
  if o in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j]
 pairs={f:[[vals[a],vals[b]] for a,b in d['complete_normalized_physical_w_coefficients'][0]] for f,d in S['families'].items()}
 point=list(map(F,P['actual_P0_physical150']))
 assert all(a[0]==point[i] and b[0]==0 for f in pairs.values() for i,(a,b) in enumerate(f))
 p1=save('phase1-source-firsts.json',{'dual_index_order':pf['dual_index_order'],'physical150_tracepairs_first9':pairs,'nu_first9':nu},live=first)
 del vals
 odd=list(map(F,Op['group30_to120']['group30_odd_origin']));positions=I['actual_group15_positions']
 rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices']
 centers=list(map(F,P['actual_CA_source19_center']));GL=pf['GL18_source_parameter_positions']
 def frame_family(alg,family,tauindex):
  alg.install();f=S['families'][family];tau=alg.tau(tauindex)
  ss=alg.dual(first[f['p_quadratic_relation']['s']]);p=H.add(tau,H.scale(ss,F(1,2)))
  G=H.ident(5)
  for name,(aa,bb) in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    index=int(name[5:]);i,j=divmod(GL[index-1],5)
    G[i][j]=H.add(G[i][j],alg.dual(first[aa]),H.mul(alg.dual(first[bb]),p))
  core=H.ident(5)
  if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
  else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
  return H.mm(G,core)
 def physical(alg,family,tauindex):
  alg.install();tau=alg.tau(tauindex)
  return [H.add(alg.dual(a),H.mul(alg.dual(b),tau)) for a,b in pairs[family]]
 def gauge_matrix(alg,ev,od):
  alg.install();G=H.ident(5);tau=alg.tau(0)
  for (i,j),a,b in zip(positions,ev,od):G[i][j]=H.add(G[i][j],alg.dual(a),H.mul(alg.dual(b),tau))
  return G
 def constraints(alg,G,HA,PA):
  alg.install();F=H.frame(H.mm(G,HA));transform=H.mm(F,G)
  moved=H.physical_push(PA,transform,H.mi(transform))
  own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in ar]);own=[q[0] for q in own]
  return [F[1][0],F[1][3],F[1][4],F[0][3],F[0][4]]+[H.add(own[j],H.c(-centers[j])) for j in fixed]
 # The sole new225 mixed coefficients, computed in length4 at the SAVEDnil center.
 ad=Algebra(0,2);ad.install();eta=ad.tau(0);delta=ad.tau(1)
 HA=H.ident(5);MA=next(r for r in I['families'] if r['family']=='CA')['source_frame5_by5_by9_first'][7]
 for i in range(5):
  for j in range(5):HA[i][j]=H.add(HA[i][j],H.scale(eta,F(MA[i][j])))
 base=H.ident(5)
 for (i,j),c in zip(positions,odd):base[i][j]=H.add(base[i][j],H.scale(eta,c))
 J0=[[F(0)]*15 for _ in range(15)];J1=[[F(0)]*15 for _ in range(15)]
 for k,(i,j) in enumerate(positions):
  G=[r[:] for r in base];G[i][j]=H.add(G[i][j],delta)
  v=constraints(ad,G,HA,[ad.c(x) for x in point])
  assert all(q[0]==q[1]==0 for q in v)
  for r,q in enumerate(v):J0[r][k]=q[2];J1[r][k]=q[3]
 assert J0==[[F(x) for x in row] for row in I['group15_unit_Jacobian']]
 Ainv,det0=inverseQ(J0);assert det0==F(-2,3)
 lower=mmQ(mmQ(Ainv,J1),Ainv);Ji=[[F(0)]*30 for _ in range(30)]
 for i in range(15):
  for j in range(15):Ji[i][j]=Ainv[i][j];Ji[15+i][15+j]=Ainv[i][j];Ji[15+i][j]=-lower[i][j]
 J=[[F(0)]*30 for _ in range(30)]
 for i in range(15):
  for j in range(15):J[i][j]=J0[i][j];J[15+i][j]=J1[i][j];J[15+i][15+j]=J0[i][j]
 assert mmQ(J,Ji)==[[F(i==j) for j in range(30)] for i in range(30)]
 p2=save('phase2-fixed30-Jacobian.json',{'A0':J0,'A1_eta_delta':J1,'J30_inverse':Ji,'det_J30':det0**2},live=pairs)
 a18=Algebra(8,1,[nu[0]]);a18.install();HA=frame_family(a18,'CA',0);PA=physical(a18,'CA',0)
 ev=[(F(0),)*9 for _ in range(15)];od=[(x,)+(F(0),)*8 for x in odd]
 v=constraints(a18,gauge_matrix(a18,ev,od),HA,PA);assert all(q[0]==q[1]==0 for q in v)
 rhs=[[v[i][2*(k+1)+b] for k in range(8)] for b in range(2) for i in range(15)]
 corr=[[-sum((Ji[i][j]*rhs[j][k] for j in range(30)),F(0)) for k in range(8)] for i in range(30)]
 ev=[(F(0),)+tuple(corr[i]) for i in range(15)];od=[(odd[i],)+tuple(corr[15+i]) for i in range(15)]
 assert all(q==a18.Z for q in constraints(a18,gauge_matrix(a18,ev,od),HA,PA))
 p3=save('phase3-gauge30-firsts.json',{'RHS30_by8':rhs,'correction30_by8':corr,'even15_first9':ev,'odd15_first9':od,'all30_constraints_exact':True},live=pairs)
 del HA,PA,base,J,Ji,J0,J1,Ainv,lower,v
 a=Algebra(8,2,nu);a.install();G=gauge_matrix(a,ev,od);HA=frame_family(a,'CA',0);HB=frame_family(a,'CB',1)
 assert H.frame(H.mm(G,HA))==H.ident(5)
 FA=H.frame(H.mm(G,HA));FB=H.frame(H.mm(G,HB));PA=physical(a,'CA',0)
 moved=H.physical_push(PA,G,H.mi(G));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in ar]);own=[x[0] for x in own]
 assert all(own[j]==H.c(centers[j]) for j in fixed)
 actual=[H.add(own[j],H.c(-centers[j])) for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
 kappa=[actual[j] for j in I['retained_source9_indices']]
 targets=[18,19,28,29,60,67,72,63,64,68,73,74]
 PB=physical(a,'CB',1);F12=[H.add(PA[j],H.neg(PB[j])) for j in targets]
 del PA,PB,moved,HA,HB,FA
 transform=H.mm(FB,G);itransform=H.mi(transform);ann=H.constmat(P['CB_annihilator17_by36']);R17=[[a.Z]*12 for _ in range(17)]
 for j,target in enumerate(targets):
  inp=[a.Z]*150;inp[target]=a.ONE
  for k,row in enumerate(E['eliminated_zero_rows4']):inp[row]=H.c(E['eliminated_four_coordinate_reconstruction'][k][rows.index(target)])
  q=H.physical_push(inp,transform,itransform);v=H.mm(ann,[[q[i]] for i in rows])
  for i in range(17):R17[i][j]=v[i][0]
 f17=[H.evalpoly(p,actual) for p in Q['equations14variable_sparse']]
 assert H.mm(R17,[[x] for x in F12])==[[x] for x in f17]
 du=H.add(a.ONE,H.neg(H.mul(kappa[1],kappa[1])));pu=H.add(H.c(5),H.mul(kappa[1],kappa[1]))
 names=['v0','v1','b','c','e'];indices=[0,1,10,11,13];roots=[]
 for name in names:
  p=R['five_explicit_rational_roots'][name];roots.append(H.mul(H.evalpoly(p['numerator'],kappa),H.power(du,-p['D_denominator_power'])))
 coeff=[H.telescope(p,actual,indices,roots) for p in Q['equations14variable_sparse']]
 asinv=H.mi(coeff[:5]);ba=H.mm(coeff[5:],asinv)
 R12=[[H.add(R17[i+5][j],H.neg(H.add(*(H.mul(ba[i][k],R17[k][j]) for k in range(5))))) for j in range(12)] for i in range(12)]
 for i,p in enumerate(R['all12_actual_remaining_and_divided_records']):
  unit=H.scale(H.power(du,p['D_denominator_power']),F(p['constant_denominator'],p['primitive_content_unit']));R12[i]=[H.mul(unit,q) for q in R12[i]]
 hpolys=[p['divided_primitive_polynomial'] for p in R['all12_actual_remaining_and_divided_records']];hv=[[H.evalpoly(p,kappa)] for p in hpolys]
 assert H.mm(R12,[[x] for x in F12])==[[H.mul(kappa[1],x[0])] for x in hv]
 p4=save('phase4-native-source-and-R12.json',{'actual14':actual,'kappa9':kappa,'physicalF12':F12,'R17':R17,'R12':R12,'all17_offintersection_identity':True,'full_first5_telescopes':True,'dH12_R12F12':True},live=pairs)
 del R17,coeff,ba,asinv,q,v,inp,transform,itransform
 core=[kappa[j] for j in [0,1,2,6]];inds=[K['source9_variable_order'].index(name) for name in K['selected5_pivot_order']];rr=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
 co=[H.telescope(p,kappa,inds,rr) for p in hpolys];perm=K['selected5_source_H_rows']+[j for j in range(12) if j not in K['selected5_source_H_rows']]
 A5=[co[j] for j in perm[:5]];B7=[co[j] for j in perm[5:]];Ai=H.mi(A5);BA=H.mm(B7,Ai);T=[[a.Z]*12 for _ in range(12)]
 for i in range(5):
  for j in range(5):T[i][perm[j]]=Ai[i][j]
 for i in range(7):
  T[i+5][perm[i+5]]=a.ONE
  for j in range(5):T[i+5][perm[j]]=H.neg(BA[i][j])
 for i,p in enumerate(K['remaining7_actual_core_polynomials']):
  lam=H.scale(H.mul(H.power(du,-p['D_denominator_power']),H.power(pu,-p['P_denominator_power'])),F(p['QQ_unit_content_numerator'],p['QQ_unit_content_denominator']));T[i+5]=[H.mul(H.inv(lam),x) for x in T[i+5]]
 E3=H.ident(12);E3[5][5]=H.inv(du);E3[6][6]=H.inv(du)
 r=H.add(core[2],H.neg(core[0]),H.neg(core[1]));s=H.add(core[0],core[2])
 E3[7][5]=H.neg(H.mul(r,H.inv(H.mul(pu,du))));E3[7][7]=H.scale(H.inv(H.mul(pu,du)),F(-2));T=H.mm(E3,T);E4=H.ident(12)
 for j in range(4):
  fac=H.power(du,[2,2,2,3][j])
  for i in range(3):E4[8+j][5+i]=H.neg(H.mul(fac,H.evalpoly(C['cofactor3_by4'][i][j],core)))
 T=H.mm(E4,T);N=H.mm(T,R12);gs=[H.add(kappa[j],H.neg(v)) for j,v in zip(inds,rr)]+[H.evalpoly(p,core) for p in prep['generators3']]
 assert H.mm(T,hv)==[[x] for x in gs]+[[a.Z] for _ in range(4)]
 assert H.mm(N,[[x] for x in F12])==[[H.mul(kappa[1],x)] for x in gs]+[[a.Z] for _ in range(4)]
 for name,matrix in [('R12',R12),('T',T),('N',N)]:
  assert [[[str(x) for x in q[:4]] for q in row] for row in matrix]==A0[name]
 assert [[str(x) for x in q[:4]] for q in gs]==A0['G8_C0']
 p5=save('phase5-native-adapter-firsts.json',{'T':T,'N':N,'G8':gs,'all_N_G8_origin_coefficients_reproduced':True,'TH12_G8':True,'NF12_dG8':True},live=pairs)
 del T,R12,co,A5,B7,Ai,BA,E3,E4,hv,actual,kappa,FB,G,own,core,rr,F12
 Sn=[[a.Z]*32 for _ in range(12)]
 for i,j,g,sign in Op['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=H.scale(gs[g],F(sign))
 for i,j,x in Op['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=H.c(x)
 native=H.mm([list(r) for r in zip(*N)],Sn);cols=[];old=origin['all48_by128_SyzQ_origins'];nnz=0
 for j in range(128):
  native_col,b=divmod(j,4);unit=[F(0)]*a.n;unit[b]=F(1);unit=tuple(unit);col=[]
  for i in range(12):
   p=H.mul(native[i][native_col],unit)
   for component in range(4):
    row=4*i+component;v=[p[4*k+component] for k in range(9)]
    assert str(v[0])==old[row][j]
    sp=[[k,x] for k,x in enumerate(v) if x]
    if sp:col.append([row,sp]);nnz+=len(sp)
  cols.append(col)
 p6=save('actual-native128-X-mod-m2.json',{'source_columns128_rows48_sparse_dual9':cols,'dual_index_order':pf['dual_index_order'],
  'native_bottom8_fullX_rows_zero':True,'all6144_origin_coefficients_reproduced':True,'physical_source_or_105_IFT_solve':False},live=[N,native,Sn,pairs])
 rep={'all_passed':True,'actual_native128_firsts_computed':True,'offintersection_full17_and_R12_T_N_identities':True,
  'fixed30_mixed_Jacobian_A1_not_assumed_zero':True,'one_exact_group30_first_correction':True,
  'all_N_G8_S_origins_match_saved_X0':True,'output_schema':'source_columns128_rows48_sparse_dual9','output_nonzero_QQ':nnz,
  'data_evidence':p6,'phase_evidence':[p1,p2,p3,p4,p5],'input_sha256':{p:H.sha(ROOT/p) for p in PATHS.values()},
  'no_new_physical_source_jet_or_105solve_or_rank_search':True,'muT_or_depth_decided':False,
  'total_saved_QQ_slots':saved,'total_output_bytes':bytes_saved,'seconds':time.monotonic()-start,
  'rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'code_sha256':H.sha(Path(__file__))}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep),flush=True)
if __name__=='__main__':main()
