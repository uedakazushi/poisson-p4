#!/usr/bin/env python3
"""ONE cache-only actual right106 cubic block; no new source/gauge coefficient."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_honest_nil_T_gauge120_cached_20261105_v63 as G
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,polynomial,pivots,invq,mv
from ambient_B18_honest_nil_T_operator_20261104_v63 import fixture_h
from ambient_B18_honest_nil_T_right106_source_preflight_20261105_v63 import PATHS as PF_PATHS
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63'
PATHS=dict(PF_PATHS)
PATHS.update({
 'preflight_report':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63/report.json',
 'preflight':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63/right106-full-first-span-and-cache-footprint.json',
 'first2':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'GL':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'filtration_peer':'reports/AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW/report.json'})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def plus(p,q,c=Q(1)):
 p=dict(p)
 for k,a in q.items():
  p[k]=p.get(k,Q(0))+c*a
  if not p[k]:del p[k]
 return p
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0;peak=0
 def save(name,d,live=()):
  nonlocal saved,peak
  q=nqq(d);peak=max(peak,saved+q+nqq(live));assert peak<150000,('AGGREGATE_GUARD',name,peak)
  blob=json.dumps(enc(d),separators=(',',':'))+'\n';assert len(blob.encode())+sum(p.stat().st_size for p in out.iterdir())<32<<20
  (out/name).write_text(blob);saved+=q;(out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ':saved,'peak_QQ':peak,'seconds':time.monotonic()-t})+'\n')
 D={k:read(ROOT/p)for k,p in PATHS.items()};PF=D['preflight'];assert D['preflight_report']['cache_gate_closed']and D['filtration_peer']['actual_m5_Tnil_zero_accepted']
 for p,h in D['preflight_report']['input_sha256'].items():assert sha(ROOT/p)==h
 assert all(not Q(c)for row in PF['all106_U_sensitivities3']+PF['all106_V_sensitivities3']for c in row)
 save('preflight.json',{'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'source_gate_QQ_upper':46320,'no_new_source_gauge_jet':True})
 MON=[tuple(m)for m in D['physics']['base25_monomials']];MI={m:i for i,m in enumerate(MON)};G.MON=MON;G.IDX=MI;nus=[[Q(c)for c in v]for v in D['physics']['nu_A_B_base25']];a=G.Finite(2,nus)
 P=D['proper'];core=[tuple(map(Q,P['kappa9_C100'][i]))for i in [0,1,2,6]];vs=[tuple(map(Q,v))for v in P['actualG8_C100'][:5]]+[polynomial(a,D['GL']['forward_GL4'][s],core)for s in ['h','c','r','l']]
 def unit(k):v=[Q(0)]*100;v[k]=Q(1);return tuple(v)
 native=[];R=D['radical']
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  pp={row:polynomial(a,p,vs)for row,p in col}
  for mask in range(4):
   v={}
   for row,q in pp.items():
    q=a.mul(q,unit(mask))
    for k in range(4):
     c={i:q[4*i+k]for i in range(25)if q[4*i+k]}
     if c:v[4*row+k]=c
   native.append(v)
 qv=[];dd=polynomial(a,R['d_polynomial'],vs)
 for pol in R['R9_generators_polynomials9']:
  q=a.mul(dd,polynomial(a,pol,vs))
  for mask in range(4):qv.append(a.mul(q,unit(mask)))
 for i in range(8):
  q=[Q(0)]*100
  if i<7 and(i,)in MI:q[4*MI[(i,)]+1]=Q(-1)
  qv.append(tuple(q))
 vg=[[tuple(v[4*i+1]for i in range(25)),tuple(v[4*i+3]for i in range(25)),tuple(v[4*i+2]-Q(11,12)*v[4*i+1]for i in range(25))]for v in qv]
 m0=[[Q(c)for c in row]for row in D['m0']['m0_24_by44']];mu1=[Q(c)for c in D['first2']['true_actual_m1_row7_r_coefficients44']];assert[v[0][MI[(7,7)]]for v in vg]==mu1
 row15=[v[1][MI[(7,7)]]for v in vg]
 Bcols=D['chain']['B48_by640_sparse_columns'];Pcols=R['actual_P8_by640_columns'];rn=R['new_Cscalar_and_R0component_nodes'];tv={i:rec(v)for i,v in D['roots2']['all1239_Broot_second_representatives']};cv={i:rec(v)for i,v in D['roots2']['all768_C0root_second_representatives']}
 AA=[]
 for row in D['X1']['new_A12_by9_C36']:
  z=[]
  for v in row:
   q=[Q(0)]*100
   for k in range(4):q[k]=Q(v[k])
   for var in range(8):
    if(var,)in MI:
     for k in range(4):q[4*MI[(var,)]+k]=Q(v[4*(var+1)+k])
   z.append(tuple(q))
  AA.append(z)
 raw=[]
 for col in D['X1']['all780_source_columns44_sparse_dual9']:raw.append({row:rec(v)for row,v in col})
 xp=D['X1']['X0_pivot_columns'];AF=[[Q(c)for c in row]for row in D['carry']['all48_actual_ResAT_row20_B11']];AF1=[[Q(c)for c in row]for row in D['first2']['all48_actual_adapter_row20_first9']]
 def xp1(p,row,var):return raw[p].get(row,{}).get(1+var,Q(0))
 def xp0(p,row):return raw[p].get(row,{}).get(0,Q(0))
 def xpscalar2(p,mon):
  if p>=640:return native[p-640].get(20,{}).get(MI[mon],Q(0))
  v=Q(0)
  for row,node in Bcols[p]:
   q=tv[node]
   if mon==(7,7):v+=AF[row][2]*q.get(0,Q(0))+AF1[row][8]*q.get(MI[(7,)],Q(0))
   else:
    i=mon[0];bi={0:7,4:9,5:10}[i];v+=AF[row][bi]*q.get(0,Q(0))+AF1[row][i+1]*q.get(MI[(7,)],Q(0))+AF1[row][8]*q.get(MI[(i,)],Q(0))
  return v
 m2xp=[-sum((mu1[j]*xp1(p,j,7)for j in range(44)),Q(0))-Q(5,2)*xpscalar2(p,(7,7))for p in xp]
 assert m2xp==list(map(Q,D['lower2']['all18_mu2_Xp0_replacement_values']))
 hp1=[[Q(5,2)*xpscalar2(p,(i,7))+sum((mu1[j]*xp1(p,j,i)for j in range(44)),Q(0))for p in xp]for i in [0,4,5]];assert hp1==[[Q(c)for c in v]for v in D['lower2']['all54_hXp1_r']]
 g3=sum(map(Q,D['C0']['cached_actual_CA_gA_all_degree3_coefficients']));assert g3==Q(125,24)
 selected=PF['selected20_contraction_basis_columns'];assert selected==list(range(18))+[53,54]
 mu2basis=m2xp+[-g3,Q(0)] # U retainedsymbol, zero LOWER representative only
 lastdict={row:list(map(Q,v))for row,v in D['lower2']['all54_last_m1_Xp0_contractions']};row23basis=lastdict[23]+[Q(11,12)*g3,Q(0)] # V likewise
 factor=[[[Q(c)for c in f]for f in row]for row in PF['all35_by8_X1_factorizations_into20_basis']]
 def dotfactor(b,p,var):return sum((c*v for c,v in zip(b,factor[p][var])),Q(0))
 def m1(row,j,var):
  sec,i=divmod(row,8);assert i<7
  if var!=7 and var>i:return Q(0)
  mon=tuple(sorted((i,var)));assert mon in MI,('UNCACHED_M1_PAIR',row,var)
  return vg[j][sec][MI[mon]]
 def mxp1(p,row,var):
  z=sum((m0[row][j]*xp1(p,j,var)for j in range(44)),Q(0))
  if row%8!=7:z+=sum((m1(row,j,var)*xp0(p,j)for j in range(44)),Q(0))
  elif var==7:return Q(0)
  return z
 projections=[(f,v,rec(pi))for f,v,pi in PF['all_nonzero_convolution_duals']];eff={(f,v):[(h,w,Q(c))for h,w,c in row]for f,v,row in PF['all_actual_effective_C1_H2_callbacks']};queries=sorted({(h,w)for row in eff.values()for h,w,c in row});pairs=list(itertools.combinations(range(8),2));HS={q:{}for q in queries}
 for sec in range(3):
  for row in range(8):
   for mi,mon in enumerate(MON):
    if len(mon)!=2 or row==7 and mon==(7,7):continue
    ex=[0]*8
    for k in mon:ex[k]+=1
    v=[{}for _ in range(8)];v[row]={tuple(ex):Q(1)};hh=fixture_h(v,[1,1,2][sec])
    for h,w in queries:
     if h//28==sec:
      e=[0]*8;e[w]=1;c=hh[h%28].get(tuple(e),Q(0))
      if c:HS[h,w][8*sec+row,mi]=c
 h0={}
 for h,w in queries:
  sec,pair=divmod(h,28);i,j=pairs[pair];h0[h]=[mxp1(p,8*sec+7,i)if j==7 else(mxp1(p,8*sec+j,i)-mxp1(p,8*sec+i,j))/2 for p in xp]
 hi=[[Q(c)for c in row]for row in D['first']['all54_actual_pivot_h0_readback']]
 C=[i for i in range(44)if i not in D['m0']['selected_source_columns24']];I=[C[j]for j in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in row]for row in D['minimal']['fixed_minor_data']['inverse']];n3=[[Q(c)for c in row]for row in D['lower2']['all140_native_row20_third_lower_surrogates']];m2=list(itertools.combinations_with_replacement(range(8),2))
 projected=[];source_input=[]
 for p,src in enumerate(PF['all35_constant_original780_source_combinations']):
  bp=[{}for _ in range(48)];pp=[{}for _ in range(8)];nc=[];x1=[[Q(0)]*8 for _ in range(44)]
  for j,c in src:
   c=Q(c)
   for row,v in raw[j].items():
    for var in range(8):x1[row][var]+=c*v.get(1+var,Q(0))
   if j<640:
    for row,node in Bcols[j]:bp[row]=plus(bp[row],tv[node],c)
    for row,node in enumerate(Pcols[j]):
     if node:pp[row]=plus(pp[row],cv[rn[node][1]],c)
   else:nc.append((j-640,c))
  assert all(not v.get(0,Q(0))for v in bp)
  bc=[]
  for i in range(12):
   v=[Q(0)]*100
   for mask in range(4):
    for mi,c in bp[4*i+mask].items():v[4*mi+mask]=c
   bc.append(tuple(v))
  x2=[{}for _ in range(44)]
  for gen in range(9):
   v=a.add(*(a.mul(AA[i][gen],bc[i])for i in range(12)if bc[i]!=a.Z))
   for mask in range(4):x2[4*gen+mask]={mi:v[4*mi+mask]for mi in range(25)if v[4*mi+mask]}
  for row in range(8):x2[36+row]=dict(pp[row])
  for j,c in nc:
   for row,v in native[j].items():x2[row]=plus(x2[row],v,c)
  for row in range(44):
   for var in [0,1,4,5,6,7]:assert x2[row].get(MI[(var,)],Q(0))==x1[row][var]
  a1=[[sum((row[k]*x1[I[k]][v]for k in range(18)),Q(0))for v in range(8)]for row in Ji]
  def prod(v,w,mon):i,j=mon;return v[i]*w[j]+(v[j]*w[i]if i!=j else Q(0))
  a2cache={}
  def a2(mon):
   if mon not in a2cache:a2cache[mon]=mv(Ji,[x2[row].get(MI[mon],Q(0))-sum((prod([xp1(z,row,v)for v in range(8)],a1[k],mon)for k,z in enumerate(xp)),Q(0))for row in I])
   return a2cache[mon]
  x3=[]
  for n,i in enumerate([0,4,5]):
   z=sum((AF1[row][i+1]*bp[row].get(MI[(7,7)],Q(0))+AF1[row][8]*bp[row].get(MI[(i,7)],Q(0))+AF[row][2]*bp[row].get(MI[(i,)],Q(0))+AF[row][[7,9,10][n]]*bp[row].get(MI[(7,)],Q(0))for row in range(48)),Q(0));z+=sum((c*n3[j][n]for j,c in nc),Q(0));x3.append(z)
  vcache={}
  def v2(row,mi):
   if(row,mi)in vcache:return vcache[row,mi]
   mon=MON[mi];assert len(mon)==2;z=sum((m0[row][j]*x2[j].get(mi,Q(0))for j in range(44)),Q(0));sec,i=divmod(row,8)
   if i==7:
    if 7 in mon:
     other=mon[0]
     z+=sum((mu1[j]*x1[j][other]for j in range(44)),Q(0))if row==7 else sum((row15[j]*x1[j][other]for j in range(44)),Q(0))if row==15 else dotfactor(row23basis,p,other)
   else:
    i,j=mon;z+=sum((m1(row,k,i)*x1[k][j]+(m1(row,k,j)*x1[k][i]if i!=j else Q(0))for k in range(44)),Q(0))
   vcache[row,mi]=z;return z
  bh={(h,w):sum((c*v2(row,mi)for(row,mi),c in HS[h,w].items()),Q(0))-sum((h0[h][k]*a1[k][w]for k in range(18)),Q(0))for h,w in queries}
  values=[]
  for f,var,pi in projections:
   val=Q(0)
   for index,c in pi.items():
    target,q=divmod(index,36);mon=m2[q]
    if target<2:
     row=[42,43][target];aa=a2(mon);v=x2[row].get(MI[mon],Q(0))-sum((xp0(z,row)*aa[k]+prod([xp1(z,row,w)for w in range(8)],a1[k],mon)for k,z in enumerate(xp)),Q(0))
    else:
     assert mon==(7,7);n=target-2;i=[0,4,5][n];aa=a2((7,7));v=Q(5,2)*x3[n]+sum((mu1[j]*x2[j].get(MI[(i,7)],Q(0))for j in range(44)),Q(0))+dotfactor(mu2basis,p,i)-sum((hi[n][k]*aa[k]+hp1[n][k]*a1[k][7]for k in range(18)),Q(0))
    val+=c*v
   val-=sum((c*bh[h,w]for h,w,c in eff[f,var]),Q(0));values.append(val)
  projected.append(values);source_input.append({'source_relation':PF['selected35_actual_first_relations'][p],'actual_scalar_X3_zir2':x3,'mu2_X1_surrogate8':[dotfactor(mu2basis,p,v)for v in range(8)],'last23_X1_surrogate8':[dotfactor(row23basis,p,v)for v in range(8)]})
 save('phase1-all35-nine-convolution-projections.json',{'convolution_labels':[[f,v]for f,v,pi in projections],'all35_projected_M2_lower_representatives':projected,'representatives_not_individually_actual_before_UV_cancel':True,'all35_source_input':source_input,'all18_mu2_replacements_recomputed_exact':True,'all54_hXp1_recomputed_exact':True,'actual_gA3':g3,'U_and_V_kept_formal_then_projected_zero':'both independent source sensitivities all106 zero in pinned preflight; lower representatives setU=V0 only for evaluation of the quotient'},live=[tv,cv])
 lookup={(f,v):j for j,(f,v,pi)in enumerate(projections)};actual=[]
 for syz in D['cubic']['all106_first_linear_syzygies280_sparse_columns']:
  v=[Q(0)]*3
  for index,c in syz:
   p,var=divmod(index,8);c=Q(c)
   for f in range(3):
    if(f,var)in lookup:v[f]+=c*projected[p][lookup[f,var]]
  actual.append(v)
 cp,rr=pivots(actual);rank=len(cp);iv,det=invq([[actual[j][i]for j in cp]for i in rr])if rank else([],Q(1));fac=[]
 for col in actual:
  q=mv(iv,[col[i]for i in rr])if rank else[];assert all(sum((actual[cp[k]][i]*q[k]for k in range(rank)),Q(0))==col[i]for i in range(3));fac.append(q)
 assert rank<=1
 exactHF=[5,5,3,2,1]if rank else None
 save('phase2-actual106-cubic-values-and-rank.json',{'all106_actual_cubic_columns3':actual,'actual_right106_rank':rank,'fixed_rank_columns':cp,'fixed_rank_rows':rr,'fixed_minor_determinant':det,'fixed_minor_inverse':iv,'all106_rank_factorizations':fac,'actual_complete917_rank_if_forced':1 if rank else None,'actual_full_HF_if_forced':exactHF,'allzero_is_partial_negative_only':not rank,'new_source_gauge_jet':False,'nilcube_not_decided':True},live=[tv,cv])
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_cache_only_right106_cubic_block','actual_right106_rank':rank,'actual_full_HF_if_forced':exactHF,'all106_source_values_and_factorizations_verified':True,'new_source_gauge_jet_or_R35_r3_callback':False,'nilcube_or_nil_generators_decided':False,'saved_QQ_slots':saved,'peak_saved_plus_live_QQ':peak,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[12]-*.json')},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
