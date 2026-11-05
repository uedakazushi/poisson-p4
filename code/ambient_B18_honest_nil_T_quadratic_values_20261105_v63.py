#!/usr/bin/env python3
"""Complete6x661 lower quadratic values, ONLY fixed cached source coefficients.
No new physical/gauge/source3 coefficient; all quotient second deviations cancel.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_honest_nil_T_gauge120_cached_20261105_v63 as G
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,polynomial,pivots,invq,mv
from ambient_B18_honest_nil_T_operator_20261104_v63 import fixture_h
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63'
PATHS={
 'gauge_report':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63/report.json',
 'physics':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63/phase1-common105-and-physical-source-base25.json',
 'proper':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63/phase3-actual-proper-G8-base25.json',
 'honest_report':'reports/AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63/report.json',
 'honest':'reports/AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63/phase3-full661-honest-second-sensitivities.json',
 'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'adapter2':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'oldcache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/phase3-oldB-common-source2-values-and-symbols.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'GL':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'first2':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'first3':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'carry':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/actual81-pivot-first-carries-and-six-H1-source-adjoints.json',
 'higher':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63/actual661-first-corrections-m2-replacement-and-native-highest-sensitivities.json',
 'free43':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_FREE43_CACHED_20261105_v63/actual-one-by661-free43-quadratic-values-and-source-sensitivities.json',
 'operator':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63/actual-complete-six-by661-quadratic-source-operator.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(d,k,c):
 if c:
  d[k]=d.get(k,Q(0))+c
  if not d[k]:del d[k]
def plus(d,v,c=Q(1)):
 d=dict(d)
 for k,a in v.items():put(d,k,c*a)
 return d
def sparse(d):return [[i,c]for i,c in sorted(d.items())if c]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0;peak=0
 def save(name,data,live=()):
  nonlocal saved,peak
  q=nqq(data);peak=max(peak,saved+q+nqq(live));assert peak<150000,('AGGREGATE_QQ_GUARD',name,peak)
  blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())+sum(p.stat().st_size for p in out.iterdir())<32<<20
  (out/name).write_text(blob);saved+=q;(out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ':saved,'saved_plus_live_peak':peak,'seconds':time.monotonic()-t})+'\n')
 D={k:read(ROOT/p)for k,p in PATHS.items()};assert D['gauge_report']['all_passed']and D['honest_report']['all661_raw_Bold2_and_P2_honest_second_sensitivities_zero']
 save('preflight.json',{'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'hard_QQ':150000,'scalar_source2_only':True,'no_new_source_gauge_jet':True})
 MON=[tuple(m)for m in D['physics']['base25_monomials']];MI={m:i for i,m in enumerate(MON)};G.MON=MON;G.IDX=MI
 def mul(p,q):return G.smul(p,q)
 def inv(p):
  c=p.get(0,Q(0));assert c;h={i:v for i,v in p.items()if i};return plus(plus({0:1/c},h,-1/c**2),mul(h,h),1/c**3)
 def compose(p,args):
  z={}
  for mon,c in zip(D['source2']['source2_monomials'],p):
   if Q(c):
    v={0:Q(c)}
    for i in mon:v=mul(v,args[i])
    z=plus(z,v)
  return z
 vars8=[{MI[(i,)]:Q(1)}if(i,)in MI else{}for i in range(8)];YA=vars8[:7]+[vars8[7],vars8[7]];YB=[]
 for lin,quad in zip(D['adapter2']['YB_in_YA_linear'],D['adapter2']['YB_in_YA_quadratic']):
  v={}
  for i,c in enumerate(lin):v=plus(v,YA[i],Q(c))
  for mon,c in quad:
   w={0:Q(c)}
   for i,e in enumerate(mon):
    for _ in range(e):w=mul(w,YA[i])
   v=plus(v,w)
  YB.append(v)
 bind={i:{j:Q(c)for j,c in enumerate(v)if Q(c)}for i,v in D['physics']['all105_actual_source2_values']}
 cbcopy=dict(D['copy']['explicit57_source_copy_bindings']);args=[{}]*7+YB[7:9]
 for i,p in D['source2']['both57_source2']['CB']['all57_source2jets']:bind[cbcopy[i]]=compose(p,args)
 def first_poly(p,full9=False):
  p=list(map(Q,p));v={0:p[0]}if p[0]else{}
  for i in range(8):
   if(i,)in MI:
    c=p[i+1]if i<7 or not full9 else p[8]+p[9]
    if c:v[MI[(i,)]]=c
  return v
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in D['source2'][key].items():bind[int(i)]=first_poly(p,True)
 for row in D['atomfirst']['all35_ordered_CA_prefix_firsts']:bind[row['node']]={i:Q(c)for i,c in row['first9_bar']if i==0 or i-1<8 and(i-1,)in MI};bind[row['node']]=first_poly([dict(row['first9_bar']).get(i,'0')for i in range(9)])
 for i,p in D['atomfirst']['all8_g_quotient_overrides']:bind[i]=first_poly([dict(p).get(j,'0')for j in range(9)])
 for i,node in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):bind[node]=vars8[i]if i<7 else{}
 R=D['radical'];rn=R['new_Cscalar_and_R0component_nodes'];Pcols=R['actual_P8_by640_columns'];tn=D['chain']['operator_scalar_nodes'];Bcols=D['chain']['B48_by640_sparse_columns'];broots={i for col in Bcols for row,i in col};proots={rn[i][1]for col in Pcols for i in col if i};tc=set();todo=list(broots);cr=set(proots)
 while todo:
  i=todo.pop()
  if i in tc:continue
  tc.add(i);n=tn[i];o=n[0]
  if o=='actual_C0_scalar_function':cr.add(n[1])
  elif o=='sum':todo.extend(j for j,a,b in n[1])
  elif o=='product':todo.extend(n[3])
  elif o in ['a','m','n','i']:todo.extend(n[1:])
  else:assert o=='c'
 cn=D['C0']['DAG_scalar_nodes'];extra=D['copy']['all57_copied_equation_indices'];cc=set();todo=list(cr)+extra;missing=[]
 while todo:
  i=todo.pop()
  if i in cc:continue
  cc.add(i)
  if i in bind:continue
  n=cn[i];o=n[0]
  if o in ['a','m']:todo.extend(n[1:3])
  elif o in ['n','i']:todo.append(n[1])
  elif o!='c':missing.append([i,n])
 assert not missing,('UNAVAILABLE_SOURCE_CALLBACK',missing)
 uses=collections.Counter()
 for i in cc:
  if i not in bind and cn[i][0]in ['a','m','n','i']:uses.update(cn[i][1:])
 uses.update(list(cr)+extra);cv={};live=0;cvpeak=0
 for i in sorted(cc):
  n=cn[i];o=n[0]
  if i in bind:v=bind[i]
  elif o=='c':v={0:Q(n[1],n[2])}if n[1]else{}
  elif o=='a':v=plus(cv[n[1]],cv[n[2]])
  elif o=='n':v={j:-c for j,c in cv[n[1]].items()}
  elif o=='m':v=mul(cv[n[1]],cv[n[2]])
  else:assert o=='i';v=inv(cv[n[1]])
  cv[i]=v;live+=len(v);cvpeak=max(cvpeak,live);assert cvpeak+saved<150000
  if i not in bind and o in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:live-=len(cv[j]);del cv[j]
 assert all(not cv[i]for i in extra)
 tv={}
 for i in sorted(tc):
  n=tn[i];o=n[0]
  if o=='c':v={0:Q(n[1],n[2])}if n[1]else{}
  elif o=='actual_C0_scalar_function':v=cv[n[1]]
  elif o=='sum':
   v={}
   for j,a,b in n[1]:v=plus(v,tv[j],Q(a,b))
  elif o=='product':
   v={0:Q(n[1],n[2])}
   for j in n[3]:v=mul(v,tv[j])
  elif o=='a':v=plus(tv[n[1]],tv[n[2]])
  elif o=='n':v={j:-c for j,c in tv[n[1]].items()}
  elif o=='m':v=mul(tv[n[1]],tv[n[2]])
  else:assert o=='i';v=inv(tv[n[1]])
  tv[i]=v
 overlap_counts={}
 for key,vals in [('C0_composite_values_and_opaque_second_gradients',cv),('typedB_root_values_and_opaque_second_gradients',tv)]:
  overlap_counts[key]=0
  for i,v,g in D['oldcache'][key]:
   if i not in vals:continue
   assert {j:Q(c)for j,c in v}=={j:c for j,c in vals[i].items()if j<20},('SOURCE2_CACHE_OVERLAP',key,i)
   overlap_counts[key]+=1
 save('phase1-all-source-root-second-representatives.json',{'base25_monomials':MON,'all1239_Broot_second_representatives':[[i,sparse(tv[i])]for i in sorted(broots)],'all768_C0root_second_representatives':[[i,sparse(cv[i])]for i in sorted(cr)],'all57_copy_equations_exact':True,'all_existing_required_B20_values_match':True,'old_cache_overlap_counts':overlap_counts,'source_stream_peak_QQ':cvpeak,'raw_representatives_not_individual_actual_honest_seconds':True},live=bind)
 # Actual proper2 coefficients and native source row20 lower-third surrogate.
 nus=[[Q(c)for c in v]for v in D['physics']['nu_A_B_base25']];a=G.Finite(2,nus);P=D['proper'];core=[tuple(map(Q,P['kappa9_C100'][i]))for i in [0,1,2,6]];vs=[tuple(map(Q,v))for v in P['actualG8_C100'][:5]]+[polynomial(a,D['GL']['forward_GL4'][s],core)for s in ['h','c','r','l']]
 def unit(alg,k):v=[Q(0)]*alg.n;v[k]=Q(1);return tuple(v)
 native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  pols={row:polynomial(a,p,vs)for row,p in col}
  for mask in range(4):
   z={}
   for row,v in pols.items():
    q=a.mul(v,unit(a,mask))
    for k in range(4):
     p={i:q[4*i+k]for i in range(25)if q[4*i+k]}
     if p:z[4*row+k]=p
   native.append(z)
 qv=[];d=polynomial(a,R['d_polynomial'],vs)
 for pol in R['R9_generators_polynomials9']:
  q=a.mul(d,polynomial(a,pol,vs))
  for mask in range(4):qv.append(a.mul(q,unit(a,mask)))
 for i in range(8):
  p=[Q(0)]*100
  if i<7 and(i,)in MI:p[4*MI[(i,)]+1]=Q(-1)
  qv.append(tuple(p))
 gamma=lambda p:[tuple(p[4*i+1]for i in range(25)),tuple(p[4*i+3]for i in range(25)),tuple(p[4*i+2]-Q(11,12)*p[4*i+1]for i in range(25))]
 vg=[gamma(q)for q in qv];m0=[[Q(c)for c in row]for row in D['m0']['m0_24_by44']]
 for j,v in enumerate(vg):
  for sec in range(3):
   for i in range(7):
    if(i,)in MI:assert v[sec][MI[(i,)]]==m0[8*sec+i][j]
   assert v[sec][MI[(7,)]if sec<2 else MI[(7,7)]]==m0[8*sec+7][j]
 mu1=[Q(c)for c in D['first2']['true_actual_m1_row7_r_coefficients44']];assert [v[0][MI[(7,7)]]for v in vg]==mu1
 triples=[(i,7,7)for i in [0,4,5]];MON28=MON+triples;G.MON=MON28;G.IDX={m:i for i,m in enumerate(MON28)}
 a3=G.Finite(2,[v+[Q(0)]*3 for v in nus])
 def lift(v):return tuple(map(Q,v))+(Q(0),)*12
 vs3=[lift(v)for v in vs];n3=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  p=next((p for row,p in col if row==5),[]);q=polynomial(a3,p,vs3)
  for mask in range(4):
   v=a3.mul(q,unit(a3,mask));n3.append([v[4*(25+i)]for i in range(3)])
 G.MON=MON;G.IDX=MI
 # Actual Res(A0+A1); omitted second Cscalar coefficients multiply B0=0.
 A1=D['X1']['new_A12_by9_C36'];AA=[]
 for row in A1:
  rr=[]
  for v in row:
   w=[Q(0)]*100
   for k in range(4):w[k]=Q(v[k])
   for variable in range(8):
    if(variable,)in MI:
     for k in range(4):w[4*MI[(variable,)]+k]=Q(v[4*(variable+1)+k])
   rr.append(tuple(w))
  AA.append(rr)
 raw1=[]
 for col in D['X1']['all780_source_columns44_sparse_dual9']:raw1.append({row:dict((k,Q(c))for k,c in sp)for row,sp in col})
 xp=D['X1']['X0_pivot_columns'];C=[i for i in range(44)if i not in D['m0']['selected_source_columns24']];I=[C[i]for i in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in row]for row in D['minimal']['fixed_minor_data']['inverse']]
 L=[[(h,v,Q(c))for h,v,c in row]for row in D['carry']['all6_actual_effective_H1_callbacks_row84_var8_sparse']];queries=sorted({(h,v)for row in L for h,v,c in row});pairs2=list(itertools.combinations(range(8),2));qm=[tuple(m)for m in D['quad']['quadratic_monomials36']];Pi=[{i:Q(c)for i,c in row}for row in D['quad']['exact_quotient_duals_rows_sparse']]
 # Canonical H2 linear stencils only for the52 actual needed carry callbacks.
 HS={q:{}for q in queries}
 for sec in range(3):
  for row in range(8):
   for mi,mon in enumerate(MON):
    if len(mon)!=2 or row==7 and mon==(7,7):continue
    e=[0]*8
    for k in mon:e[k]+=1
    v=[{}for _ in range(8)];v[row]={tuple(e):Q(1)};hh=fixture_h(v,[1,1,2][sec])
    for h,var in queries:
     if h//28==sec:
      ee=[0]*8;ee[var]=1;c=hh[h%28].get(tuple(ee),Q(0))
      if c:HS[h,var][sec*8+row,mi]=c
 def m1(row,j,var):
  sec,i=divmod(row,8)
  if i==7:raise AssertionError('LAST_M1_VALUE_MUST_BE_REPLACED_BY_SOURCE_CONTRACTION')
  if var!=7 and var>i:return Q(0)
  mon=tuple(sorted((i,var)));assert mon in MI,('UNAVAILABLE_M1_PAIR',row,j,var,mon)
  return vg[j][sec][MI[mon]]
 def xpfirst(p,row,var):return raw1[p].get(row,{}).get(1+var,Q(0))
 def xpzero(p,row):return raw1[p].get(row,{}).get(0,Q(0))
 def mxp1(p,row,var):
  z=sum((m0[row][j]*xpfirst(p,j,var)for j in range(44)),Q(0))
  if row%8==7:
   if var==7:return Q(0) # exact pure-r source cycle
  else:z+=sum((m1(row,j,var)*xpzero(p,j)for j in range(44)),Q(0))
  return z
 h0={}
 for h,var in queries:
  sec,pairid=divmod(h,28);i,j=pairs2[pairid]
  h0[h]=[mxp1(p,8*sec+7,i)if j==7 else(mxp1(p,8*sec+j,i)-mxp1(p,8*sec+i,j))/2 for p in xp]
 hi=[[Q(c)for c in row]for row in D['first3']['all54_actual_pivot_h0_readback']]
 for k,i in enumerate([0,4,5]):
  h=pairs2.index((i,7));assert h0.get(h,[mxp1(p,7,i)for p in xp])==hi[k]
 lastcontracts={row:[-sum((m0[row][j]*xpfirst(p,j,7)for j in range(44)),Q(0))for p in xp]for row in [7,15,23]}
 # Exact individual Xp row20 second callbacks need no honest second values.
 AF=[[Q(c)for c in row]for row in D['carry']['all48_actual_ResAT_row20_B11']];AF1=[[Q(c)for c in row]for row in D['first2']['all48_actual_adapter_row20_first9']]
 def old_scalar2(p,mon):
  z=Q(0)
  for row,node in Bcols[p]:
   q=tv[node]
   if mon==(7,7):z+=AF[row][2]*q.get(0,Q(0))+AF1[row][8]*q.get(MI[(7,)],Q(0))
   else:
    i=mon[0];bi={0:7,4:9,5:10}[i];z+=AF[row][bi]*q.get(0,Q(0))+AF1[row][1+i]*q.get(MI[(7,)],Q(0))+AF1[row][8]*q.get(MI[(i,)],Q(0))
  return z
 def xpscalar2(p,mon):return old_scalar2(p,mon)if p<640 else native[p-640].get(20,{}).get(MI[mon],Q(0))
 hp1=[[Q(5,2)*xpscalar2(p,(i,7))+sum((mu1[j]*xpfirst(p,j,i)for j in range(44)),Q(0))for p in xp]for i in [0,4,5]]
 mu2p=[-sum((mu1[j]*xpfirst(p,j,7)for j in range(44)),Q(0))-Q(5,2)*xpscalar2(p,(7,7))for p in xp]
 save('phase2-actual-native-and-m1-source-inputs.json',{'all140_native_row20_third_lower_surrogates':n3,'raw_native_third_values_not_claimed_individually_actual':True,'actual_m1_row7_r':mu1,'all_needed_hXp0':[[h,v]for h,v in sorted(h0.items())],'all54_hXp1_r':hp1,'all18_mu2_Xp0_replacement_values':mu2p,'all54_last_m1_Xp0_contractions':[[row,v]for row,v in lastcontracts.items()],'all_constant_H2_stencils':[[h,var,[[row,mi,c]for(row,mi),c in sp.items()]]for(h,var),sp in HS.items()],'all_actual_m0_readbacks_verified':True,'no_last_m1_or_m2_higher_value_used':True},live=[tv,cv])
 allvalues=[];alldirect=[];allcarry=[];allx3=[];sources=D['free43']['all661_constant_original780_source_combinations'];a1s=D['higher']['all661_a1_18_by8_sparse'];original=D['operator']['all661_source_columns_in_original849'];d3=D['Koszul']['all8_differentials_sparse_columns'][2]
 for sc,source in enumerate(sources):
  bp=[{}for _ in range(48)];pp=[{}for _ in range(8)];nc=[];x1=[[Q(0)]*8 for _ in range(44)]
  for j,c in source:
   c=Q(c)
   for row,v in raw1[j].items():
    for var in range(8):x1[row][var]+=c*v.get(1+var,Q(0))
   if j<640:
    for row,node in Bcols[j]:bp[row]=plus(bp[row],tv[node],c)
    for row,node in enumerate(Pcols[j]):
     if node:pp[row]=plus(pp[row],cv[rn[node][1]],c)
   else:nc.append((j-640,c))
  assert all(not q.get(0,Q(0))for q in bp),('B0_SOURCE_NOT_ZERO',sc)
  bc=[]
  for i in range(12):
   q=[Q(0)]*100
   for mask in range(4):
    for mi,c in bp[4*i+mask].items():q[4*mi+mask]=c
   bc.append(tuple(q))
  x2=[{}for _ in range(44)]
  for gen in range(9):
   v=a.add(*(a.mul(AA[i][gen],bc[i])for i in range(12)if bc[i]!=a.Z))
   for mask in range(4):x2[4*gen+mask]={mi:v[4*mi+mask]for mi in range(25)if v[4*mi+mask]}
  for i in range(8):x2[36+i]=dict(pp[i])
  for j,c in nc:
   for row,p in native[j].items():x2[row]=plus(x2[row],p,c)
  # Same actual X1 must be recovered in every active direction.
  for row in range(44):
   for var in [0,1,4,5,6,7]:assert x2[row].get(MI[(var,)],Q(0))==x1[row][var],('ACTUAL_FIRST_SOURCE_READBACK',sc,row,var)
  a1=[[Q(0)]*8 for _ in range(18)]
  for index,c in a1s[sc]:a1[index//8][index%8]=Q(c)
  for row in range(44):
   for var in range(8):assert x1[row][var]==sum((xpzero(p,row)*a1[k][var]for k,p in enumerate(xp)),Q(0))
  def prod11(v,w,mon):
   i,j=mon;return v[i]*w[j]+(v[j]*w[i]if i!=j else Q(0))
  def a2(mon):
   z=[x2[row].get(MI[mon],Q(0))-sum((prod11([xpfirst(p,row,v)for v in range(8)],a1[k],mon)for k,p in enumerate(xp)),Q(0))for row in I]
   return mv(Ji,z)
  a2cache={}
  def get_a2(mon):
   if mon not in a2cache:a2cache[mon]=a2(mon)
   return a2cache[mon]
  row3=[]
  for n,i in enumerate([0,4,5]):
   z=sum((AF1[row][i+1]*bp[row].get(MI[(7,7)],Q(0))+AF1[row][8]*bp[row].get(MI[(i,7)],Q(0))+AF[row][2]*bp[row].get(MI[(i,)],Q(0))+AF[row][[7,9,10][n]]*bp[row].get(MI[(7,)],Q(0))for row in range(48)),Q(0))
   z+=sum((c*n3[j][n]for j,c in nc),Q(0));row3.append(z)
  direct=[Q(0)]*6
  for f in range(6):
   for index,c in Pi[f].items():
    target,mi=divmod(index,36);mon=qm[mi]
    if target<2:
     row=[42,43][target];aa=get_a2(mon);v=x2[row].get(MI[mon],Q(0))-sum((xpzero(p,row)*aa[k]+prod11([xpfirst(p,row,z)for z in range(8)],a1[k],mon)for k,p in enumerate(xp)),Q(0))
    else:
     assert mon==(7,7);n=target-2;i=[0,4,5][n];aa=get_a2((7,7));v=Q(5,2)*row3[n]+sum((mu1[j]*x2[j].get(MI[(i,7)],Q(0))for j in range(44)),Q(0))+sum((mu2p[k]*a1[k][i]for k in range(18)),Q(0))-sum((hi[n][k]*aa[k]+hp1[n][k]*a1[k][7]for k in range(18)),Q(0))
    direct[f]+=c*v
  cachev={}
  def v2(row,mi):
   key=row,mi
   if key in cachev:return cachev[key]
   mon=MON[mi];assert len(mon)==2
   z=sum((m0[row][j]*x2[j].get(mi,Q(0))for j in range(44)),Q(0));sec,i=divmod(row,8)
   if i==7:
    if 7 in mon:
     other=mon[0];z+=sum((lastcontracts[row][k]*a1[k][other]for k in range(18)),Q(0))
   else:
    v,w=mon
    z+=sum((m1(row,j,v)*x1[j][w]+(m1(row,j,w)*x1[j][v]if v!=w else Q(0))for j in range(44)),Q(0))
   cachev[key]=z;return z
  b1={}
  for h,var in queries:
   z=sum((c*v2(row,mi)for(row,mi),c in HS[h,var].items()),Q(0))-sum((h0[h][k]*a1[k][var]for k in range(18)),Q(0))
   for j,c in original[sc]:
    if j>=681:
     for rr,mon,sign in d3[j-681]:
      if rr==h and sum(mon)==1 and mon[var]==1:z+=Q(c)*Q(sign)
   b1[h,var]=z
  for row in [7,15,23]:assert v2(row,MI[(7,7)])==0,('LAST_CYCLE_R2_NOT_ZERO',sc,row)
  carry=[-sum((c*b1[h,var]for h,var,c in row),Q(0))for row in L];final=[v+w for v,w in zip(direct,carry)]
  assert final[3]==Q(D['free43']['all661_actual_quadratic_dual_values'][sc]),('FREE43_INDEPENDENT_ROW',sc)
  allvalues.append(final);alldirect.append(direct);allcarry.append(carry);allx3.append(row3)
 cp,cr=pivots(allvalues);rank=len(cp);iv,det=invq([[allvalues[j][i]for j in cp]for i in cr])if rank else([],Q(1));factor=[]
 for col in allvalues:
  q=mv(iv,[col[i]for i in cr])if rank else[];assert all(sum((allvalues[cp[k]][i]*q[k]for k in range(rank)),Q(0))==col[i]for i in range(6));factor.append(q)
 assert 1<=rank<=4
 save('phase3-all661-actual-six-quadratic-values-and-rank.json',{'all661_direct_six_values':alldirect,'all661_actual_carry_six_values':allcarry,'all661_actual_quadratic_six_values':allvalues,'all661_row20_third_functionals3':allx3,'all661_source_paths_include_15mixed_d3':True,'all661_free43_values_independently_match':True,'all1983_last_component_pure_r2_cycles_verified':True,'actual_operator_rank':rank,'fixed_rank_columns':cp,'fixed_rank_rows':cr,'fixed_minor_determinant':det,'fixed_minor_inverse':iv,'all661_rank_factorizations':factor,'actual_HF2_Tnil':6-rank,'new_source_or_gauge_value':False},live=[tv,cv])
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_complete_quadratic_Loewy_layer_Tnil_from_fixed_caches','actual_complete_operator_shape':[6,661],'actual_rank':rank,'actual_HF_Tnil_through2':[5,5,6-rank],'all_mixed_d3_and_C1_carries_retained':True,'new_physical_source_or_gauge_jet':False,'no_original_nil_generators_or_cube_claim':True,'saved_QQ_slots':saved,'peak_saved_plus_live_QQ':peak,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[123]-*.json')},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
