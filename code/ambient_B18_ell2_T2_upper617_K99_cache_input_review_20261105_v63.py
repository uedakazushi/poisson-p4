#!/usr/bin/env python3
"""Fresh P5/native-frame/positive/custody source peer; no signed K responses."""
from fractions import Fraction as Q
from pathlib import Path
import json,math,resource,signal,subprocess,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_CACHE_INPUT_GATE_V2_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CACHE_INPUT_GATE_V2_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='d3bf63533e2ab2a2882a51d09639ccba98d9d966d00edffd73a3f80e54c26bc1'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'pruned99-Kfull-P5-native23-and-prior396-input.json');assert count(d)==4458
 source=rd('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63/actual99-selected-source-K-callback-roots-and-literal-ancestry.json');order=rd('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63/all99-old-plane-order-pruning-and-nativeS0-source-gate.json')
 for p,h in [('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63_REVIEW/report.json','3704c0773a0170d258d4ae4f94f22365547a2bb0c924b3a12ae6a44595576573'),('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63_REVIEW/report.json','9a664530e5075d9cc6dc4f97dab7ff17fba198faec78dfe22878535dadbc9d68'),('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW/report.json','76e1aa8842efe0d6d5af78572641c23437ebc902d612c8424876d0058624a48e')]:assert rd(p)['all_passed'] and pins[p]==h
 assert d['selected99_original768_constant_combinations']==source['selected99_original768_constant_combinations'] and d['typed_program50']==source['typed_program'] and d['old80_root_columns']==source['needed_original_old_columns'] and d['Gamma_two']==source['Gamma_two']
 assert d['source_order_skip_seven_products']==sorted(i for i,n,k in order['only_seven_honest_parent_products']) and d['honest327554_opaque_typed_leaf']==94 and d['gB_plane_order3_typed_leaf']==107
 assert d['all99_weighted_nativeS0zero_source_order_peer_pinned']
 cache='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json';raw=fld(cache,['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se'])
 for k in ['native_L0','native_Lr','native_Le']:assert d[k]==raw[k]
 L0,Lr,Le,S0,Sr,Se=[list(map(sp,raw[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se']]
 nids=sorted({j-640 for c in d['selected99_original768_constant_combinations'] for j,x in c if j>=640});assert len(nids)==23
 for k,A in [('native23_Sr',Sr),('native23_Se',Se)]:assert [(j,sp(v)) for j,v in d[k]]==[(j,A[j]) for j in nids]
 s11=fld('reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',['B11_basis_order','actual_S_source_columns128_rows48_sparse_B11'])
 assert s11['B11_basis_order']==['1','r','r2','z0','z1','z4','z5','z0r','z1r','z4r','z5r'];rr={}
 for j,col in enumerate(s11['actual_S_source_columns128_rows48_sparse_B11']):
  layers=[{} for _ in range(3)]
  for i,entries in col:
   for k,x in sp(entries).items():
    if k<3:layers[k][i]=x
  assert layers[0]==S0[j] and layers[1]==Sr[j]
  if j in nids:rr[j]=layers[2]
 assert [(j,sp(v)) for j,v in d['native23_S2rr']]==sorted(rr.items())
 c16=fld('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',['intrinsic_S2rell_columns','actual_L1r_columns','actual_L1ell_columns'])
 assert Lr==list(map(sp,c16['actual_L1r_columns'])) and Le==list(map(sp,c16['actual_L1ell_columns']))
 re={j:sp(c16['intrinsic_S2rell_columns'][j]) for j in nids};assert [(j,sp(v)) for j,v in d['native23_S2rell']]==sorted(re.items())
 C=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json',['Cfull_two_by2_by56'])['Cfull_two_by2_by56'];assert C==d['Cfull_two_by2_by56'];CF=[[[Q(x) for x in row] for row in ch] for ch in C];crows={i for ch in CF for row in ch for i,x in enumerate(row) if x};assert all(i<48 for i in crows)
 pv=fld('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values'];wanted=order['all24_remaining_physical_boundary_ids'];leaf={j:[str(sp(v).get(k,Q(0))) for k in [0,8,9,24,25]] for j,v in pv if j in wanted};assert sorted(leaf.items())==[(j,v) for j,v in d['physical24_plane5_values']]
 oldids=[j for j,v in d['old80_root_columns']];old=fld('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',['source_columns640_rows56_sparse_dual9'])['source_columns640_rows56_sparse_dual9'];oldr={j:{i:sp(v)[8] for i,v in old[j] if i in crows and sp(v).get(8)} for j in oldids}
 oldell=fld('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns'];olde={j:{i:x for i,x in sp(oldell[j]).items() if i in crows} for j in oldids}
 assert [(j,sp(v)) for j,v in d['old80_r_first_columns']]==sorted(oldr.items()) and [(j,sp(v)) for j,v in d['old80_ell_first_columns']]==sorted(olde.items())
 # Independent nonnegative rational envelopes, with explicit SUM at P/Q.
 maxbits=0
 def mark(a):
  nonlocal maxbits
  maxbits=max(maxbits,a[0].bit_length(),a[1].bit_length());assert maxbits<=147;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 Z=(1,0)
 def bdot(row,col):
  q=Z
  for i,x in row.items():q=plus(q,mul(ex(x),col.get(i,Z)))
  return q
 G=list(map(sp,d['Gamma_two']));R=[];P=[];E=[]
 for t0 in range(2):
  cr={i:x for i,x in enumerate(CF[t0][0]) if x};ce={i:x for i,x in enumerate(CF[t0][1]) if x}
  R.append([bdot(G[t0],{i:ex(x) for i,x in c.items()}) for c in L0])
  P.append([plus(bdot(G[t0],{i:ex(x) for i,x in a.items()}),bdot(cr,{i:ex(x) for i,x in b.items()})) for a,b in zip(Lr,L0)])
  E.append([plus(bdot(G[t0],{i:ex(x) for i,x in a.items()}),bdot(ce,{i:ex(x) for i,x in b.items()})) for a,b in zip(Le,L0)])
 mon=[(0,0),(1,0),(0,1),(2,0),(1,1)];ix={v:i for i,v in enumerate(mon)};jets={};skip=set(d['source_order_skip_seven_products'])
 for i,n in d['typed_program50']:
  if i==94:jets[i]=None
  elif i==107 or i in skip:jets[i]=[Z]*5
  elif n[0]=='actual_C0_scalar_function':assert n[1] in leaf;jets[i]=list(map(ex,leaf[n[1]]))
  elif n[0]=='sum':
   q=[Z]*5
   for k,a,denom in n[1]:
    assert jets[k] is not None
    for u,x in enumerate(jets[k]):q[u]=plus(q[u],mul(ex(Q(a,denom)),x))
   jets[i]=q
  elif n[0]=='product':
   q=[ex(Q(n[1],n[2])),Z,Z,Z,Z]
   for k in n[3]:
    assert jets[k] is not None;v=[Z]*5
    for u,x in enumerate(q):
     for w,y in enumerate(jets[k]):
      m=(mon[u][0]+mon[w][0],mon[u][1]+mon[w][1])
      if m in ix:v[ix[m]]=plus(v[ix[m]],mul(x,y))
    q=v
   jets[i]=q
  else:raise AssertionError(n)
 oldg={}
 for j,col in d['old80_root_columns']:
  v=[[Z,Z],[Z,Z]]
  for row,(kind,k) in col:
   if kind=='C0':
    if k==0:continue
    assert k in leaf;q=list(map(ex,leaf[k]))
   else:assert kind=='typed' and jets[k] is not None;q=jets[k]
   for t0 in range(2):
    for u in range(2):v[t0][u]=plus(v[t0][u],mul(ex(G[t0].get(row,Q(0))),q[3+u]))
  oldg[j]=v
 def weighted(A,c):
  v={}
  for j,x in c.items():
   for i,y in A[j].items():v[i]=plus(v.get(i,Z),mul(ex(x),ex(y)))
  return v
 def rowdot(row,v):
  q=Z
  for i,x in v.items():q=plus(q,mul(row[i],x))
  return q
 env=[]
 for c in d['selected99_original768_constant_combinations']:
  co={j:Q(x) for j,x in c if j<640};cn={j-640:Q(x) for j,x in c if j>=640};sr,se,srr,sre,ar,ae=[weighted(A,cc) for A,cc in [(Sr,cn),(Se,cn),(rr,cn),(re,cn),(oldr,co),(olde,co)]];ch=[]
  for t0 in range(2):
   q=[Z,Z]
   for j,x in co.items():
    for u in range(2):q[u]=plus(q[u],mul(ex(x),oldg[j][t0][u]))
   cr={i:x for i,x in enumerate(CF[t0][0]) if x};ce={i:x for i,x in enumerate(CF[t0][1]) if x}
   q[0]=plus(q[0],plus(bdot(cr,ar),plus(rowdot(R[t0],srr),rowdot(P[t0],sr))))
   q[1]=plus(q[1],plus(plus(bdot(cr,ae),bdot(ce,ar)),plus(rowdot(R[t0],sre),plus(rowdot(P[t0],se),rowdot(E[t0],sr)))))
   ch.append([[str(a),str(b)] for a,b in q])
  env.append(ch)
 assert env==d['positive396_Kfull_envelopes'] and maxbits==147
 custody='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_PREPARED_CUSTODY_20261105_v63/report.json';cu=rd(custody);assert pins[custody]=='c8c9fc9eb54623a457dc5fd689461d15a342909c2791efe2c5882c613cc7eb6e' and cu['all_passed']
 for path,h in cu['literal_input_sha256'].items():assert sha(ROOT/path)==h
 worker=ROOT/cu['worker'];assert sha(worker)==cu['worker_sha256']=='36494eeaabafd15e9c758d8077529eb0268edc676f5ae8de40fce0bc896d09f0'
 helper=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helper.returncode==0 and hashlib_sha(helper.stdout)==cu['safe_help_stdout_sha256']
 parsed=count([G,CF,L0,Lr,Le,{j:Sr[j] for j in nids},{j:Se[j] for j in nids},rr,re,oldr,olde,{j:list(map(Q,v)) for j,v in leaf.items()}]);live=2*count(d)+parsed+1152+300+400+384+200+4752+1000
 assert live==cu['explicit_live_upper']==20036<r['future_value_numeric_bound']==23302
 assert cu['derived_committed_bits']==147 and cu['derived_temporary_bits']==442 and cu['hard_guards']==[384,1153]
 ledger=rd('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_PREPARED_CUSTODY_20261105_v63/explicit-live-RAW-and-prior-ownership.json')
 assert sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_PREPARED_CUSTODY_20261105_v63/explicit-live-RAW-and-prior-ownership.json')==cu['ownership_ledger_sha256']
 rawset=set(cu['literal_input_sha256'])|{cu['worker'],custody,b+'report.json',b+'pruned99-Kfull-P5-native23-and-prior396-input.json'};rawprior=sum((ROOT/path).stat().st_size for path in rawset)+1000000+2000000;assert rawprior<33554432
 proof={'all_original24physical_P5_and80first_fields_match':True,'all128_B11_S0_Sr_match_current_branch':True,'C16_Lr_Le_and23_intrinsic_rell_match':True,'all99_nativeS0zero_inherited3704':True,'all50AST_and_explicit_C0_constant_or_physical_routes':True,'honest94_opaqueNone_with_only7_source_order_skips':True,'truegB107_order3_notrawf63':True,'R_P_Q_semantics':'R=GammaL0;P=GammaLr+C_rL0;Q=GammaLe+C_ellL0. AllbranchcontributionsareSUMS. Whole99 rr/rell effectiveformulas only; no individualpivot substitution.','all396_positive_envelopes_reconstructed_before_signedK':True,'committed147_temp442':True,'explicit_worker_live':live,'numeric_cap23302':True,'whole_RAW_independent_upper':rawprior,'custody_whole_RAW3168312':cu['whole_RAW_upper'],'worker_nonehelpers_safehelp_and_literal_inputs':True,'oldV1_before_signedK_failure_bytes_preserved':True,'no_K_secondresponse_b1_rho_Z2_upper_or_higher_values_evaluated':True}
 ep='independent-K396-source-frame-positive-custody-and-whole-bound.json';put(ep,proof)
 owned=count([r,d,source,order,raw,L0,Lr,Le,S0,Sr,Se,s11,rr,c16,re,C,CF,pv,leaf,old,oldr,oldell,olde,G,R,P,E,jets,oldg,env,cu,ledger,proof,pins])+18000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_code_custody_positive_only_no_signedK':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
def hashlib_sha(b):
 import hashlib
 return hashlib.sha256(b).hexdigest()
if __name__=='__main__':main()
