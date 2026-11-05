#!/usr/bin/env python3
"""Full661 honest-second source sensitivities, including actual carry.
Only origin derivatives with respect to true opaque quotient leaves; no jet solve.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,polynomial
from ambient_B18_honest_nil_T_operator_20261104_v63 import fixture_h
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63'
PATHS={
 'free43':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_FREE43_CACHED_20261105_v63/actual-one-by661-free43-quadratic-values-and-source-sensitivities.json',
 'carry':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/actual81-pivot-first-carries-and-six-H1-source-adjoints.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'first2':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'first3':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'rootcache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/phase3-oldB-common-source2-values-and-symbols.json',
 'sourcecache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'J':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def put(a,k,c):
 if c:
  a[k]=a.get(k,Q(0))+c
  if not a[k]:del a[k]
def plus(a,b,c=Q(1)):
 a=dict(a)
 for k,v in b.items():put(a,k,c*v)
 return a
def sp(a):return [[i,c]for i,c in sorted(a.items())if c]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0
 def save(name,d):
  nonlocal saved
  saved+=nqq(d);assert saved<150000,('VISIBLE_SLOT_GUARD',saved)
  blob=json.dumps(enc(d),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20
  (out/name).write_text(blob);(out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ':saved,'seconds':time.monotonic()-t})+'\n')
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 save('preflight.json',{'task':'origin-gradient extension, full661 sensitivity only','no_actual_unknown_coefficient_evaluated':True,'hard_QQ_saved_plus_live':150000,'inputs':{p:sha(ROOT/p)for p in PATHS.values()}})
 # Actual origins and derivative IDs of previously source-proved coefficient functions.
 stops={}
 for key in ['all105_actual_source2_values','all57_actual_copy_source2_values']:
  for i,v in D['sourcecache'][key]:stops[i]=(dict((j,Q(c))for j,c in v).get(0,Q(0)),{})
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,v in D['J'][key].items():stops[int(i)]=(Q(v[0]),{int(i):Q(1)})
 for row in D['atomfirst']['all35_ordered_CA_prefix_firsts']:
  stops[row['node']]=(dict((j,Q(c))for j,c in row['first9_bar']).get(0,Q(0)),{row['node']:Q(1)})
 for i,v in D['atomfirst']['all8_g_quotient_overrides']:stops[i]=(dict((j,Q(c))for j,c in v).get(0,Q(0)),{i:Q(1)})
 for i,v,g in D['rootcache']['C0_composite_values_and_opaque_second_gradients']:
  stops[i]=(dict((j,Q(c))for j,c in v).get(0,Q(0)),dict((j,Q(c))for j,c in g))
 rn=D['radical']['new_Cscalar_and_R0component_nodes'];Pcols=D['radical']['actual_P8_by640_columns'];proots={i for col in Pcols for i in col if i}
 assert all(rn[i][0]=='R0_original_honest_C0_scalar'for i in proots)
 cnodes=D['C0']['DAG_scalar_nodes'];tnodes=D['chain']['operator_scalar_nodes'];Bcols=D['chain']['B48_by640_sparse_columns'];broots={i for col in Bcols for row,i in col}
 oldtyped={i:(dict((j,Q(c))for j,c in v).get(0,Q(0)),dict((j,Q(c))for j,c in g))for i,v,g in D['rootcache']['typedB_root_values_and_opaque_second_gradients']}
 tc=set();cr={rn[i][1]for i in proots};todo=list(broots)
 while todo:
  i=todo.pop()
  if i in tc or i in oldtyped:continue
  tc.add(i);n=tnodes[i];o=n[0]
  if o=='actual_C0_scalar_function':cr.add(n[1])
  elif o=='sum':todo.extend(j for j,a,b in n[1])
  elif o=='product':todo.extend(n[3])
  elif o in ['a','m','n','i']:todo.extend(n[1:])
  else:assert o=='c',(i,n)
 cc=set();todo=list(cr);missing=[]
 while todo:
  i=todo.pop()
  if i in cc or i in stops:continue
  cc.add(i);n=cnodes[i];o=n[0]
  if o in ['a','m']:todo.extend(n[1:3])
  elif o in ['n','i']:todo.append(n[1])
  elif o!='c':missing.append([i,n])
 save('phase1-full-root-ancestor-profile.json',{'B_roots':len(broots),'P_roots':len(proots),'new_typed_ancestors':len(tc),'C0_roots':len(cr),'new_C0_ancestors':len(cc),'missing_actual_source_leaves':missing})
 assert not missing,('UNAVAILABLE_ORIGIN_OR_QUOTIENT_FIRST',missing)
 val=dict(stops);peak=0
 def add(x,y):return x[0]+y[0],plus(x[1],y[1])
 def neg(x):return -x[0],{k:-v for k,v in x[1].items()}
 def mul(x,y):return x[0]*y[0],plus({k:v*y[0]for k,v in x[1].items()},y[1],x[0])
 def inv(x):
  assert x[0];return 1/x[0],{k:-v/x[0]**2 for k,v in x[1].items()}
 for i in sorted(cc):
  n=cnodes[i];o=n[0]
  val[i]=(Q(n[1],n[2]),{})if o=='c'else add(val[n[1]],val[n[2]])if o=='a'else neg(val[n[1]])if o=='n'else mul(val[n[1]],val[n[2]])if o=='m'else inv(val[n[1]])
  peak=max(peak,sum(1+len(g)for z,g in val.values()));assert peak+saved<150000
 tv=dict(oldtyped)
 for i in sorted(tc):
  n=tnodes[i];o=n[0]
  if o=='c':v=(Q(n[1],n[2]),{})
  elif o=='actual_C0_scalar_function':v=val[n[1]]
  elif o=='sum':
   v=(Q(0),{})
   for j,a,b in n[1]:v=add(v,mul((Q(a,b),{}),tv[j]))
  elif o=='product':
   v=(Q(n[1],n[2]),{})
   for j in n[3]:v=mul(v,tv[j])
  elif o=='a':v=add(tv[n[1]],tv[n[2]])
  elif o=='n':v=neg(tv[n[1]])
  elif o=='m':v=mul(tv[n[1]],tv[n[2]])
  else:assert o=='i';v=inv(tv[n[1]])
  tv[i]=v
 save('phase2-allB-and-P-true-honest-second-gradients.json',{'all1239_B_root_origin_and_gradients':[[i,tv[i][0],sp(tv[i][1])]for i in sorted(broots)],'all723_P_root_origin_and_gradients':[[i,rn[i][1],val[rn[i][1]][0],sp(val[rn[i][1]][1])]for i in sorted(proots)],'peak_origin_gradient_QQ':peak,'new_source_or_jet':False})
 # Exact ordinary coefficient adjoints to deltaX2 and deltaBold2.
 mons=list(itertools.combinations_with_replacement(range(8),2));pi=[]
 for row in D['quad']['exact_quotient_duals_rows_sparse']:pi.append(dict((j,Q(c))for j,c in row))
 M=D['m0'];m0=[[Q(c)for c in row]for row in M['m0_24_by44']];C=[i for i in range(44)if i not in M['selected_source_columns24']];I=[C[j]for j in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in row]for row in D['minimal']['fixed_minor_data']['inverse']]
 xp=D['X1']['X0_pivot_columns'];raw1=[]
 for col in D['X1']['all780_source_columns44_sparse_dual9']:
  raw1.append({row:dict((i,Q(c))for i,c in v)for row,v in col})
 mu1=[Q(c)for c in D['first2']['true_actual_m1_row7_r_coefficients44']];hi=[[Q(c)for c in row]for row in D['first3']['all54_actual_pivot_h0_readback']]
 AF=[[Q(c)for c in row]for row in D['first2']['all48_actual_adapter_row20_first9']];L=[[(h,v,Q(c))for h,v,c in row]for row in D['carry']['all6_actual_effective_H1_callbacks_row84_var8_sparse']]
 pair=list(itertools.combinations(range(8),2));AD=[[[Q(0)]*44 for q in mons]for f in range(6)];AB=[[[Q(0)]*48 for q in mons]for f in range(6)]
 # H2 coefficient stencil, linear extension on the proven last-component pure-r zero subspace.
 H={}
 for sec in range(3):
  for row in range(8):
   for q,mon in enumerate(mons):
    if row==7 and mon==(7,7):continue
    e=[0]*8
    for k in mon:e[k]+=1
    v=[{}for _ in range(8)];v[row]={tuple(e):Q(1)};hh=fixture_h(v,[1,1,2][sec]);z={}
    for h,var,c in sum(L,[]):
     if h//28!=sec:continue
     ee=[0]*8;ee[var]=1;a=hh[h%28].get(tuple(ee),Q(0))
     if a:z[h,var]=a
    H[sec,row,q]=z
 for f in range(6):
  for q,mon in enumerate(mons):
   for j in range(44):
    aI=[Ji[k][I.index(j)]if j in I else Q(0)for k in range(18)]
    v=Q(0)
    for tt,outrow in enumerate([42,43]):
     coef=pi[f].get(36*tt+q,Q(0))
     if coef:v+=coef*(Q(j==outrow)-sum((raw1[p].get(outrow,{}).get(0,Q(0))*aI[k]for k,p in enumerate(xp)),Q(0)))
    for n,zi in enumerate([0,4,5]):
     coef=pi[f].get(36*(n+2)+35,Q(0))
     if mon==(zi,7):v+=coef*mu1[j]
     if mon==(7,7):v-=coef*sum((hi[n][k]*aI[k]for k in range(18)),Q(0))
    for sec in range(3):
     for row in range(8):
      if m0[8*sec+row][j]:v-=m0[8*sec+row][j]*sum((c*H.get((sec,row,q),{}).get((h,var),Q(0))for h,var,c in L[f]),Q(0))
    AD[f][q][j]=v
   for row in range(48):
    for n,zi in enumerate([0,4,5]):
     coef=Q(5,2)*pi[f].get(36*(n+2)+35,Q(0))
     if mon==(7,7):AB[f][q][row]+=coef*AF[row][zi+1]
     if mon==(zi,7):AB[f][q][row]+=coef*AF[row][8]
 # The actual coefficient-ring restriction Res(A0), with both monics zero at origin.
 Avec=D['X1']['new_A12_by9_C36'];RA=[[Q(0)]*48 for _ in range(36)]
 for oldrow in range(48):
  i,mask=divmod(oldrow,4)
  for gen in range(9):
   for k,c in enumerate(Avec[i][gen][:4]):
    if not k&mask:RA[4*gen+(k|mask)][oldrow]+=Q(c)
 Bgrad=[];Dgrad=[]
 for col,pcol in zip(Bcols,Pcols):
  b={};d={}
  for row,node in col:
   for leaf,c in tv[node][1].items():b.setdefault(leaf,{})[row]=b.setdefault(leaf,{}).get(row,Q(0))+c
  for leaf,v in b.items():
   d[leaf]={}
   for j in range(36):
    c=sum((RA[j][r]*a for r,a in v.items()if RA[j][r]),Q(0))
    if c:d[leaf][j]=c
  for k,node in enumerate(pcol):
   if node:
    for leaf,c in val[rn[node][1]][1].items():put(d.setdefault(leaf,{}),36+k,c)
  Bgrad.append(b);Dgrad.append(d)
 gknown={i for i,v in D['g2']['actual_eight_regular_g_second_values']};nonzero=[];needed=set();tested=0;knownresponses=[]
 for sc,source in enumerate(D['free43']['all661_constant_original780_source_combinations']):
  bg={};dg={}
  for j,c in source:
   if j>=640:continue
   c=Q(c)
   for leaf,v in Bgrad[j].items():bg[leaf]=plus(bg.get(leaf,{}),v,c)
   for leaf,v in Dgrad[j].items():dg[leaf]=plus(dg.get(leaf,{}),v,c)
  for leaf in set(bg)|set(dg):
   b=bg.get(leaf,{});d=dg.get(leaf,{})
   if not b and not d:continue
   for q,mon in enumerate(mons):
    ans=[]
    for f in range(6):ans.append(sum((AD[f][q][j]*c for j,c in d.items()),Q(0))+sum((AB[f][q][j]*c for j,c in b.items()),Q(0)))
    tested+=6
    if any(ans):
     rec=[sc,leaf,q,sp(dict(enumerate(ans)))]
     if mon==(7,7)and leaf in gknown:knownresponses.append(rec)
     else:nonzero.append(rec);needed.add((leaf,q))
    assert sum(len(v[3])for v in nonzero)+saved<150000,('SOURCE_SENSITIVITY_GUARD',sc)
 save('phase3-full661-honest-second-sensitivities.json',{'monomial_order36':mons,'adjoint_deltaX2_6_by36_by44_sparse':[[sp(dict(enumerate(v)))for v in row]for row in AD],'adjoint_deltaBold2_6_by36_by48_sparse':[[sp(dict(enumerate(v)))for v in row]for row in AB],'all_remaining_nonzero_formal_sensitivities':nonzero,'all_known_eight_g2_r2_responses':knownresponses,'unknown_callback_union':[[i,q]for i,q in sorted(needed)],'all_formal_sensitivities_tested':tested,'all_unknown_honest_seconds_cancel':not nonzero,'all661_raw_Bold2_and_P2_honest_second_sensitivities_zero':tested==0,'stronger_cancellation_scope':'tested=0 means the complete raw Bold2 and deltaX2 source gradients themselves vanish, before any six-adjoint contraction; deltaP2 then also vanishes','source_cycle_last_component_pure_r2_is_exact_zero':'H2 ignores only this coefficient, certified by d1*m*X=0; no actual coefficient setzero by approximation','m1_m2_values_not_evaluated':True})
 r={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_full661_honest_second_sensitivities','all_unknown_honest_second_sensitivities_cancel':not nonzero,'all661_raw_Bold2_and_P2_honest_second_sensitivities_zero':tested==0,'remaining_nonzero_response_blocks':len(nonzero),'remaining_callback_count':len(needed),'all_original_B_and_P_roots_source_bound':True,'new_source_or_gauge_jet_or_actual_unknown_second_value':False,'saved_QQ_support':saved,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[123]-*.json')},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
