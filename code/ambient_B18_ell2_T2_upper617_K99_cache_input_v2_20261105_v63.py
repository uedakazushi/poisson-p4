#!/usr/bin/env python3
"""Prior-only P5/first/native cache gate for396 Kfull coefficients."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_CACHE_INPUT_GATE_V2_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};maximum_selected_field_numeric=0
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal maximum_selected_field_numeric
  pins[p]=H.sha(ROOT/p);v=H.selected(p,k);maximum_selected_field_numeric=max(maximum_selected_field_numeric,count(v));return v
 cb='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63/actual99-selected-source-K-callback-roots-and-literal-ancestry.json';src=read(cb)
 op='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63/all99-old-plane-order-pruning-and-nativeS0-source-gate.json';order=read(op)
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json';C=field(cp,'Cfull_two_by2_by56');CF=[[[F(x) for x in v] for v in ch] for ch in C]
 assert len(CF)==2 and all(len(ch)==2 and all(len(v)==56 for v in ch) for ch in CF)
 Crows={i for ch in CF for v in ch for i,x in enumerate(v) if x};assert all(i<48 for i in Crows)
 combs=src['selected99_original768_constant_combinations'];oldids={j for j,col in src['needed_original_old_columns']};nids={j-640 for c in combs for j,x in c if j>=640};assert len(oldids)==80 and len(nids)==23
 cache='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 raw={k:field(cache,k) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se']};L0=list(map(H.parse,raw['native_L0']));Lr=list(map(H.parse,raw['native_Lr']));Le=list(map(H.parse,raw['native_Le']));S0=list(map(H.parse,raw['native_S0']));Sr=list(map(H.parse,raw['native_Sr']));Se=list(map(H.parse,raw['native_Se']))
 b11='reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json';assert field(b11,'B11_basis_order')==['1','r','r2','z0','z1','z4','z5','z0r','z1r','z4r','z5r']
 s11=field(b11,'actual_S_source_columns128_rows48_sparse_B11');assert len(s11)==128
 rr={};lower_checks=0
 for j,col in enumerate(s11):
  a={};b={};q={}
  for i,p in col:
   v=dict(p)
   if v.get(0):a[i]=F(v[0])
   if v.get(1):b[i]=F(v[1])
   if j in nids and v.get(2):q[i]=F(v[2])
  assert a==S0[j] and b==Sr[j];lower_checks+=96
  if j in nids:rr[j]=q
 del s11
 c16='reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json'
 reall=field(c16,'intrinsic_S2rell_columns');assert len(reall)==128
 re={j:H.parse(reall[j]) for j in nids};del reall
 assert Lr==list(map(H.parse,field(c16,'actual_L1r_columns'))) and Le==list(map(H.parse,field(c16,'actual_L1ell_columns')))
 physical='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';wanted=set(order['all24_remaining_physical_boundary_ids']);assert len(wanted)==24
 full=field(physical,'all_requested164_root_values');vals={i:H.parse(v) for i,v in full if i in wanted};del full
 assert set(vals)==wanted
 leaves={i:[str(v.get(k,F(0))) for k in [0,8,9,24,25]] for i,v in vals.items()};del vals
 oldr={};olde={}
 cols=field('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9')
 for j in oldids:
  v={}
  for i,p in cols[j]:
   if i in Crows:
    x=H.parse(p).get(8,F(0))
    if x:v[i]=x
  oldr[j]=v
 del cols
 cols=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json','actual_first640_ell_columns')
 for j in oldids:olde[j]={i:x for i,x in H.parse(cols[j]).items() if i in Crows}
 del cols
 peers=[
 'reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json']
 for p in peers:assert read(p)['all_passed']
 # All numerical source fields retained here are existing cached values. The
 # following20k reserve covers positive typed jets/rows/396 output envelopes.
 sourcebase=count([src,order,C,CF,raw,L0,Lr,Le,S0,Sr,Se,rr,re,leaves,oldr,olde,pins]);sourceowned=sourcebase+2*maximum_selected_field_numeric+20000;assert sourceowned<150000
 maximum=0
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());assert maximum<=384;return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 zero=(1,0)
 def bdot(row,col):
  q=zero
  for i,x in row.items():q=add(q,mul(ex(x),col.get(i,zero)))
  return q
 def left(row,cols):return [bdot(row,{i:ex(x) for i,x in col.items()}) for col in cols]
 G=list(map(H.parse,src['Gamma_two']));R=[];P=[];Q=[]
 for t in range(2):
  cr={i:x for i,x in enumerate(CF[t][0]) if x};ce={i:x for i,x in enumerate(CF[t][1]) if x}
  R.append(left(G[t],L0));P.append([add(a,b) for a,b in zip(left(G[t],Lr),left(cr,L0))]);Q.append([add(a,b) for a,b in zip(left(G[t],Le),left(ce,L0))])
 mon=[(0,0),(1,0),(0,1),(2,0),(1,1)];idx={m:i for i,m in enumerate(mon)};table={(i,j):idx[(a[0]+b[0],a[1]+b[1])] for i,a in enumerate(mon) for j,b in enumerate(mon) if (a[0]+b[0],a[1]+b[1]) in idx}
 def jadd(a,b):return [add(x,y) for x,y in zip(a,b)]
 def jscale(a,c):return [mul(x,ex(c)) for x in a]
 def jmul(a,b):
  q=[zero]*5
  for (i,j),k in table.items():q[k]=add(q[k],mul(a[i],b[j]))
  return q
 skip={i for i,n,o in order['only_seven_honest_parent_products']};honest=order['honest_function327554_typed_leaf'];gb=order['gB299650_typed_leaf'];jets={}
 for i,n in src['typed_program']:
  if i==honest:jets[i]=None
  elif i in skip or i==gb:jets[i]=[zero]*5
  elif n[0]=='actual_C0_scalar_function':assert n[1] in leaves;jets[i]=list(map(ex,leaves[n[1]]))
  elif n[0]=='sum':
   v=[zero]*5
   for k,a,b in n[1]:assert jets[k] is not None;v=jadd(v,jscale(jets[k],F(a,b)))
   jets[i]=v
  elif n[0]=='product':
   v=[ex(F(n[1],n[2])),zero,zero,zero,zero]
   for k in n[3]:assert jets[k] is not None;v=jmul(v,jets[k])
   jets[i]=v
  else:raise AssertionError(n)
 oldg={}
 for j,col in src['needed_original_old_columns']:
  q=[[[zero,zero][k] for k in range(2)] for t in range(2)]
  for row,(kind,root) in col:
   if kind=='C0':
    if root==0:continue
    assert root in leaves,'directC0rootrequiresliteralphysical24binding'
    jet=list(map(ex,leaves[root]))
   else:
    assert kind=='typed' and root in jets and jets[root] is not None
    jet=jets[root]
   for t in range(2):
    for k,mi in enumerate([3,4]):q[t][k]=add(q[t][k],mul(ex(G[t].get(row,F(0))),jet[mi]))
  oldg[j]=q
 def weighted(cols,c):
  v={}
  for j,a in c.items():
   for i,x in cols[j].items():v[i]=add(v.get(i,zero),mul(ex(a),ex(x)))
  return v
 def rowbound(row,v):
  q=zero
  for i,x in v.items():q=add(q,mul(row[i],x))
  return q
 envelopes=[]
 for n,c0 in enumerate(combs):
  co={j:F(x) for j,x in c0 if j<640};cn={j-640:F(x) for j,x in c0 if j>=640}
  sr=weighted(Sr,cn);se=weighted(Se,cn);srr=weighted(rr,cn);sre=weighted(re,cn);ar=weighted(oldr,co);ae=weighted(olde,co);ch=[]
  for t in range(2):
   q=[zero,zero]
   for j,a in co.items():
    for k in range(2):q[k]=add(q[k],mul(ex(a),oldg[j][t][k]))
   cr={i:x for i,x in enumerate(CF[t][0]) if x};ce={i:x for i,x in enumerate(CF[t][1]) if x}
   q[0]=add(q[0],add(bdot(cr,ar),add(rowbound(R[t],srr),rowbound(P[t],sr))))
   q[1]=add(q[1],add(add(bdot(cr,ae),bdot(ce,ar)),add(rowbound(R[t],sre),add(rowbound(P[t],se),rowbound(Q[t],sr)))))
   ch.append([[str(D),str(U)] for D,U in q])
  envelopes.append(ch)
 def rec(v):return [[i,str(x)] for i,x in sorted(v.items())]
 payload={'selected99_original768_constant_combinations':combs,'Gamma_two':src['Gamma_two'],'Cfull_two_by2_by56':C,'typed_program50':src['typed_program'],'old80_root_columns':src['needed_original_old_columns'],'physical24_plane5_values':sorted(leaves.items()),'source_order_skip_seven_products':sorted(skip),'honest327554_opaque_typed_leaf':honest,'gB_plane_order3_typed_leaf':gb,'old80_r_first_columns':[[j,rec(v)] for j,v in sorted(oldr.items())],'old80_ell_first_columns':[[j,rec(v)] for j,v in sorted(olde.items())],'native_L0':raw['native_L0'],'native_Lr':raw['native_Lr'],'native_Le':raw['native_Le'],'native23_Sr':[[j,rec(Sr[j])] for j in sorted(nids)],'native23_Se':[[j,rec(Se[j])] for j in sorted(nids)],'native23_S2rr':[[j,rec(rr[j])] for j in sorted(nids)],'native23_S2rell':[[j,rec(re[j])] for j in sorted(nids)],'positive396_Kfull_envelopes':envelopes,'all128_B11_intrinsic_S0_Sr_entries_match_current':True,'all48_Lr_Le_C16_match_current':True,'all99_weighted_nativeS0zero_source_order_peer_pinned':True,'formula':'R=GammaL0;P=GammaLr+C_rL0;Q=GammaLe+C_ellL0. Krr=OLDGamma2rr+C_rOLDfirst_r+R*S2rr+P*Sr;Krell=OLDGamma2rell+C_rOLDfirst_ell+C_ellOLDfirst_r+R*S2rell+P*Se+Q*Sr. Allweightedselected99wholecombinations.','actual_K_or_upper_values_not_evaluated':True}
 numeric=count(payload);future=2*numeric+count([CF,G,L0,Lr,Le,R,P,Q,rr,re,oldr,olde,jets,oldg])+10000;assert future<60000
 rawbound=2_000_000
 ep=out/'pruned99-Kfull-P5-native23-and-prior396-input.json';H.dump(ep,payload)
 report={'id':STEM,'all_passed':True,'source_cache_prior_only_no_actual_K':True,'all128_B11_S0_Sr_matches':True,'all_Lr_Le_C16_matches':True,'old_physical_boundaries24':24,'native_source_columns23':23,'positive_committed_prior_bits':maximum,'temporary_prior_bits':3*maximum+1,'source_builder_owned_bound':sourceowned,'source_maximum_selected_field_numeric':maximum_selected_field_numeric,'future_value_numeric_bound':future,'pruned_numeric_records':numeric,'future_RAW_output_bound':rawbound,'actual_K_secondresponse_or_upper_values_not_computed':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_committed_prior_bits','temporary_prior_bits','source_builder_owned_bound','future_value_numeric_bound','pruned_numeric_records','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
