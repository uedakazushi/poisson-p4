#!/usr/bin/env python3
"""Source order pruning for all99 OLD K callbacks; point kernels only."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 p='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63/';r=read(p+'report.json');assert r['all_passed'] and not r['unknown_boundaries'] and not r['inverse_nodes']
 d=read(p+'actual99-selected-source-K-callback-roots-and-literal-ancestry.json');ops=dict(d['typed_program'])
 tail='reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63/';t=read(tail+'literal-gB-tail-orders-and-existing-source3-source2-selectors.json');assert read(tail+'report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 z={v['zB_node'] for v in t['ordered_source_rows']};assert z=={27897,27948,28203,28254,28437,28692,28794}
 assert all(v['zB_point_r_ell_exact0'] and v['h_point_exact0'] for v in t['ordered_source_rows'])
 old='reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63/';o=read(old+'case508-original595-source-order-and-plane-cut-proof.json');assert o['gB_plane_exact_node']==299650 and o['both_zero7_conductor_functions_have_zero_all6_source2_coefficients']
 assert read('reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json')['all_passed']
 honest=next(i for i,n in ops.items() if n==['actual_C0_scalar_function',327554]);gb=next(i for i,n in ops.items() if n==['actual_C0_scalar_function',299650]);zt={i for i,n in ops.items() if n[0]=='actual_C0_scalar_function' and n[1] in z}
 parents=[]
 for i,n in ops.items():
  if n[0]=='product' and honest in n[3]:
   assert len(n[3])==2 and len(set(n[3])&zt)==1;parents.append(i)
  elif n[0]=='sum':assert honest not in {j for j,a,b in n[1]}
 assert len(parents)==7
 top_roots={arg for j,col in d['needed_original_old_columns'] for row,(kind,arg) in col if kind=='typed'}
 assert honest not in top_roots and d['literal_C0_program']==[[0,['c',0,1]]]
 order={}
 for i,n in sorted(ops.items()):
  if n[0]=='actual_C0_scalar_function':order[i]=1 if i==honest else 3 if i==gb else 2 if i in zt else 0
  elif n[0]=='product':order[i]=sum(order[j] for j in n[3])
  elif n[0]=='sum':order[i]=min((order[j] for j,a,b in n[1]),default=100)
  else:raise AssertionError(('unexpectedsourcekind',i,n))
 assert all(order[i]>=3 for i in parents) and order[gb]>=3
 # This point-source verification does not evaluate any second response.
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';S0=list(map(H.parse,field(mug,'native_S0_columns')));assert len(S0)==128
 maximum=0
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());assert maximum<=384;return q
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def exact(x):return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 native=[{j-640:F(c) for j,c in v if j>=640} for v in d['selected99_original768_constant_combinations']]
 for c in native:
  bounds={}
  for j,x in c.items():
   for i,y in S0[j].items():bounds[i]=add(bounds.get(i,(1,0)),mul(exact(x),exact(y)))
 H.dump(out/'phase.json',{'phase':'source_orders_and_native_point_positive_bound_before_exact_pointzero','prior_point_bits':maximum,'no_K_second_or_first_response_evaluated':True})
 for c in native:
  v={}
  for j,x in c.items():
   for i,y in S0[j].items():v[i]=v.get(i,F(0))+x*y;assert max(abs(v[i].numerator).bit_length(),v[i].denominator.bit_length())<=max(1,maximum)
  assert not any(v.values())
 literal=[i for i,kind,n in d['literal_cache_boundaries'] if kind=='source164'];assert len(literal)==24
 qp='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';mons=field(qp,'plane9_monomials');pi=field(qp,'raw33_plane9_map_columns');selectors=[]
 for rawi,p9i,mon in [(0,0,[0,0]),(8,1,[1,0]),(9,2,[0,1]),(24,3,[2,0]),(25,4,[1,1])]:
  assert mons[p9i]==mon and pi[rawi]==[[p9i,'1']] and [j for j,col in enumerate(pi) if any(k==p9i and x!='0' for k,x in col)]==[rawi];selectors.append([rawi,p9i,mon])
 assert read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63_REVIEW/report.json')['all_passed']
 ep=out/'all99-old-plane-order-pruning-and-nativeS0-source-gate.json';ev={'honest_function327554_typed_leaf':honest,'only_seven_honest_parent_products':[[i,ops[i],order[i]] for i in parents],'gB299650_typed_leaf':gb,'gB_plane_order_at_least3_source_bound':True,'h327554_order_at_least1_and_zB7_order_at_least2':True,'all99_intrinsic_native_full48_S0_combinations_exact0':True,'native_point_positive_prior_bits':maximum,'all24_remaining_physical_boundary_ids':literal,'raw33_to_plane_degree2_selector_map':selectors,'plane5_basis':[[0,0],[1,0],[0,1],[2,0],[1,1]],'allOLD99_secondresponses_descend_to_physical24_plane5_with_source_order_zero_pruning':True,'no_actual_second_K_b1_or_upperresponse_values':True,'scope':'Sourcefunctionorder+pointkernelsONLY;futureactualKcachebundleandall918rightsolveSTILLseparate.'};H.dump(ep,ev)
 owned=count([d,ops,t,o,S0,native,order,literal,mons,pi,ev,pins])+15000;assert owned<60000
 report={'id':STEM,'all_passed':True,'source_orders_and_native_origin_only':True,'all99_nativeS0_exact0':True,'honest327554_plane2_independence_by_order':True,'gB_plane2_exact0_by_order':True,'actual_K_or_secondresponses_computed':False,'source_owned_bound':owned,'point_positive_prior_bits':maximum,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','all99_nativeS0_exact0','source_owned_bound','point_positive_prior_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
