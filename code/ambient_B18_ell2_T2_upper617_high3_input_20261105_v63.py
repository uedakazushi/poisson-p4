#!/usr/bin/env python3
"""Positive-only directOLD highestthird617 gate; no signed gradients/values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};maximum=0;cursor={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 sup=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_SUPPORT_PREFLIGHT_20261105_v63/factored617-lower-kernel-direct-third-support-and-source-frontier.json');p8=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_P8_ANCESTOR_SOURCE_20261105_v63/literal86-P8-composite-program-and-exact-boundary-catalog.json');assert not p8['uncached_or_unrecognized_opaque_boundaries'] and not p8['source_inverse_nodes_remain_symbolic']
 typed=dict(sup['pruned_typed_source_closure']);c0=dict(p8['selected_original_C0_arithmetic_nodes']);tids={n[1] for n in typed.values() if n[0]=='actual_C0_scalar_function'};leaves=(tids|{i for i,*x in p8['exact_literal_cache_boundary_nodes']})-set(c0)
 fc=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640');a79={i:H.parse(v).get(0,F(0)) for i,v in fc['all79_atom_values_full10']};del fc
 vals=field('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values');a164={i:H.parse(v).get(0,F(0)) for i,v in vals};del vals
 origins={};classes=[]
 for i in sorted(leaves):
  assert i in a79 or i in a164
  if i in a79 and i in a164:assert a79[i]==a164[i]
  origins[i]=a79[i] if i in a79 else a164[i];classes.append([i,'source164' if i in a164 else 'C0_79_only',str(origins[i])])
 zero=(1,0)
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length())
  if maximum>384:
   H.dump(out/'positive-guard-frontier.json',{'cursor':cursor,'positive_uncommitted_bound':list(map(str,q)),'bits':maximum,'no_signed_or_source_response_values':True});raise RuntimeError('sourcepositive384guard')
  return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else zero
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return mark((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else zero
 def add(a,b,c=(1,1)):
  z=dict(a)
  for k,x in b.items():
   q=plus(z.get(k,zero),mul(c,x))
   if q[1]:z[k]=q
  return z
 def scale(a,c):return {i:mul(x,c) for i,x in a.items() if x[1] and c[1]}
 co={};cg={}
 def C(i):return (co[i],cg[i]) if i in co else (ex(origins[i]),{i:ex(1)})
 for i,n in sorted(c0.items()):
  cursor={'stage':'positive_C0_origin_dual','node':i};op=n[0]
  if op=='c':o=ex(F(n[1],n[2]));g={}
  elif op=='n':o,g=C(n[1]);g=dict(g)
  elif op=='a':a,u=C(n[1]);b,v=C(n[2]);o=plus(a,b);g=add(u,v)
  else:
   assert op=='m';a,u=C(n[1]);b,v=C(n[2]);o=mul(a,b);g=add(scale(u,b),v,a)
  co[i]=o;cg[i]=g
 to={};tg={}
 for i,n in sorted(typed.items()):
  cursor={'stage':'positive_typed_origin_dual','node':i};op=n[0]
  if op=='actual_C0_scalar_function':o,g=C(n[1]);g=dict(g)
  elif op=='c':o=ex(F(n[1],n[2]));g={}
  elif op=='sum':
   o=zero;g={}
   for j,a,b in n[1]:s=ex(F(a,b));o=plus(o,mul(s,to[j]));g=add(g,tg[j],s)
  else:
   assert op=='product';o=ex(F(n[1],n[2]));g={}
   for j in n[3]:g=add(scale(g,to[j]),tg[j],o);o=mul(o,to[j])
  to[i]=o;tg[i]=g
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';B=[dict(v) for v in field(mug,'old640_source_columns')];P=field(H.PATHS['C0'],'honest_D0_relation_columns');gamma=[H.parse(v) for v in sup['Gamma_two']];needed=sorted({j for row in sup['all617_source_and_direct_third_support_recipes'] for j in row['original_old640_support']});pruned={};column_bounds={}
 for j in needed:
  cursor={'stage':'positive_original_old_column_gradient','source':j};roots={};gg=[]
  for ch,row in enumerate(gamma):
   z={}
   for i,x in row.items():
    if i<48:
     k=B[j].get(i)
     if k is None:continue
     roots[i]=['typed',k];v=tg.get(k,{})
    else:k=P[j][i-48+1];roots[i]=['C0',k];v=C(k)[1]
    z=add(z,v,ex(x))
   gg.append(z)
  pruned[j]=sorted(roots.items());column_bounds[j]=gg
 del B,P
 sourceC=sup['original619_768_combinations']
 def case(case,ch):
  z={}
  for j,x in sourceC[case]:
   if j<640:z=add(z,column_bounds[j][ch],ex(x))
  return z
 basis={j:[case(j,ch) for ch in range(2)] for j in [563,609]};globalbounds=[zero,zero];supportentries=0
 def union(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return mark((d,max(a[1]*(d//a[0]),b[1]*(d//b[0]))))
 for rec in sup['all617_source_and_direct_third_support_recipes']:
  cursor={'stage':'positive617_residual_gradient','case':rec['free_case619']}
  for ch in range(2):
   z={}
   for f in rec['constant_case619_source_factors']:
    caseid=f['case619'];v=basis[caseid][ch] if caseid in basis else case(caseid,ch);z=add(z,v,ex(f['scale_operator'][1]))
   supportentries+=len(z)
   for q in z.values():globalbounds[ch]=union(globalbounds[ch],q)
 assert 3*maximum+1<=1153
 data={'case_source_recipes617':[{'free_case619':v['free_case619'],'constant_case619_source_factors':v['constant_case619_source_factors']} for v in sup['all617_source_and_direct_third_support_recipes']],
 'original619_source768_combinations':sourceC,'pruned_original_old_columns_roots':sorted(pruned.items()),'Gamma_two':sup['Gamma_two'],'typed_source_program':sorted(typed.items()),'literal_C0_program':sorted(c0.items()),'literal_boundary_origins_and_membership':classes,
 'positive_common_two_row_envelopes':[[str(x) for x in q] for q in globalbounds],'positive_prior_committed_bits':maximum,'temporary_prior_bits':3*maximum+1,'positive_sparse_node_gradient_entries':sum(map(len,cg.values()))+sum(map(len,tg.values())),'positive_sparse_old_column_gradient_entries':sum(sum(map(len,g)) for g in column_bounds.values()),'positive617_two_row_gradient_entries':supportentries,
 'formal_scope_ONLY':'directOLDhighestthirdindependentsymbolweightsAFTERactual563609sourcesubtraction. Native,b1,A2,M1,knowncrosspartsNOTevaluated.79-onlyboundarythirdsremainunknown;noRAW33availabilityguess.',
 'no_signed_weighted_gradient_or_scalar_response_values':True,'no_higher_source_order_necessity_claim':True}
 node_entries=data['positive_sparse_node_gradient_entries'];old_entries=data['positive_sparse_old_column_gradient_entries'];inp=count(data);future=2*(inp+1000)+4*node_entries+4*old_entries+4*supportentries+12000;raw=2_000_000+future*(2*math.ceil(maximum*math.log10(2))+40)
 sourceowned=count([sup,p8,typed,c0,origins,co,cg,to,tg,column_bounds,data])+25000
 data.update({'future_formal_numeric_bound':future,'future_formal_RAW_bound':raw,'metadata_reserve':1000,'future_includes_input_doublecopy_node_and_column_duals_outputandserialization':True})
 ep=out/'617-direct-old-highestthird-coherent-source-positive-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'positive_source_only':True,'actual_scalar_jet_or_upper_values_not_queried':True,'C0_arithmetic_nodes':len(c0),'typed_nodes':len(typed),'literal_independent_third_boundaries':len(leaves),'boundary79_only_count':sum(v[1]=='C0_79_only' for v in classes),'positive_prior_bits':maximum,'temporary_bits':3*maximum+1,'future_formal_numeric_bound':future,'future_formal_RAW_bound':raw,'source_owned_bound':sourceowned,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start,'future_formal_worker_NOTauthorized_by_this_gate':True};assert sourceowned<150000 and future<150000 and raw<32*1024*1024;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','literal_independent_third_boundaries','boundary79_only_count','positive_prior_bits','future_formal_numeric_bound','future_formal_RAW_bound','source_owned_bound']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
