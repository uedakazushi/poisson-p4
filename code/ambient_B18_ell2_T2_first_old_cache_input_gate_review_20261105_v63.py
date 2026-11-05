#!/usr/bin/env python3
"""Independent source-availability gate for old rell callbacks; no second evaluation."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual104-old-C0-ancestor-cache-bindings-and-Rell4-gate.json',
'formal':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json',
'formal_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63_REVIEW/report.json',
'prev':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
'prev_peer':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json',
'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
'common':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json',
'base':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json',
'common_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json',
'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
'quartic':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
'quartic_peer':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json',
'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'LITERAL_SOURCE_AVAILABILITY_AND_RING_QUOTIENT_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k not in {'FJ','C0','quartic'}}
 for k in ['primary','formal_peer','prev_peer','common_peer','quartic_peer']:assert D[k]['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];ids=D['formal']['unknown_second_label_list'];catalog={i:(v,g) for i,v,g in D['prev']['all108_composite_values_and_opaque_second_gradients']}
 assert sorted(set(ids)-set(catalog))==[0] and all(not catalog[i][1] for i in ids if i) and E['formal_original_C0_second_labels104']==ids
 bmon=[tuple(m) for m in D['common']['basis56']];assert bmon==list(map(tuple,D['base']['basis56']));ix={m:i for i,m in enumerate(bmon)}
 target=[(0,0),(1,0),(0,1),(1,1)];index={m:i for i,m in enumerate(target)}
 images=[index[(m[7],m[8])] if not any(m[:7]) and (m[7],m[8]) in index else None for m in bmon]
 actual=[next(i for i,v in enumerate(images) if v==j) for j in range(4)];assert actual==E['B56_to_Rell4_coefficient_indices']
 for i,a in enumerate(bmon):
  for j,b in enumerate(bmon):
   prod=ix.get(tuple(x+y for x,y in zip(a,b)));left=images[prod] if prod is not None else None
   right=None
   if images[i] is not None and images[j] is not None:
    aa,bb=target[images[i]],target[images[j]];right=index.get((aa[0]+bb[0],aa[1]+bb[1]))
   assert left==right
 fj=fields(ROOT/P['FJ'],['both57_source2','all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']);CB=fj['both57_source2']['CB'];bind=dict(D['copy']['explicit57_source_copy_bindings']);assert set(bind)==set(CB['unknown57']) and D['copy']['all57_copied_equations_literal_match']
 vals={i:sp(v).get(0,Q(0)) for i,v in D['base']['all164_original_source_roots']}
 def put(i,c):
  if i in vals:assert vals[i]==c
  vals[i]=c
 for i,v in zip(D['common']['unknown105_order'],D['common']['all105_cached_values_B56']):put(i,sp(v).get(0,Q(0)))
 for i,p in CB['all57_source2jets']:put(bind[i],Q(p[0]))
 opaque={int(i):Q(p[0]) for f in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA'] for i,p in fj[f].items()}
 for b in D['mu_gate']['all35_ordered_prefix_bindings']:opaque[b['node']]=sp(b['origin_and_ds']).get(0,Q(0))
 regular=fields(ROOT/P['quartic'],['original_regular_function_nodes_unmodified'])['original_regular_function_nodes_unmodified']
 for i in regular['old_CB']+regular['new_CA']:opaque[i]=Q(0)
 C=fields(ROOT/P['C0'],['DAG_scalar_nodes','actual_C0_basis8_scalar_nodes']);nodes=C['DAG_scalar_nodes'];assert nodes[0]==['c',0,1]
 for i in C['actual_C0_basis8_scalar_nodes']:put(i,Q(0))
 stops=set(vals)|set(opaque);roots=ids+CB['F9_nodes'][7:9]+D['copy']['all57_copied_equation_indices'];need=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in stops:continue
  n=nodes[i]
  if n[0] in ('a','m'):todo.extend(n[1:3])
  elif n[0] in ('n','i'):todo.append(n[1])
  else:assert n[0]=='c',(i,n)
 assert [[i,nodes[i]] for i in sorted(need)]==E['source_node_closure'] and len(need)==14541
 values={};deps={};inverses=[]
 for i in sorted(need):
  n=nodes[i]
  if i in vals:v=vals[i];d=set()
  elif i in opaque:v=opaque[i];d={i}
  elif n[0]=='c':v=Q(n[1],n[2]);d=set()
  elif n[0]=='a':v=values[n[1]]+values[n[2]];d=deps[n[1]]|deps[n[2]]
  elif n[0]=='n':v=-values[n[1]];d=deps[n[1]]
  elif n[0]=='i':assert values[n[1]];v=1/values[n[1]];d=deps[n[1]];inverses.append([i,n[1],str(values[n[1]])])
  else:
   assert n[0]=='m';a,b=n[1:];v=values[a]*values[b];d=(deps[a] if values[b] else set())|(deps[b] if values[a] else set())
  values[i]=v;deps[i]=d
 assert inverses==E['all_unit_inverse_origins_nonzero'] and [[i,sorted(deps[i])] for i in sorted(need)]==E['source_second_dependency_labels']
 expected={'selected_original_DAG_closure4coeff':4*len(need),'sparse_opaque_second_support':sum(map(len,deps.values())),'all105_and164_B56_input_nonzeros':sum(len(v) for v in D['common']['all105_cached_values_B56'])+sum(len(v) for i,v in D['base']['all164_original_source_roots']),'ownCB57_source2_input_dense':57*55,'selectedFULL9firstcache_max':len(need)*3,'all619_two_old_knownparts_and_zero_sens_upper':1478*2+2*619,'equations_copy_closure_and_metadata':12000}
 assert expected==E['slot_breakdown'] and sum(expected.values())==123117==E['whole_slot_bound']<150000
 ep=o/'independent104-source-callbacks-and-Rell4-input-gate.json';dump(o,ep.name,{'inputs':inputs,'all103_literal_catalog_matches_and_zero_opaque_second_gradients':True,'remaining_label0_literal_zero':True,'all3136_B56_to_Rell4_product_identities_verified':True,'all14541_original_ID_ancestor_nodes_and_all_unit_inverse_origins_verified':True,'actual_common105_and_zero7_copiedCB57_same_source_binding':True,'all_unknown_honest_second_coefficients_must_remain_symbolic_in_future_worker_until_gradients_zero':True,'whole_slot_bound':123117,'no_second_response_or_source_implicit_gauge_or_native_values_computed':True,'pair_independent_gradient_proof':'For a regular source expression f(u,h), the coefficient of each unknown h_rell in f_rell is partial_h f at the same center. It is unchanged when first directions are changed. Known cross-terms must be recomputed with actual full9 firsts.'})
 report={'id':STEM,'all_passed':True,'source_availability_gate_for_old_rell_cache_composition_accepted':True,'actual_second_or_HF1_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
