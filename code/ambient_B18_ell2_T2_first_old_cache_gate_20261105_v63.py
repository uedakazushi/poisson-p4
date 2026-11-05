#!/usr/bin/env python3
"""Literal ancestor/source gate for all104 formal old rell callbacks."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 formal=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json');roots=formal['unknown_second_label_list']
 prev=read('reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json');old={i:(v,g) for i,v,g in prev['all108_composite_values_and_opaque_second_gradients']}
 assert set(roots)-set(old)=={0} and all(not old[i][1] for i in roots if i)
 copy=read('reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json')
 read('reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json')
 common=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json')
 base=read(H.PATHS['B56']);bmon=[tuple(m) for m in common['basis56']];assert common['basis56']==base['basis56']
 # B56 -> k[r,ell]/(r²,ell²,z7) is a literal monomial quotient.
 wanted=[]
 for r,e in [(0,0),(1,0),(0,1),(1,1)]:m=[0]*9;m[7]=r;m[8]=e;wanted.append(bmon.index(tuple(m)))
 vals={i:H.parse(v).get(0,F(0)) for i,v in base['all164_original_source_roots']}
 vals.update({i:H.parse(v).get(0,F(0)) for i,v in zip(common['unknown105_order'],common['all105_cached_values_B56'])})
 hs[H.PATHS['firstjets']]=H.sha(ROOT/H.PATHS['firstjets']);J=H.selected(H.PATHS['firstjets'],'both57_source2')['CB'];bindings=dict(copy['explicit57_source_copy_bindings'])
 assert set(bindings)==set(J['unknown57'])
 for i,p in J['all57_source2jets']:vals[bindings[i]]=F(p[0])
 opaque={}
 for field in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  opaque.update({int(i):F(p[0]) for i,p in H.selected(H.PATHS['firstjets'],field).items()})
 g=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json')
 for b in g['all35_ordered_prefix_bindings']:opaque[b['node']]=H.parse(b['origin_and_ds']).get(0,F(0))
 hs[H.PATHS['quartic']]=H.sha(ROOT/H.PATHS['quartic']);gids=H.selected(H.PATHS['quartic'],'original_regular_function_nodes_unmodified')
 for i in gids['old_CB']+gids['new_CA']:opaque[i]=F(0)
 hs[H.PATHS['C0']]=H.sha(ROOT/H.PATHS['C0']);basis=H.selected(H.PATHS['C0'],'actual_C0_basis8_scalar_nodes')
 for i in basis:vals[i]=F(0)
 nodes=H.selected(H.PATHS['C0'],'DAG_scalar_nodes');assert nodes[0]==['c',0,1]
 stops=set(vals)|set(opaque);need=set();todo=roots+J['F9_nodes'][7:9]+copy['all57_copied_equation_indices']
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in stops:continue
  n=nodes[i]
  if n[0] in ['a','m']:todo.extend(n[1:3])
  elif n[0] in ['n','i']:todo.append(n[1])
  else:assert n[0]=='c',('exact_missing_noncached_source_function',i,n)
 value={};deps={};inverse_origins=[]
 for i in sorted(need):
  n=nodes[i]
  if i in vals:v=vals[i];dep=set()
  elif i in opaque:v=opaque[i];dep={i}
  elif n[0]=='c':v=F(n[1],n[2]);dep=set()
  elif n[0]=='a':v=value[n[1]]+value[n[2]];dep=deps[n[1]]|deps[n[2]]
  elif n[0]=='n':v=-value[n[1]];dep=deps[n[1]]
  elif n[0]=='i':
   assert value[n[1]];v=1/value[n[1]];dep=deps[n[1]];inverse_origins.append([i,n[1],str(value[n[1]])])
  else:
   assert n[0]=='m';v=value[n[1]]*value[n[2]];dep=(deps[n[1]] if value[n[2]] else set())|(deps[n[2]] if value[n[1]] else set())
  H.guard({0:v} if v else {});value[i]=v;deps[i]=dep
 slots={'selected_original_DAG_closure4coeff':4*len(need),'sparse_opaque_second_support':sum(map(len,deps.values())),
  'all105_and164_B56_input_nonzeros':sum(len(v) for v in common['all105_cached_values_B56'])+sum(len(v) for i,v in base['all164_original_source_roots']),
  'ownCB57_source2_input_dense':57*55,'selectedFULL9firstcache_max':len(need)*3,
  'all619_two_old_knownparts_and_zero_sens_upper':1478*2+2*619,
  'equations_copy_closure_and_metadata':12000}
 total=sum(slots.values());assert total<150000
 data={'formal_original_C0_second_labels104':roots,'all103_in_accepted108_catalog_gradient_zero':True,'literal_zero_node0':nodes[0],
  'B56_to_Rell4_coefficient_indices':wanted,'Rell4_basis':['1','rA','ell','rA*ell'],
  'source_node_closure':[[i,nodes[i]] for i in sorted(need)],'source_second_dependency_labels':[[i,sorted(deps[i])] for i in sorted(need)],
  'all_unit_inverse_origins_nonzero':inverse_origins,'all_unit_or_honest_leaves_bound_by_literal_existing_ID':True,
  'copied57_exact_source_bindings':copy['explicit57_source_copy_bindings'],'copied57_equations':copy['all57_copied_equation_indices'],
  'copy_rule':'ownCB57source2 composedwithactualYB7/8fromcommon105, zero7retainedinverseinstance; nonlinearadaptationnotfrozenlinear',
  'native_contributions_remain_external':True,
  'source_theorem':'opaque secondgradient depends onlycenters, independentofr/ellfirstdirections. Reconstructall104rootswithsameoriginalDAG andliteralstops, retainallopaqueh2untilgradzeroisverifiedagain. Newrellknownpartsusefull9firsts, notoldB11numericalvalues.',
  'whole_slot_bound':total,'slot_breakdown':slots,'actual_second_cache_composition_not_started':True,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384}}
 ep=out/'actual104-old-C0-ancestor-cache-bindings-and-Rell4-gate.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'actual_second_cache_composition_started':False,'source_closure_nodes':len(need),
  'whole_slot_bound':total,'all104_literal_leaves_source_bound':True,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_closure_nodes','whole_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
