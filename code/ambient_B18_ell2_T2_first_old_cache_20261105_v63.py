#!/usr/bin/env python3
"""One existing-cache Rell4 ancestor composition, no source solve."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63'
def plus(a,b):return [x+y for x,y in zip(a,b)]
def scal(a,c):return [c*x for x in a]
def mul(a,b):return [a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0],a[0]*b[3]+a[3]*b[0]+a[1]*b[2]+a[2]*b[1]]
def inv(a):
 assert a[0];c=a[0];v=[1/c,-a[1]/c**2,-a[2]/c**2,2*a[1]*a[2]/c**3-a[3]/c**2];assert mul(a,v)==[F(1),F(0),F(0),F(0)];return v
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorizedRell4cachecomposition_noimplicit'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 gr=read(gate);assert gr['all_passed'] and not gr['actual_second_cache_composition_started'];ip=str(Path(gate).parent/'actual104-old-C0-ancestor-cache-bindings-and-Rell4-gate.json');D=read(ip);assert hs[ip]==gr['evidence_sha256'][Path(ip).name]
 peer=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_INPUT_GATE_20261105_v63_REVIEW/report.json');assert peer['all_passed']
 for p,h in gr['input_sha256'].items():assert H.sha(ROOT/p)==h;hs[p]=h
 B=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json');B2=read(H.PATHS['B56']);mons=[tuple(m) for m in B['basis56']];index={m:i for i,m in enumerate(mons)}
 ringmons=[(0,0),(1,0),(0,1),(1,1)];ri={m:i for i,m in enumerate(ringmons)};images=[]
 for m in mons:
  v=[F(0)]*4
  if not any(m[:7]) and (m[7],m[8]) in ri:v[ri[m[7],m[8]]]=F(1)
  images.append(v)
 for i,m in enumerate(mons):
  for j,n in enumerate(mons):
   k=index.get(tuple(a+b for a,b in zip(m,n)));assert mul(images[i],images[j])==(images[k] if k is not None else [F(0)]*4)
 def project(p):
  p=H.parse(p);return [p.get(i,F(0)) for i in D['B56_to_Rell4_coefficient_indices']]
 known={i:project(p) for i,p in B2['all164_original_source_roots']};known.update({i:project(p) for i,p in zip(B['unknown105_order'],B['all105_cached_values_B56'])})
 assert B['all9_shared_physical_inverse_source2_exact'] and B['all105_commonplane_lower_overlaps_exact']
 nodes=dict(D['source_node_closure']);YB=H.selected(H.PATHS['firstjets'],'both57_source2')['CB'];hs[H.PATHS['firstjets']]=H.sha(ROOT/H.PATHS['firstjets']);words=H.selected(H.PATHS['firstjets'],'source2_monomials')
 def evaluate(roots,bind,opaque=None):
  opaque=opaque or {};need=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i in need:continue
   need.add(i)
   if i in bind or i in opaque:continue
   n=nodes[i]
   if n[0] in ['a','m']:todo.extend(n[1:3])
   elif n[0] in ['n','i']:todo.append(n[1])
   else:assert n[0]=='c'
  values={};grad={}
  for i in sorted(need):
   if i in bind:v=bind[i];g={}
   elif i in opaque:v=opaque[i];g={i:F(1)}
   else:
    n=nodes[i];op=n[0]
    if op=='c':v=[F(n[1],n[2]),F(0),F(0),F(0)];g={}
    elif op=='a':v=plus(values[n[1]],values[n[2]]);g=add(grad[n[1]],grad[n[2]])
    elif op=='n':v=scal(values[n[1]],-1);g=scale(grad[n[1]],-1)
    elif op=='m':v=mul(values[n[1]],values[n[2]]);g=add(scale(grad[n[1]],values[n[2]][0]),grad[n[2]],values[n[1]][0])
    else:assert op=='i';v=inv(values[n[1]]);g=scale(grad[n[1]],-1/values[n[1]][0]**2)
   H.guard({k:c for k,c in enumerate(v) if c});H.guard(g);values[i]=v;grad[i]=g
  assert sum(4+len(grad[i]) for i in need)+20000<150000
  return {i:values[i] for i in roots},{i:grad[i] for i in roots},len(need)
 y,_,yn=evaluate(YB['F9_nodes'][7:9],known);args=[[F(0)]*4 for _ in range(7)]+[y[i] for i in YB['F9_nodes'][7:9]];bindcopy=dict(D['copied57_exact_source_bindings'])
 for i,p in YB['all57_source2jets']:
  v=[F(0)]*4
  for word,c in zip(words,p):
   c=F(c)
   if not c:continue
   term=[c,F(0),F(0),F(0)]
   for j in word:term=mul(term,args[j])
   v=plus(v,term)
  known[bindcopy[i]]=v
 ceq,_,cn=evaluate(D['copied57_equations'],known);assert all(not any(v) for v in ceq.values())
 # FULL9 honest firsts, never old B11 numeric values in the new direction.
 opaque={}
 for field in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in H.selected(H.PATHS['firstjets'],field).items():opaque[int(i)]=[F(p[0]),F(p[8])+F(p[9]),F(p[9])/5,F(0)]
 atom=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json');mu=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json');prefix={p['node']:H.parse(p['origin_and_ds']) for p in mu['all35_ordered_prefix_bindings']}
 for p in atom['all35_ordered_CA_prefix_firsts']:
  v=H.parse(p['first9_bar']);opaque[p['node']]=[v.get(0,F(0)),v.get(8,F(0)),prefix[p['node']].get(1,F(0))/5,F(0)]
 gids=H.selected(H.PATHS['quartic'],'original_regular_function_nodes_unmodified');gradients=H.selected(H.PATHS['quartic'],'all8_regular_function_first_gradients_in_FULL9_CA_F9')
 for i,g in zip(gids['old_CB']+gids['new_CA'],gradients):opaque[i]=[F(0),F(g[7])+F(g[8]),F(g[8])/5,F(0)]
 for i in H.selected(H.PATHS['C0'],'actual_C0_basis8_scalar_nodes'):known[i]=[F(0)]*4
 values,grads,n=evaluate(D['formal_original_C0_second_labels104'],known,opaque);assert all(not g for g in grads.values())
 formal=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json');actual=[]
 for row in formal['all619_results']:
  v=[]
  for part in row['two_formal_old_functionals']:
   a=F(part['known_lower_first_response_only']);g={}
   for i,c in part['literal_original_C0_second_rell_sensitivity']:a+=F(c)*values[i][3];g=add(g,grads[i],F(c))
   assert not g;H.guard({0:a} if a else {});v.append(a)
  actual.append(v)
 ep=out/'actual104-source-rell-values-and-two-by619-old-knownparts.json';H.dump(ep,{'Rell4_basis':['1','rA','ell','rA*ell'],
  'all3136_B56_ring_projection_products_verified':True,'all105_equations_zero_by_literal_ring_homomorphism_from_accepted_B56':True,
  'all57_copy_equations_zero_directly_verified':True,'actual_YB7_8_from_common105':[[i,list(map(str,y[i]))] for i in YB['F9_nodes'][7:9]],
  'all104_original_C0_values_and_opaque_gradients':[[i,list(map(str,values[i])),srec(grads[i])] for i in sorted(values)],
  'all619_two_old_knownparts_including_frozen_pivot_first_correction':[list(map(str,v)) for v in actual],
  'all_honest_second_rell_sensitivities_exactly_zero':True,'native_second_and_L1ell_Omega_r_terms_still_external':True,
  'actual_complete_two_by619_HF_rank_not_claimed':True,'no_new_physical_implicit_or_gauge_coefficients':True})
 count=4*n+sum(map(len,grads.values()))+sum(sum(bool(x) for x in v) for v in actual)+sum(sum(bool(c) for c in v) for v in known.values())+20000;assert count<150000
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'all_honest_second_rell_sensitivities_zero':True,'all104_composite_values_source_bound':True,
  'actual_complete_T2_HF1_or_length_computed':False,'native_values_incomplete':True,
  'new_physical_implicit_or_gauge_coefficients':False,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},
  'owned_QQ_slot_bound':count,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','all_honest_second_rell_sensitivities_zero','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)
