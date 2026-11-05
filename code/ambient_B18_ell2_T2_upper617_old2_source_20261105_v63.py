#!/usr/bin/env python3
"""Combined OLD second source inventory only; no gradient or response values."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 dp='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63/actual-D-coefficients-two-by10-by56.json';D=field(dp,'D_two_by10_by56');pairs=field(dp,'active_pairs10');assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63_REVIEW/report.json')['all_passed']
 rows={j for ch in D for v in ch for j,x in enumerate(v) if x!='0'};assert rows=={0,1,2,3,4,6,17,18,19,21,22,23,28,29,30,32,37,38} and all(j<48 for j in rows)
 src='reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63/617-direct-old-highestthird-coherent-source-positive-input.json';c619=field(src,'original619_source768_combinations');assert len(c619)==619
 kern='reports/AMBIENT_B18_ELL2_T2_LOWER_RANK2_HF_TAIL_SOURCE_20261105_v63/complete1012-kernel617-and-filtered-r-tail-source-proof.json';recipes=field(kern,'all617_kernel_lower_constant_source_recipes');assert len(recipes)==617
 oldids=sorted({j for c in c619 for j,x in c if j<640});mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';B=[dict(v) for v in field(mug,'old640_source_columns')];allops=dict(field(mug,'old640_typed_closure_nodes'));rootmap=[];roots=set()
 for j in oldids:
  v=[[row,B[j][row]] for row in sorted(rows) if row in B[j]];rootmap.append([j,v]);roots.update(i for row,i in v)
 del B
 outer='reports/AMBIENT_B18_ELL2_T2_UPPER617_PIVOT4_SECOND_SOURCE_20261105_v63/literal-four-OLDpivot-Gamma-plane2-callbacks-and-frontier.json';piv=read(outer);roots.update(i for j,v in piv['four_source_root_columns'] for row,(kind,i) in v if kind=='typed')
 keep=H.closure(allops,roots);ops={i:allops[i] for i in keep};del allops
 classes={i:kind for i,kind,n in piv['cached_boundaries']};known79=set(field(mug,'old79_atom_ids'));raw='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';known164={i for i,v in field(raw,'all_requested164_root_values')};leafids={n[1] for n in ops.values() if n[0]=='actual_C0_scalar_function'}
 unbound=sorted(leafids-known79-known164);aliases=dict(field(mug,'old_exact84_original_arithmetic_nodes'));alias_seen={};todo=list(unbound)
 while todo:
  i=todo.pop()
  if i in known79 or i in known164 or i in alias_seen:continue
  if i not in aliases:continue
  n=aliases[i];alias_seen[i]=n
  if n[0] in ['a','m']:todo.extend(n[1:3])
  elif n[0] in ['n','i']:todo.extend(n[1:2])
  else:assert n[0]=='c'
 remaining=sorted({i for i in unbound if i not in alias_seen}|{j for n in alias_seen.values() for j in (n[1:3] if n[0] in ['a','m'] else n[1:2] if n[0] in ['n','i'] else []) if j not in alias_seen and j not in known79 and j not in known164})
 a1path='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json';a1=field(a1path,'all_zero_bottom_and_Y0_source_combinations');rawa1=[]
 for j,r in enumerate(a1):
  v=[[i,x] for i,x in r['actual_a1_27_times9'] if i//9<4 and i%9 in [7,8] and x!='0']
  if v:rawa1.append([j,v])
 del a1
 assert sum(len(v) for j,v in rawa1)==28
 p8new=[i for i,n in piv['unrecognized_opaque_frontier']];assert p8new==[566400,566401,566402,566403]
 leafclass=[[i,'source164' if i in known164 else 'C0_79_only' if i in known79 else 'exact_alias84' if i in alias_seen else 'UNBOUND'] for i in sorted(leafids)]
 ep=out/'combined-OLD2-D-and-four-pivot-source-factor-inventory.json';ev={'D_two_by10_by56':D,'pair_order10':pairs,'original619_c0_source_combinations':c619,'kernel617_source_recipes':recipes,'D_old_root_columns':rootmap,'typed_program':sorted(ops.items()),'typed_actual_C0_function_classification':leafclass,'exact84_alias_restoration_subset':sorted(alias_seen.items()),'unresolved_function_boundary_ids':remaining,'outer_four_pivot_source':piv,'raw619_old4_plane_a1_nonempty28coefficients':rawa1,'newg4_second_symbols_not_zero_filled':p8new,'source_composition_factor_program':{'D_Xc2':'forrecipej,a,b:sumcases(j,1),(563,-a),(609,-b),thenoriginalold640constantcomboandD_t,pair,row;foldIDENTICALtyped/originalfunctionrootsBEFOREAD.','outer_GammaXp2a1':'SAMErecipeonactualrawa1;subtractGammaXp_rr*a1_ell+GammaXp_rell*a1_r,OLD4onlyhere.NATIVE23L2S0carryseparate.','unknown_second_symbols':'originalactualregularfunctionID+ordinarypair;allg/hindependentuntilwhole617fold;newg4fromC0goperatorareindependent,notformaldivision0/0.','knownfirst_only_terms':'retain-DXp1a1andallnonlinearfirstcross,NOTcomputedinthisformalstage.','direct_old_cubic':'acceptedglobalJ1killedgenuinehonesth2/h3DIRECTOLDONLY;remainingphysical/Gconductorshandledseparately.'},'no_signed_rootweights_gradients_secondvalues_or_upperresponses':True}
 H.dump(ep,ev);owned=count([D,c619,recipes,rootmap,ops,piv,leafids,leafclass,alias_seen,remaining,rawa1,ev,pins])+10000;assert owned<150000
 report={'id':STEM,'all_passed':True,'metadata_source_only_no_formal_gradient':True,'old_source_columns':len(oldids),'rootmap_scalar_roots':sum(len(v) for j,v in rootmap),'typed_nodes':len(ops),'typed_function_boundaries':len(leafids),'only79_function_boundaries':sum(v[1]=='C0_79_only' for v in leafclass),'alias84_nodes':len(alias_seen),'unresolved_boundaries':len(remaining),'outer_newg_opaque4_retained':True,'raw_a1_plane_nonzero_coefficients':28,'source_owned_bound':owned,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','old_source_columns','rootmap_scalar_roots','typed_nodes','typed_function_boundaries','only79_function_boundaries','alias84_nodes','unresolved_boundaries','source_owned_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
