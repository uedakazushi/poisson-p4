#!/usr/bin/env python3
"""Source-order theorem for the one old case outside first-only548."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 for stem in ['AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW','AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW','AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW','AMBIENT_B18_ELL2_T2_MIXED619_OLD_ALIAS_ORIGIN_20261105_v63']:
  assert read('reports/'+stem+'/report.json')['all_passed']
 c0=read('reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',['actual_C0_basis8_scalar_nodes','all35_regular_coordinate_divisions'])
 g=read('reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',['CB_g_restriction_node','CB_f63_node'])
 chain=read('reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',['accepted_actual_regular_function_identities'])
 tail=next(x for x in chain['accepted_actual_regular_function_identities'] if x['kind']=='accepted_actual_CB_f63_tail')
 assert tail['gB']==g['CB_g_restriction_node']==299650 and tail['function']==g['CB_f63_node']==28539
 f=read('reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',['source2_monomials','both57_source2'])
 cut=[j for j,mon in enumerate(f['source2_monomials']) if set(mon)<=set([7,8])];assert len(cut)==6
 first=c0['all35_regular_coordinate_divisions'][0];assert first['restriction_r_s_node']==c0['actual_C0_basis8_scalar_nodes'][7]==505938
 assert first['original_pure_CA_source_function']==f['both57_source2']['CA']['f63']
 assert tail['function']==f['both57_source2']['CB']['f63']
 for family in ['CA','CB']:
  data=f['both57_source2'][family];jet=list(map(Q,data['pure_function2jets'][str(data['f63'])]));assert all(jet[j]==0 for j in cut)
 # Accepted pointed source identities identify these exact zero7 functions;
 # centered full9 pullback preserves their order>=3.
 coord=set(c0['actual_C0_basis8_scalar_nodes'][:7]);high={299650,505938}
 m=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes']);B=list(map(dict,m['old640_source_columns']));ops=dict(m['old640_typed_closure_nodes'])
 gg=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json');Gamma=list(map(sp,gg['Gamma4_on_source56_times8']));pairs=gg['ordinary_second_pair_order8']
 cc=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'][508]
 assert cc['original768_combination']==[[595,'1']] and not cc['actual_a1_27_times9']
 memo={}
 def order(k,plane):
  key=(k,plane)
  if key in memo:return memo[key]
  n=ops[k];op=n[0]
  if op=='c':v=0 if n[1] else 999
  elif op=='actual_C0_scalar_function':v=999 if plane and n[1] in coord else 3 if n[1] in high else 1
  elif op=='product':v=sum(order(j,plane) for j in n[3]) if n[1] else 999
  else:assert op=='sum';v=min([order(j,plane) for j,a,b in n[1] if a]+[999])
  memo[key]=min(v,999);return memo[key]
 receipts=[]
 for target,row in enumerate(Gamma):
  for slot,c in row.items():
   i,pair=divmod(slot,8);root=B[595].get(i)
   if root is None:continue
   plane=all(v in [7,8] for v in pairs[pair]);bound=order(root,plane)
   assert bound>=3,(target,i,pairs[pair],root,bound)
   receipts.append([target,i,pair,root,'zero7-plane' if plane else 'ordinary-total',bound])
 ordinary29=order(B[595][29],False);plane29=order(B[595][29],True)
 assert ordinary29==2 and plane29>=3 and all(not any(slot//8==29 for slot in Gamma[t]) for t in [2,3])
 nodes=sorted({k for k,p in memo});owned=count([c0,g,chain,f,m,gg,cc,memo,receipts])+6000;assert owned<150000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'case508-original595-source-order-and-plane-cut-proof.json'
 ev={'case619':508,'original_old640_source':595,'no_a1_pivot_carry':True,'gA_plane_exact_node':505938,'gB_plane_exact_node':299650,'old_f63_distinct_node':28539,
 'both_zero7_conductor_functions_have_zero_all6_source2_coefficients':True,'source_coordinate_nodes7':sorted(coord),
 'all_Gamma_slot_order_receipts':receipts,'selected_typed_nodes':[[k,ops[k]] for k in nodes],
 'special_row29_root281_total_order_lower':ordinary29,'special_row29_zero7_plane_order_lower':plane29,
 'all_four_lower_coefficients_zero_by_source_order':True,'no_numerical_quadratic_image_or_newjet_job':True,
 'proof':'EveryGamma-supportedsecondcoefficientvanishesbyorder>=3. Theonlytotal-order2root281isusedonlyatrr/rell;itszero7restrictionhasorder>=3becauseallothertermsarezi-timesregularfunctions. gA/gBownzero7order>=3followsall6savedsource2cutcoefficientsandactualpointedrestrictionbindings. Centeredpullbackpreservesorder. a1=0sofirstcarryabsent.',
 'scope':'Onlyfourlowerchannelsforcase508. Withseparatelypeered548lowerzeros,thisclosesall549pureOLDlowerchannels;upperchannelsand70native/mixedcyclesstayopen.'}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all_four_lower_coefficients_zero_by_order':True,'no_new_numerical_image_or_sourcejet':True,'source_owned_bound':owned,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'receipt_count':len(receipts),'nodes':len(nodes),'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
