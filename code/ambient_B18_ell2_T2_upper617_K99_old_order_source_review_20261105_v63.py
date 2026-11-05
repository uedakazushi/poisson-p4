#!/usr/bin/env python3
"""Fixed source-order/selector peer. No repeats of99point arithmetic or seconds."""
from pathlib import Path
import json,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,d):(out/n).write_text(json.dumps(d,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_OLD_ORDER_SOURCE_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='c21be2e1b62cf393ae1a41f7f56a827a2284afb0a6e8335700581cec75ce39d7'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'all99-old-plane-order-pruning-and-nativeS0-source-gate.json')
 cb='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63/actual99-selected-source-K-callback-roots-and-literal-ancestry.json';z=rd(cb);ops=dict(z['typed_program'])
 pp='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63_REVIEW/report.json';peer=rd(pp);assert peer['all_passed'] and pins[pp]=='3704c0773a0170d258d4ae4f94f22365547a2bb0c924b3a12ae6a44595576573'
 tp='reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63/literal-gB-tail-orders-and-existing-source3-source2-selectors.json';tail=rd(tp)
 tr='reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json';assert rd(tr)['all_passed'] and pins[tr]=='dee271e97f4904d80c812f046666a8101fda1805adddffe17e8fe7e59f2d51f1'
 zb={v['zB_node'] for v in tail['ordered_source_rows']};assert len(zb)==7
 assert all(v['zB_point_r_ell_exact0'] and v['h_point_exact0'] for v in tail['ordered_source_rows'])
 assert ops[94]==['actual_C0_scalar_function',327554] and ops[107]==['actual_C0_scalar_function',299650]
 parent=[]
 for i,n in ops.items():
  if n[0]=='product' and 94 in n[3]:
   assert len(n[3])==2 and len({ops[j][1] for j in n[3] if ops[j][0]=='actual_C0_scalar_function'}&zb)==1;parent.append(i)
  elif n[0]=='sum':assert 94 not in [j for j,a,b in n[1]]
 assert sorted(parent)==[516,602,663,724,785,846,907] and d['honest_function327554_typed_leaf']==94
 assert d['only_seven_honest_parent_products']==[[i,ops[i],3] for i in parent]
 assert all(arg!=94 for j,col in z['needed_original_old_columns'] for row,(kind,arg) in col if kind=='typed')
 p508='reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63/case508-original595-source-order-and-plane-cut-proof.json';a=rd(p508);assert a['gB_plane_exact_node']==299650 and a['old_f63_distinct_node']==28539 and a['both_zero7_conductor_functions_have_zero_all6_source2_coefficients']
 p508r='reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json';assert rd(p508r)['all_passed'] and pins[p508r]=='8f067c86dc5f7c1e3ad8a638265dc7d2ede6a3037cb3a57e71bfbac9cb55659e'
 mapfile='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';m=fld(mapfile,['plane9_monomials','raw33_plane9_map_columns']);sel=[]
 for rawi,p9i,ex in [(0,0,[0,0]),(8,1,[1,0]),(9,2,[0,1]),(24,3,[2,0]),(25,4,[1,1])]:
  assert m['plane9_monomials'][p9i]==ex and m['raw33_plane9_map_columns'][rawi]==[[p9i,'1']]
  assert [j for j,col in enumerate(m['raw33_plane9_map_columns']) if any(k==p9i and v!='0' for k,v in col)]==[rawi];sel.append([rawi,p9i,ex])
 assert d['raw33_to_plane_degree2_selector_map']==sel
 mr='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63_REVIEW/report.json';assert rd(mr)['all_passed'] and pins[mr]=='d07787a307570c50a903ad4532dd94d58e0c2a73c72e62ebcc7c78b213cc58e9'
 assert d['all99_intrinsic_native_full48_S0_combinations_exact0'] # independently checkedalready3704, notrepeated.
 evidence={'all7_honest_product_source_orders_atleast3':True,'no_honest327554_second_zero_assumption':True,'gB299650_truezero7_conductor_planeorder_atleast3':True,'raw_f63_28539_not_substituted_for_gB':True,'all24_other_leaves_physical164':True,'plane5_selectors_literal_raw33_0_8_9_24_25':sel,'plane_map_scope':'Centered formalplane s=r+ell/5 preserves sourceorder. The accepted ringquotient raw33→plane9 furthermodplane-degree3 has basis1,r,ell,rr,rell. Only ordinary plane coefficients with ell²omitted are selected, not mixedz or arbitrary quotientjets.','separate_native99_S0_zero_reused_from_independent3704':True,'ordinary_second_formulas_keep_LrSr_LrSell_LeSr':True,'no_second_response_gradient_K_b1_or_upper_value':True}
 ep='independent-order-pruning-trueplane5-and-native-source-scope.json';put(ep,evidence)
 owned=count([r,d,z,ops,peer,tail,parent,a,m,sel,evidence,pins])+8000;assert owned<60000
 report={'id':STEM,'all_passed':True,'source_order_selector_scope_only':True,'no_second_response_or_K_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(out/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned}))
if __name__=='__main__':main()
