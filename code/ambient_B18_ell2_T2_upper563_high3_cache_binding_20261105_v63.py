#!/usr/bin/env python3
"""Literal selector/frame binding ONLY; do not extract any third coefficient."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_CACHE_BINDING_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);hs={}
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 p8=read('reports/AMBIENT_B18_ELL2_T2_UPPER563_P8_ANCESTOR_SOURCE_20261105_v63/literal-two-P8-composite-program-and-exact-boundary-catalog.json');assert dict(p8['selected_original_C0_arithmetic_nodes'])=={327602:['n',40892],327603:['a',20114,327602],327626:['n',40916],327627:['a',20138,327626]}
 phys='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json';t=field(phys,'physical12_coefficients');records=[]
 for idx,row,a,b,node in [(0,18,20114,40892,327603),(4,60,20138,40916,327627)]:
  v=t[idx];assert v['physical_row']==row and v['roots']['aA_node']==a and v['roots']['aB_node']==b and v['roots']['delta_operator']==['subtract',a,b]
  records.append({'literal_C0_function':node,'physical_row':row,'physical12_record_index':idx,'field':'delta','source_roots':[a,b],'coefficient_field':'R0 scalar physical coefficient, SAME trace-free frame','raw33_coefficient_selector':30})
 raw='reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json';basis=field(raw,'actual33_basis');assert basis[30]==[0,2,1]
 qp='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';mon=field(qp,'plane9_monomials');pi=field(qp,'raw33_plane9_map_columns');assert mon[6]==[2,1] and pi[30]==[[6,'1']]
 assert [i for i,v in enumerate(pi) if any(j==6 and x!='0' for j,x in v)]==[30]
 prefix='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';b8=field(prefix,'actual_C0_basis8_scalar_nodes');div=field(prefix,'all35_regular_coordinate_divisions');assert b8[7]==505938 and div[0]['original_pure_CA_source_function']==7790 and div[0]['restriction_r_s_node']==505938
 vals='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';catalog=field(vals,'all_requested164_root_values');ids={i for i,v in catalog};assert {20114,40892,20138,40916,7790}<=ids
 # No coefficient at index30 (or any other nonconstant index) is read here.
 del catalog
 premises={
 'physical105_and12_source':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33_to_plane9':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63_REVIEW/report.json',
 'prefix0_unit57_clone':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'prefix_and_164_readback_values':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'literal_P8_AST':'reports/AMBIENT_B18_ELL2_T2_UPPER563_P8_ANCESTOR_SOURCE_20261105_v63_REVIEW/report.json'}
 for p in premises.values():assert read(p)['all_passed']
 ep=out/'exact-three-direct-high3-existing-cache-source-selectors.json';data={'P8_delta_source_records':records,'gA_source_record':{'literal_C0_function':505938,'actual_zero7_prefix_original_function':7790,'equality_in_plane_quotient':'7790-505938 in(zA0,...,zA6) byactualregularcoordinateunit57prefix','physical164_member':True,'raw33_coefficient_selector':30},'actual_raw33_basis_index30':[0,2,1],'actual_plane9_output_coordinate6':[2,1],'all_other_raw33_basis_columns_have_zero_output_coordinate6':True,'ordinary_Taylor_coefficient_no_factorial':True,'source_premises':premises,'numerical_third_coefficients_NOT_extracted':True,'actual_two_linear_combinations_or_upper_response_NOT_evaluated':True,'all_other_upper_carries_still_retained':True};H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'source_selector_and_frame_binding_only':True,'actual_third_or_upper_values_queried':False,'source_premises':premises,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'owned_numeric_bound':count([p8,t,basis,mon,pi,b8,div,records,data])+10000,'seconds':time.monotonic()-start};assert report['owned_numeric_bound']<60000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_selector_and_frame_binding_only','actual_third_or_upper_values_queried','owned_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
