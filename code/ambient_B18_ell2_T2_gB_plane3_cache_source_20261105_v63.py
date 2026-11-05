#!/usr/bin/env python3
"""Allorder gB tail/sourceorders/selectors only; no source2/3 extraction."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 chain='reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/';assert read(chain+'report.json')['all_passed']
 tail=next(v for v in field(chain+'actual-UP-QB-640-typed-operator.json','accepted_actual_regular_function_identities') if v['kind']=='accepted_actual_CB_f63_tail');assert tail['function']==28539 and tail['gB']==299650
 rel='reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json';Y=field(rel,'CB_coordinate_Y9_nodes');h=field(rel,'CB_f63_zero7_telescoping_coefficients');assert Y[:7]==[27897,27948,28203,28254,28437,28692,28794] and h==tail['coordinate_quotient_nodes']==[327551,327554,327557,327560,327563,327566,327569]
 assert field(rel,'CB_f63_node')==28539 and field(rel,'CB_g_restriction_node')==299650
 div={v['node']:v for v in field(rel,'all_coordinate_divisions')}
 for z,q in zip(Y[:7],h,strict=True):assert div[q]['coordinate']==z and 'paired unique unit57' in div[q]['proof']
 fc=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640');first={i:H.parse(v) for i,v in fc['all79_atom_values_full10']};rows=[]
 for z,q in zip(Y[:7],h,strict=True):
  assert z in first and q in first and all(not first[z].get(i) for i in [0,8,9]) and not first[q].get(0)
  rows.append({'zB_node':z,'h_node':q,'zB_point_r_ell_exact0':True,'h_point_exact0':True,'known_h_r_first':str(first[q].get(8,0)),'known_h_ell_first':str(first[q].get(9,0))})
 raw='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';ids={i for i,v in field(raw,'all_requested164_root_values')};assert {28539,*Y[:7]}<=ids
 basis=field('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json','actual33_basis');qp='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';mon=field(qp,'plane9_monomials');pi=field(qp,'raw33_plane9_map_columns')
 selectors=[]
 for rawi,p9i,m in [(24,3,[2,0]),(25,4,[1,1]),(30,6,[2,1])]:
  assert basis[rawi]==[0,*m] and mon[p9i]==m and pi[rawi]==[[p9i,'1']] and [i for i,v in enumerate(pi) if any(j==p9i and x!='0' for j,x in v)]==[rawi]
  selectors.append([m,rawi,p9i])
 premises=['reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json']
 for p in premises:assert read(p)['all_passed']
 cp=read(premises[0]);assert cp['input_sha256'][chain+'report.json']==hs[chain+'report.json'] and cp['input_sha256'][chain+'actual-UP-QB-640-typed-operator.json']==hs[chain+'actual-UP-QB-640-typed-operator.json']
 ep=out/'literal-gB-tail-orders-and-existing-source3-source2-selectors.json';data={'accepted_chain_stem':'NEGATIVE_BOUNDARY','source_allorder_identity':'28539-299650=sum zB_i*h_i withliteralpairedunit57telescopes','ordered_source_rows':rows,'fCB_literal164_node':28539,'zB7_literal164_nodes':Y[:7],'ordinary_gB_r2ell_formula':'fCB[30]-sum_i(zB_i[24]*h_i_ell+zB_i[25]*h_i_r)','source_selector_triplets':selectors,'zB3_h0_and_all_h2_h3_terms_zero_by_source_order':True,'no_higher_source_coefficient_numbers_extracted_or_multiplied':True,'no_actual_gB3_or_OLD_or_upper_response':True,'source_premises':premises};H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'allorder_source_orders_and_cache_selectors_only':True,'new_source2_3_or_gauge_or_implicit_coefficients':False,'actual_gB3_or_relation_values_queried':False,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'owned_numeric_bound':count([tail,Y,h,div,first,basis,mon,pi,rows,data])+12000,'seconds':time.monotonic()-start};assert report['owned_numeric_bound']<60000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_gB3_or_relation_values_queried','owned_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
