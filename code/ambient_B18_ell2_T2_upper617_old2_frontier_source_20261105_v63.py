#!/usr/bin/env python3
"""Literal opaque OLD2 frontier/cache inventory only, no response values."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(exist_ok=False,parents=True);pins={};rawreads=0
 def read(p):
  nonlocal rawreads
  pins[p]=H.sha(ROOT/p);rawreads+=2*(ROOT/p).stat().st_size;assert rawreads<30*1024*1024
  return json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/'
 r=read(base+'report.json');assert r['all_passed'] and pins[base+'report.json']=='643f72906d5db1b74cbf100ed401160336d8c94619f3fd35f5b7d4e8872415e0'
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json')['all_passed']
 x=read(base+'actual617-formal-combined-OLD2-second-symbol-coefficients.json')
 supports={}
 for case,channels in x['all617_combined_OLD2_independent_symbol_coefficients']:
  for ch,v in enumerate(channels):
   for p,f,q in v:
    assert q!='0';supports.setdefault(f,set()).add(p)
 symbols=sorted(f for f in supports if f>=299000 and f not in [299650,505938])
 oldg=[330164,330289,330524,330594];oldh=[330184,330309,330544,330614];newg=list(range(566400,566404))
 assert symbols==sorted(oldg+oldh+newg) and all(supports[f]=={8,9} for f in symbols)
 invp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';inv=read(invp)
 old4=inv['outer_four_pivot_source']['original4_old_pivots'];assert old4==[604,613,622,631]
 assert all(not any(j in old4 for j,q in c) for c in inv['original619_c0_source_combinations'])
 # Expand arithmetic only. r/g are retained as the actual regular quotient
 # constructors, never replaced by division in an Artin coefficient ring.
 cp=H.PATHS['C0'];pins[cp]=H.sha(ROOT/cp);rawreads+=2*(ROOT/cp).stat().st_size;assert rawreads<30*1024*1024
 text=(ROOT/cp).read_text();begin=text.index('[',text.index('"DAG_scalar_nodes"'))+1;dec=json.JSONDecoder();front=set(symbols);nodes={};boundaries={};layers=[]
 known79={i for i,k,q in x['all71_literal_function_classifications_and_zero_origins'] if k!='source164'}-{*oldg,*oldh,*newg}
 # The old f63 telescope coefficients are independent accepted79 boundaries.
 known79.update([327551,327554,327557,327560,327563,327566,327569])
 while front:
  want=set(front);front=set();found={};pos=begin
  for i in range(max(want)+1):
   while text[pos].isspace() or text[pos]==',':pos+=1
   n,pos=dec.raw_decode(text,pos)
   if i in want:found[i]=n
  assert set(found)==want
  for i,n in found.items():
   if i in known79:boundaries[i]=n;continue
   nodes[i]=n
   if n[0] in ['r','g']:continue
   assert n[0] in ['a','m','n','c']
   kids=n[1:3] if n[0] in ['a','m'] else n[1:2] if n[0]=='n' else []
   for j in kids:
    assert j<i
    if j not in nodes and j not in boundaries:front.add(j)
  layers.append(sorted(want));assert len(nodes)<200 and len(layers)<20
 del text
 quotientnodes=[[i,n] for i,n in sorted(nodes.items()) if n[0] in ['r','g']]
 assert all(nodes[i][0]=='g' and nodes[i][2]==299650 for i in oldg)
 assert all(nodes[i][0]=='g' and nodes[i][2]==505938 for i in newg)
 rawr=[i for i,n in quotientnodes if n[0]=='r'];assert len(rawr)==4 and all(nodes[i][2]==27948 for i in rawr)
 # Bind the *presence and labels* of accepted rr values, not their numerics.
 ubase='reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/'
 ur=read(ubase+'report.json');assert ur['all_passed']
 up=ubase+'actual-gdiff-second-and-complete619-pure-cubic-values.json';u=read(up)
 accepted_rr_ids=[v[0] for v in u['actual_eight_regular_g_second_values']]
 assert accepted_rr_ids==oldg+newg
 del u
 upr=read('reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json');assert upr['all_passed']
 orderp='reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json';assert read(orderp)['all_passed']
 ev={'all42_surviving_formal_symbol_pair_support':[[f,sorted(p)] for f,p in sorted(supports.items())],'physical164_surviving_ids28':sorted(f for f in supports if f<299000),'conductors_second_order_zero_by_accepted_order3_source':[299650,505938],'opaque_remaining_functions12':symbols,'opaque24_scalar_symbols':[[f,p] for f in symbols for p in [8,9]],'opaque_labels_occurring_only_in_outer_old4':True,'all619_c0_oldpivot_coefficients_literally_absent':old4,'literal_original_C0_arithmetic_and_regular_quotient_nodes':sorted(nodes.items()),'literal_quotient_frontier':quotientnodes,'literal_higher_arithmetic_boundaries':sorted(boundaries.items()),'ancestry_layers':layers,'accepted_existing_rr_only_cache_binding':{'report':ubase+'report.json','data':up,'field':'actual_eight_regular_g_second_values','functionIDs':accepted_rr_ids,'pair_index10':8,'no_numerical_rr_readback_copied_or_contracted':True},'remaining_individual_scalar_frontier_after_rr_cache_presence':[[f,p] for f in symbols for p in [8,9] if not(p==8 and f in accepted_rr_ids)],'outer_response_functional_recipe':{'for_each_oldpivot_k_and_each_pair_rr_re':'Gamma4*restored_h1_k + (Gamma28+Gamma55)*oldg_k - Gamma55*newg_k','same_two_channels_not_individual_new_coeff_requests':True,'coupling_to_617':'multiply_by_MINUS_actual_corrected_a1_ell_for_rr_or_a1_r_for_re;same563609sourcecorrection','do_not_assume_independent_symbols_are_unavoidable_new_jets':True},'future_source_questions':['Expand each restored_h1 actual arithmetic identity, preserving its raw r-cofactor and known first-cross restoration.','Bind each r/g constructor to original unit57/zero7 telescope source coordinates before a source-order/request claim.','Retain rr/rell as ordinary coefficients. Existing eight g rr values are not mixed-rell values.','Whole physical/known-first/native and actual upper values need a separate input/height/ownership gate.'],'no_new_gradient_coefficient_second_or_upper_response_values':True}
 ep=out/'literal12-opaque-function-and24-plane-second-frontier-with-existing-rr-labels.json';H.dump(ep,ev)
 owned=count([inv,x,supports,nodes,boundaries,layers,ev,pins])+10000;assert owned<60000 and rawreads+sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'metadata_ancestry_cache_presence_only':True,'actual12functions_only_rr_rell':True,'all8_existing_g_rr_labels_bound_without_values':True,'remaining_individual_scalar_labels16':True,'literal_nodes':len(nodes),'literal_r_frontiers':rawr,'literal_g_frontiers':oldg+newg,'source_owned_bound':owned,'whole_RAW_reads_upper':rawreads,'no_actual_response_or_new_gradient_values':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','literal_nodes','literal_r_frontiers','source_owned_bound','whole_RAW_reads_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
