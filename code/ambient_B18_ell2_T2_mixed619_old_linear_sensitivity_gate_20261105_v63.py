#!/usr/bin/env python3
"""Factored formal sensitivity input, no weighted signed gradient or values."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_GATE_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 alias=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_ALIAS_ORIGIN_20261105_v63/report.json');assert alias['all_passed'] and alias['all62_leaf_origins_zero']
 prof=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_CALLBACK_PROFILE_20261105_v63/literal-oldB-native-X2-callback-and-firstcarry-profile.json')
 m=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes'])
 B=list(map(dict,m['old640_source_columns']));ops=dict(m['old640_typed_closure_nodes'])
 g=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json')
 C=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations']
 zero={};linear={}
 for k in prof['reachable_typed_node_ids']:
  n=ops[k];op=n[0]
  if op=='actual_C0_scalar_function':zero[k]=True;terms=[[n[1],[]]]
  elif op=='c':zero[k]=n[1]==0;terms=[]
  elif op=='sum':
   zero[k]=all(zero[j] for j,a,b in n[1]);terms=[]
   for j,a,b in n[1]:terms.extend([[label,fac+[str(Q(a,b))]] for label,fac in linear[j]])
  else:
   assert op=='product';zero[k]=any(zero[j] for j in n[3]);terms=[]
   for pos,j in enumerate(n[3]):
    others=[v for q,v in enumerate(n[3]) if q!=pos]
    if any(zero[v] for v in others):continue
    assert all(ops[v][0]=='c' for v in others)
    factors=[str(Q(n[1],n[2]))]+[str(Q(ops[v][1],ops[v][2])) for v in others]
    terms.extend([[label,fac+factors] for label,fac in linear[j]])
  linear[k]=terms
  assert op=='c' or zero[k],('nonconstant_origin_not_certified_zero',k)
 assert sum(bool(linear[k]) for k in prof['old_distinct_roots'])==13
 source_input=count([prof,m,g,C,linear,zero])+6000;assert source_input<60000
 # Stream the original factors. Stop source construction before exceeding a
 # separate1500-record limit; this is not a dynamic arithmetic-height claim.
 records=[]
 for case,rec in enumerate(C):
  for source,c in rec['original768_combination']:
   if source>=640:continue
   for channel,row in enumerate(g['Gamma4_on_source56_times8']):
    for flat,a in row:
     target,pair=divmod(flat,8);root=B[source].get(target)
     if root is None or root not in linear:continue
     for label,factors in linear[root]:
      records.append({'case619':case,'channel4':channel,'ordinary_pair8':pair,'original_C0_symbol':label,'literal_old640_source':source,'literal_typed_root':root,'rational_factors':[c,a]+factors})
      assert len(records)<=1500,('source_record_guard_before_weighted_AD',len(records))
 maximum=0;bound={}
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 for rec in records:
  D,U=1,1
  for x in rec['rational_factors']:
   x=Q(x);D*=x.denominator;U*=abs(x.numerator);mark((D,U))
  key=(rec['case619'],rec['channel4'],rec['ordinary_pair8'],rec['original_C0_symbol']);dd,uu=bound.get(key,(1,0));nd=math.lcm(D,dd);bound[key]=mark((nd,uu*(nd//dd)+U*(nd//D)))
 assert maximum<=384
 labels=sorted({r['original_C0_symbol'] for r in records});oldroots=sorted(k for k in prof['old_distinct_roots'] if linear[k])
 data={'all62_leaf_origin_zero_source_pinned':True,'all668_nonconstant_typed_nodes_certified_zero_origin':True,'nonzero_linear_parts_only13_raw_roots':oldroots,
 'typed_root_linear_part_factor_recipes':[[k,v] for k,v in sorted(linear.items()) if v],
 'all619_four_channel_old_second_factor_records':records,'all_independent_C0_second_symbols':labels,
 'all_positive_symbol_coefficient_envelopes':[[list(k),str(D),str(U)] for k,(D,U) in sorted(bound.items())],
 'ordinary_second_pair_order8':g['ordinary_second_pair_order8'],'scope':'FORMALoldsecondlinearcoefficientsONLY.Noactualhonestsecond/native/gauge/Phi/firstcross/carryresponsequeried. All knownfirstcrossproducts ofnonlinearrootsremainforfutureactualvalues.',
 'no_signed_weighted_gradient_or_callback_value':True,'positive_prior_bits':maximum,'temporary_prior_bits':3*maximum+1,
 'zero_case_channel_coefficients_interpretation':'Caseswithoutfactorrecords haveformallyzerooldsecondlinearcoefficient,notzerorelationvalue. Every619case/fourchannelswastraversed.'}
 future=2*count(data)+2*(len(bound)*6)+6000;sourceowned=source_input+count(data)+4000
 raw=2_000_000+future*(2*math.ceil(maximum*math.log10(2))+40)
 assert future<60000 and sourceowned<150000 and raw<32*1024*1024
 data.update({'future_formal_fold_numeric_bound':future,'future_RAW_bound':raw,'source_owned_bound':sourceowned})
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'factored619-old-honest-second-linear-input-and-positive-bound.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'formal_factor_records':len(records),'distinct_weighted_symbol_slots':len(bound),'surviving_original_symbols':labels,'positive_prior_bits':maximum,'future_formal_numeric_bound':future,'future_RAW_bound':raw,'no_signed_weighted_gradient_or_callback_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'records':len(records),'slots':len(bound),'symbols':labels,'prior':maximum,'future_numeric':future,'RAW':raw,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
