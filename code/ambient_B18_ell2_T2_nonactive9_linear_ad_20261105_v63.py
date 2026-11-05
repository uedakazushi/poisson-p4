#!/usr/bin/env python3
"""Prepared formal linear old-second/source fold only; no function values."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_20261105_v63'
def main(out,gate,peer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/gate;g=json.loads(gp.read_text());pp=ROOT/peer;p=json.loads(pp.read_text());assert g['all_passed'] and p['all_passed'] and H.sha(gp) in json.dumps(p)
 assert not g['signed_sensitivities_or_values_evaluated'];ip=gp.parent/'eighteen-factored-old-linear-second-and-native-source-inputs.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 for path,h in g['input_sha256'].items():assert H.sha(ROOT/path)==h
 assert count([g,p])<600 and 2*count(d)+4000<g['future_whole_numeric_bound']<60000
 cap=g['positive_prior_bits'];assert cap<=384 and g['temporary_prior_bits']==3*cap+1<=1153
 H.dump(out/'phase.json',{'phase':'authorized_formal_linear_sensitivity_only','source_gate_sha256':H.sha(gp),'source_peer_sha256':H.sha(pp),'no_function_values':True})
 maximum=0
 def guard(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=cap;return x
 results=[]
 for r in d['all18_factored_linear_second_sensitivity_inputs']:
  grad={};native=[{},{}]
  for label,factors in r['old_literal_second_terms_as_four_factor_products']:
   q=F(1)
   for x in factors:q=guard(q*F(x))
   v=guard(grad.get(label,F(0))+q)
   if v:grad[label]=v
   else:grad.pop(label,None)
  for ch,index,factors in r['native128_terms_as_two_factor_products']:
   q=F(1)
   for x in factors:q=guard(q*F(x))
   v=guard(native[ch].get(index,F(0))+q)
   if v:native[ch][index]=v
   else:native[ch].pop(index,None)
  for label,ds,us in r['positive_old_symbol_coefficient_D_U']:
   x=grad.get(label,F(0));D,U=int(ds),int(us);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  for (ch,index),ds,us in r['positive_native_coefficient_D_U']:
   x=native[ch].get(index,F(0));D,U=int(ds),int(us);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  results.append({'original418_source':r['original418_source'],'contraction':r['contraction'],'exact_independent_original_C0_second_sensitivity':srec(grad),
   'exact_original_native128_combinations_two_channels':[srec(v) for v in native],
   'native_S0_zero_reason':'linearcombinationofindependentlyverifiedselected99S0zeros; literalinputfactorrecipe preserved',
   'M1_pivot_G2_and_known_lower_terms_remain_unevaluated':True})
 assert count([d,results])+3000<g['future_whole_numeric_bound']
 ep=out/'actual-eighteen-formal-old-second-sensitivities-and-native-source-combinations.json';H.dump(ep,{'all18_formal_source_results':results,'all_old_symbol_positive_envelopes_verified':True,'old_second_function_values_not_queried':True,'native_S2_L1_values_not_queried':True,
  'M1_pivot_G2_remaining_formal':True,'actual_T2_HF2_or_length_not_computed':True,
  'native_eligibility_scope':'weightedsourcecombinationsonly,notindividualeffectivepivots; L1rOmegaEll/L1ellOmegaRstillrequiredbylatervaluegate'})
 nz=[(r['original418_source'],r['contraction']) for r in results if r['exact_independent_original_C0_second_sensitivity']]
 report={'id':STEM,'all_passed':True,'formal_old_second_nonzero_contraction_count':len(nz),'all_formal_old_second_gradients_zero':not nz,'no_physical_gauge_implicit_or_honest_second_value_evaluated':True,'actual_HF2_or_length_not_computed':True,
  'input_sha256':{gate:H.sha(gp),peer:H.sha(pp),str(ip.relative_to(ROOT)):H.sha(ip)},'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'actual_maximum_bits':maximum,'source_prior_bits':cap,'future_whole_numeric_bound':g['future_whole_numeric_bound'],'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);assert sum(v.stat().st_size for v in out.iterdir())<g['future_RAW_bound'];print(json.dumps({k:report[k] for k in ['all_passed','formal_old_second_nonzero_contraction_count','all_formal_old_second_gradients_zero','actual_maximum_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
