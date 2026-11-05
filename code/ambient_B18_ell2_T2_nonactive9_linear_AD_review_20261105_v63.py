#!/usr/bin/env python3
"""Independent formal factor sums; no unknown callback response is read."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and r['no_physical_gauge_implicit_or_honest_second_value_evaluated']
 assert r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_linear_ad_20261105_v63.py')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 source='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_INPUT_GATE_V2_20261105_v63'
 own=read(source+'_REVIEW/report.json');assert own['all_passed']
 d=read(source+'/eighteen-factored-old-linear-second-and-native-source-inputs.json')
 ep=f'reports/{BASE}/actual-eighteen-formal-old-second-sensitivities-and-native-source-combinations.json';v=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 assert count(d)==2572 and 2*count(d)+4000==9144<11726
 assert len(v['all18_formal_source_results'])==18;nonzero=[]
 for rec,result in zip(d['all18_factored_linear_second_sensitivity_inputs'],v['all18_formal_source_results'],strict=True):
  assert (rec['original418_source'],rec['contraction'])==(result['original418_source'],result['contraction'])
  # Group by exact symbolic slot first, then sum its rational products. No
  # value is assigned even to literal label0 or any coordinate label.
  old={};native=[{},{}]
  for label,factors in rec['old_literal_second_terms_as_four_factor_products']:
   q=Q(1)
   for x in factors:q*=Q(x)
   old.setdefault(label,[]).append(q)
  grad={label:sum(xs,Q(0)) for label,xs in old.items() if sum(xs,Q(0))}
  for ch,idx,factors in rec['native128_terms_as_two_factor_products']:
   q=Q(1)
   for x in factors:q*=Q(x)
   native[ch][idx]=native[ch].get(idx,Q(0))+q
  native=[{i:x for i,x in vv.items() if x} for vv in native]
  assert grad==sp(result['exact_independent_original_C0_second_sensitivity'])
  assert native==list(map(sp,result['exact_original_native128_combinations_two_channels']))
  for label,ds,us in rec['positive_old_symbol_coefficient_D_U']:
   x=grad.get(label,Q(0));D,U=int(ds),int(us);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  for (ch,idx),ds,us in rec['positive_native_coefficient_D_U']:
   x=native[ch].get(idx,Q(0));D,U=int(ds),int(us);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  assert result['M1_pivot_G2_and_known_lower_terms_remain_unevaluated']
  if grad:nonzero.append([rec['original418_source'],rec['contraction']])
 assert len(nonzero)==r['formal_old_second_nonzero_contraction_count']==12
 assert not r['all_formal_old_second_gradients_zero']
 assert v['old_second_function_values_not_queried'] and v['native_S2_L1_values_not_queried'] and v['M1_pivot_G2_remaining_formal']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent18-symbolic-factor-sums-and-positive-envelope-checks.json';e.write_text(json.dumps({'all18_old_and_native_formal_sums_verified':True,'nonzero_old_gradient_contractions':nonzero,'nonzero_count':12,'zero_count':6,'all10_symbols_retained_independent':True,'no_honestsecond_or_native_response_values':True,'scope':'Formal source weights only. Nonzero symbolic coefficients are not evidence of unavailable physical data or nonzero relations; M1/pivot/G2/canonicalframe terms remain for a separate source/cache gate.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'formal_nonzero_count':12,'formal_zero_count':6,'no_unknown_function_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'formal_nonzero':12}))
if __name__=='__main__':main()
