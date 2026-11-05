#!/usr/bin/env python3
"""Prepared formal coefficient fold only; no function/cross/native response."""
import argparse,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63'
def main(out,peer):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_GATE_20261105_v63/report.json';g=json.loads(gp.read_text());pp=ROOT/peer;p=json.loads(pp.read_text());assert g['all_passed'] and p['all_passed'] and sha(gp) in json.dumps(p)
 ep=gp.parent/'factored619-old-honest-second-linear-input-and-positive-bound.json';assert sha(ep)==g['evidence_sha256'][ep.name];d=json.loads(ep.read_text())
 assert g['no_signed_weighted_gradient_or_callback_values'] and d['no_signed_weighted_gradient_or_callback_value']
 assert count([g,p])<800 and 2*count(d)+10000<g['future_formal_numeric_bound']==31822
 for _,h in g['input_sha256'].items():assert sha(ROOT/h['path'])==h['sha256']
 cap=d['positive_prior_bits'];assert cap==46 and d['temporary_prior_bits']==139
 (out/'phase.json').write_text(json.dumps({'phase':'authorized_FORMAL_oldsecond_weights_only','gate_sha256':sha(gp),'source_peer_sha256':sha(pp),'input_numeric_including_appended_metadata':count(d),'decoded_raw_overlap_plus10000':2*count(d)+10000,'no_function_values':True},separators=(',',':'))+'\n')
 maximum=0
 def guard(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=cap
  return x
 grad={}
 for rec in d['all619_four_channel_old_second_factor_records']:
  key=(rec['case619'],rec['channel4'],rec['ordinary_pair8'],rec['original_C0_symbol']);q=Q(1)
  for x in rec['rational_factors']:q=guard(q*Q(x))
  v=guard(grad.get(key,Q(0))+q)
  if v:grad[key]=v
  else:grad.pop(key,None)
 for key,ds,us in d['all_positive_symbol_coefficient_envelopes']:
  v=grad.get(tuple(key),Q(0));D,U=int(ds),int(us);assert D%v.denominator==0 and abs(v.numerator)*(D//v.denominator)<=U
 labels=sorted({k[3] for k in grad});nonzero_cases=sorted({k[0] for k in grad});nonzero_channels=sorted({(k[0],k[1]) for k in grad})
 result={'exact_nonzero_old_second_symbol_coefficients':[[list(k),str(v)] for k,v in sorted(grad.items())],
 'nonzero_original_C0_symbols':labels,'nonzero_case619_indices':nonzero_cases,'nonzero_case_channel_pairs':[list(x) for x in nonzero_channels],
 'all_positive_envelopes_checked':True,'all619_and_all_four_channels_traversed':True,'ordinary_second_pair_order8':d['ordinary_second_pair_order8'],
 'no_actual_honestsecond_native_gauge_firstcross_or_pivotcarry_response':True,
 'scope':'Onlyformaloldsecondlinearcoefficients. Absentslotsarezeroweights,notzerorelationvalues. Firstcrossproducts/native/pivotcarryandactualcachebindingremainseparategates.'}
 assert count([d,grad,result])+3000<g['future_formal_numeric_bound']
 path=out/'actual-formal619-oldsecond-weights-and-symbol-support.json';path.write_text(json.dumps(result,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'formal_nonzero_symbol_slots':len(grad),'formal_nonzero_case_channels':len(nonzero_channels),'formal_nonzero_cases':len(nonzero_cases),'nonzero_original_symbols':labels,'no_unknown_function_response_values':True,
 'source_prior_bits':cap,'actual_bits':maximum,'input_sha256':{str(gp.relative_to(ROOT)):sha(gp),str(pp.relative_to(ROOT)):sha(pp),str(ep.relative_to(ROOT)):sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir())<g['future_RAW_bound'];print(json.dumps({'all_passed':True,'formal_nonzero_slots':len(grad),'symbols':labels,'actual_bits':maximum,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.source_peer)
