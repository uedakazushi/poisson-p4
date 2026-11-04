#!/usr/bin/env python3
"""Source-only resource frontier. No global matrix/minor/field computation."""
import hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='S61_NONSQUARE_P100_LOW98_POLYNOMIAL_COST_FRONTIER_20261104_v63';O=R/'reports'/S
paths=['reports/S61_NONSQUARE_P100_LOW98_FIXED_SOURCE_PREFLIGHT_20261104_v63/source-inventory.json','reports/S61_NONSQUARE_P100_PURE70_SOURCE_PREFLIGHT_20261104_v63/certificate.json']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 O.mkdir(exist_ok=True);a,b=[json.loads((R/p).read_text()) for p in paths]
 f=b['operator_support']['seven_numeric98x84_operator_forecasts'];total=sum(x['all_numeric_operator_RF_terms_upper'] for x in f)
 assert b['all70_outputs_if_retained_RF_upper']==454974 and b['max_one_row_live_upper']==56407 and total==134032
 e={'actual_global98_matrix_evaluated':False,'all70_kernel_products_retained_forecast':454974,'cap':150000,'straight_resident_forecast_passes':False,'max_streamed_one_row_forecast':56407,'one_row_generation_bound_passes':True,'seven_materialized_multiplication_operator_forecasts_total':total,'raw98_det_state_bound_not_certified':True,'possible_next_method':'Compute one actual low-source global determinant D(z), with known curve unit denominator and prime-integral monic numerator; if degree and every old-prime factor are bounded, a finite quotient CRT/Nakayama test could close the entire candidate support. No such determinant or factor bound has been computed.','no_claim_all_algorithms_infeasible':True,'global239_still_open':True}
 p=O/'resource-frontier.json';p.write_text(json.dumps(e,indent=2)+'\n');report={'id':S,'source_only':True,'all_metadata_checks_passed':True,'no_arithmetic_coefficient_rank_GB_point_or_newprime_job':True,'input_sha256':{q:sha(R/q) for q in paths},'evidence_sha256':{p.name:sha(p)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md')};(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
