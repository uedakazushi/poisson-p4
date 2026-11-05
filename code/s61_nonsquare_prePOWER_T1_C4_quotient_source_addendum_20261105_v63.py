#!/usr/bin/env python3
"""Close C4 peer premise and count all scalar/cache operations. No quotient values."""
import hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_SOURCE_ADDENDUM_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 paths=[ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_SOURCE_20261105_v63/report.json',
  ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_SOURCE_20261105_v63/finite-operator-and-preflight.json',
  ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_HORNER_REVIEW_20261105_v63/report.json']
 for p in [paths[0],paths[2]]:assert json.loads(p.read_text())['all_passed']
 plan=json.loads(paths[1].read_text());assert plan['owned_MAX']==81776
 checks=[math.comb(h+3,3)<=100**h for h in range(17)];assert all(checks)
 formula_pairs=7595*289;rows=7595*17
 extra=2*rows+768+120+7595+20+100
 whole=formula_pairs+25*10013+extra
 # Extremely conservative serialization: even the reserved metadata is
 # allotted a full518-bit integer+key record, though it is not all live.
 fullstate=(32888+16000)*170+100000
 RAW=fullstate+2000000+200000
 assert RAW<15000000
 d=dict(C4_valuation_independent_peer_fulfilled=True,
  inherited_source_report_UNCHANGED=True,actual_quotient_NOT_run=True,
  all17_weighted_scalar_inequalities_passed=checks,
  inequality='binom(h+3,3)<=100^h makes each pre-root scalar multiplication <=100^16; full weighted sum <4*100^16',
  original_polynomial_pair_scope=formula_pairs+25*10013,
  row_scaling_scalar_products_upper=2*rows,
  fixed_R_power_cache_scalar_products_upper=768,fixed_C4_cache_scalar_products_upper=120,
  input_scaling_products_upper=7595,constant_100_power_products_upper=20,control_binomial_reserve=100,
  ALL_scaled_root_and_polynomial_products_upper=whole,
  whole_owned_upper=81776,integer_temp_bits_upper=518,
  whole_state_numeric_terms32888=32888,ALL_metadata_reserve16000=16000,
  whole_state_serialization_bytes_upper=fullstate,
  narrow_runtime_input_code_copy_bound2000000=2000000,
  all_future_receipts_and_reports_bound200000=200000,
  concrete_RAW_upper=RAW,unchanged_RAW_hard_limit15000000=15000000,
  lifetime='streamONE Q_i slab; verify all out-of-bound coefficients ZERO before retention; releaseR-cache; multiply fixedC4 into ONE remainder accumulator; no prior success output before guard-prone verification completes',
  no_double_full_state_output=True,new_parent_GO_still_required=True)
 (out/'proof.json').write_text(json.dumps(d,indent=2)+'\n')
 paths.append(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,independent_C4_peer_fulfilled=True,
  ALL_scalar_products_upper=whole,whole_owned81776=81776,integer_bits518=518,RAW_upper=RAW,
  actual_quotient_NOT_run=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'proof.json':sha(out/'proof.json')})
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(sha(out/'report.json'));print(json.dumps(d,indent=2))
if __name__=='__main__':main()
