#!/usr/bin/env python3
"""Serialization/accounting addendum only; never reruns the mu computation."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_MU_RESOURCE_ADDENDUM_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 primary=ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/report.json'
 assert sha(primary)=='4f98f4cc228369865a1b444c09b91f045322f04e37733a4fc430719297c1ce95'
 r=json.loads(primary.read_text());p=primary.parent/'phase1-actual-first640-ell-from-existing-full9-cache.json';d=json.loads(p.read_text())
 # The original parser retains its dense JSON X0 table during compiler phase.
 # Its97 nonzero terms were counted, but43008 dense source records were not.
 # Add the ENTIRE dense table conservatively without subtracting those97.
 dense=56*768;corrected=max(r['peak_owned_QQ_slot_bound'],d['source_compiler_live_slot_bound']+dense)
 assert corrected<=136665<150000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 result={'id':STEM,'all_passed':True,'primary_sha256':sha(primary),'primary_path':str(primary.relative_to(ROOT)),
  'first640_phase_sha256':sha(p),'original_reported_peak':r['peak_owned_QQ_slot_bound'],
  'omitted_dense_X0_source_records':dense,'corrected_conservative_peak_QQ_source_records':corrected,
  'frozen_source_gate_bound':136665,'hard_limit':150000,
  'scope':'accounting-only; all original source/value/report bytes remain immutable; no arithmetic rerun',
  'primary_mathematical_result_unchanged':{'actual_Y0_rank':99,'actual_mu_B_T2':6},
  'code_sha256':sha(Path(__file__))}
 (out/'report.json').write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps(result))
if __name__=='__main__':main()
