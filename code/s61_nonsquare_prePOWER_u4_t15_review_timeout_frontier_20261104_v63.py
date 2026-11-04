#!/usr/bin/env python3
"""Record the stopped independent check; does not launch arithmetic."""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_PREPOWER_DIRECT_SCHUR_ENTRY_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def main():
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists()
 work=ROOT/'work/prePOWER-direct-one-entry-independent-review'
 assert (work/'native.log').read_text()=='48 second guard\n'
 phase=json.loads((work/'phase.json').read_text());assert phase['phase']=='independent source/type/input checks complete; before arithmetic'
 assert sha(str((work/'input.txt').relative_to(ROOT)))==phase['input_sha256']=='b8b29b2e9239548811f64e45d828b49d0e3a6e15969d369422e56e1ecd5a3e74'
 for src,dst in [('native.log','native-stop.log'),('failed-checker.py','failed-checker.py'),('failed-checker.cpp','failed-checker.cpp')]:shutil.copyfile(work/src,out/dst)
 paths=['reports/S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63/report.json','reports/S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63/report.json','reports/S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW/report.json','research/notes/'+STEM+'.md',str(Path(__file__).relative_to(ROOT))]
 result={'id':STEM,'all_passed':False,'status':'STOPPED_independent_fixed_arithmetic_at48s_guard_no_retry',
  'source_type_preflight_completed':True,'typed_chart_peer_separately_final':True,
  'one_entry7022_QQ_coefficients_accepted':False,'partial_numeric_identity_acceptance':False,
  'arithmetic_stage_receipt_before_stop':'none; do not infer whether the source DAG or border determinant completed',
  'resource_bound_increased_or_extra_worker_launched':False,'no_other_entry_minor_point_prime_or_search':True,
  'ignored_input':{'path':str((work/'input.txt').relative_to(ROOT)),'sha256':phase['input_sha256'],'bytes':(work/'input.txt').stat().st_size},
  'input_sha256':{p:sha(p) for p in paths},
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('input_sha256','evidence_sha256')}))
if __name__=='__main__':main()
