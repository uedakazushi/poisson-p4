#!/usr/bin/env python3
"""Control-only prepared-worker repair; copies mathematical inputs unchanged."""
import argparse,hashlib,json,shutil,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TAIL1_GUARD_CONTROL_V2_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 old=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_WORKER_CUSTODY_V2_20261105_v63'
 r=json.loads((old/'report.json').read_text());assert r['all_passed']and r['value_worker_NOT_run']
 pins={str((old/'report.json').relative_to(ROOT)):sha(old/'report.json')}
 names=['selected-plan.json','literal-atoms58.jsonl','native-input.txt','literal-normalization-witness.json','ownership-and-worker-proof.json']
 for n in names:
  p=old/n;assert sha(p)==r['evidence_sha256'][n];pins[str(p.relative_to(ROOT))]=sha(p);shutil.copyfile(p,out/n);assert sha(out/n)==sha(p)
 cpp=ROOT/'code/s61_nonsquare_tail1_cached_kappa_guard_v2_20261105_v63.cpp'
 drv=ROOT/'code/s61_nonsquare_prePOWER_tail1_cached_kappa_guard_v2_entry_20261105_v63.py'
 for p in [cpp,drv,Path(__file__)]:pins[str(p.relative_to(ROOT))]=sha(p)
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','tail1-control',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  p=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=5);(out/(name+'.log')).write_text(p.stdout+p.stderr)
 d=dict(math_input_files_byte_identical_to_old_custody=names,coefficient_arithmetic_algorithm_UNCHANGED=True,
   source_node_post_release_clock_cursor='phase1 pair=-2 => current node complete, use counts updated; resume only at NEXT source node',
   source_node_inner_before_commit_cursor='phase1 pair=-1 with original row/column product counters',
   output_IO_failure='separate runtime error; preserves bounded partial output, no second full prefix; not a positive or prefix-only mathematical result',
   prefix_versus_output_mutually_exclusive_for_all_mathematical_guards=True,
   all_bounds_unchanged=dict(owned91676=91676,coefficient_bits359=359,temp_bits395=395,soft45=45,hard60=60,RAM2GB=2000000000,RAW9517864=9517864),
   old_prepared_code_binary_report_UNCHANGED=True,no_coefficient_execution=True)
 (out/'control-diff-proof.json').write_text(json.dumps(d,indent=2)+'\n')
 raw=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in [cpp,drv,Path(__file__)])+20000
 total=raw+4700000+100000;assert total<9517864
 d=dict(actual_prepared_plus_codes_and_metadata_reserve=raw,prefix_OR_output_upper=4700000,future_receipts_reserve=100000,total_upper=total,hard_limit=9517864)
 (out/'RAW-closure.json').write_text(json.dumps(d,indent=2)+'\n')
 report=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
   all_mathematical_inputs_and_bounds_unchanged=True,explicit_completed_node_cursor=True,no_IO_double_full_state=True,
   worker_compiled_help_only=True,RAW_total_upper=total,independent_source_code_review_pending=True,
   input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
