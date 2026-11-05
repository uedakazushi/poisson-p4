#!/usr/bin/env python3
"""Immutable worker/helper/source custody and safe --help, no arithmetic."""
import argparse,ast,hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CUSTODY_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 worker=ROOT/'code/ambient_B18_ell2_T2_active4_cache_bundle_20261105_v63.py';seen={}
 def visit(p):
  key=str(p.relative_to(ROOT))
  if key in seen:return
  tree=ast.parse(p.read_text());seen[key]=sha(p)
  for n in ast.walk(tree):
   names=([n.module] if isinstance(n,ast.ImportFrom) and n.module else [x.name for x in n.names] if isinstance(n,ast.Import) else [])
   for name in names:
    q=ROOT/'code'/(name.replace('.','/')+'.py')
    if q.exists():visit(q)
 visit(worker);before=dict(seen)
 gate='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/report.json';g=json.loads((ROOT/gate).read_text());assert g['all_passed']
 assert sha(ROOT/gate)=='e196d39bab86c6f15073bd655eec3a49f4947bd2ebf6c2dc329fa63336bf3672'
 old='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_INPUT_GATE_20261105_v63/report.json';o=json.loads((ROOT/old).read_text());assert sha(ROOT/old)=='7b6a485ad19ffee2002c2c099133e56e4a4b7609f0179c3bd2c3f180874bb88b'
 assert sha(ROOT/'code/ambient_B18_ell2_T2_active4_cache_bundle_input_20261105_v63.py')==o['code_sha256']
 err='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_WHOLE_HEIGHT_ERRATUM_20261105_v63/report.json';assert sha(ROOT/err)=='be755cc7b5c5ac2b3131231274b4a5f1f58a5170efa50477d331e7b86f6b95f6'
 env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
 result=subprocess.run([sys.executable,str(worker),'--help'],cwd=ROOT,env=env,text=True,capture_output=True,timeout=10)
 assert result.returncode==0 and not result.stderr
 assert all(x in result.stdout for x in ['usage:','--output','--gate','--source-peer','--custody'])
 assert all(sha(ROOT/p)==h for p,h in seen.items())
 hp=out/'safe-help.stdout.txt';hp.write_text(result.stdout)
 report={'id':STEM,'all_passed':True,'safe_help_literal_PASS':True,'safe_help_exit_status':result.returncode,
  'worker_path':str(worker.relative_to(ROOT)),'worker_sha256':sha(worker),'local_worker_and_helper_sha256':seen,
  'worker_note_sha256':sha(ROOT/'research/notes/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63.md'),
  'corrected_gate_path':gate,'corrected_gate_sha256':sha(ROOT/gate),'corrected_gate_producer_code_sha256':g['code_sha256'],
  'literal_corrected_bounds_PASS':g['whole_owned_numeric_bound']==22891 and g['whole_raw_byte_bound']==5131577 and g['whole_value_height_bound_bits']==948,
  'source_temporary_bound_1231_under2048':True,'frozen_withdrawn_old_gate':{'path':old,'sha256':sha(ROOT/old)},
  'withdrawal_erratum':{'path':err,'sha256':sha(ROOT/err)},'safe_help_stdout_sha256':sha(hp),
  'safe_help_had_no_arithmetic_input_arguments':True,'actual_value_worker_executed':False,'fresh_corrected_SOURCE_peer_still_separately_required':True,
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 assert report['literal_corrected_bounds_PASS']
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','safe_help_literal_PASS','worker_sha256','corrected_gate_sha256','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
