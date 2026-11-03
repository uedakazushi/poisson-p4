#!/usr/bin/env python3
"""Freeze confirmed saved-prefix replay and an interrupted continuation."""
import argparse,gzip,hashlib,json,re,tarfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);t=time.monotonic()
 w=ROOT/'work/nonsquare-J0-pid-v63/continuation600-v1';stage=w/'replay-1';continued=w/'continued-row-gcd'
 receipt=json.loads((stage/'receipt.json').read_text());log=(stage/'solver.log').read_text();matrix=stage/'state-after-12.sing'
 assert receipt['status']=='completed'and not receipt['arithmetic_error']and receipt['exported_steps']==[4,8,11,12]
 assert 'EXACT_SAVED_12_STEP_PREFIX_REPLAY_PASS'in log and'SAVED_STATE_12_EXPORT_COMPLETE'in log and not re.search(r'(?m)^\s*\?',log)
 with matrix.open()as f:assert f.readline().strip()=='ring R=(0,e),(a),dp;matrix M[24][70];'
 assert not(w/'report.json').exists()and not(continued/'receipt.json').exists()
 more=(continued/'solver.log').read_text();assert 'EXACT_ROW_GCD_PIVOT_e 1 LEFT 24'in more and not re.search(r'(?m)^\s*\?',more)
 members=[stage/n for n in ['run.sing','solver.log','receipt.json','state-after-12.sing']]
 members +=[continued/n for n in ['run.sing','replay-library.sing','solver.log','proof.sing']]
 members +=[w/'continued-row-gcd.stdout',w.parent/'continuation600-v1.stdout']
 archive=args.output/'confirmed-state.tar.gz'
 with archive.open('wb')as stream:
  with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:
   with tarfile.open(fileobj=zipped,mode='w')as tar:
    for p in members:
     info=tar.gettarinfo(str(p),arcname=str(p.relative_to(ROOT)));info.uid=info.gid=0;info.uname=info.gname='';info.mtime=0;info.mode=0o644
     with p.open('rb')as f:tar.addfile(info,f)
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_J0_pid_resume_20261103_v63.py',ROOT/'code/s61_nonsquare_J0_pid_rowgcd_20261103_v63.py',ROOT/'reports/S61_NONSQUARE_J0_PID_20261103_v63/report.json',ROOT/'reports/S61_NONSQUARE_J0_PID_20261103_v63/manifest.json']
 report={'id':'S61_NONSQUARE_J0_PID_RESUME_20261103_V63','all_saved_prefix_replay_checks_passed':True,
  'field':'QQ(e)[a], known-unit localization a,F,Bs,C,3e+2','confirmed_saved_steps':12,'confirmed_residual_shape':[24,70],
  'confirmed_replay_and_export_seconds':receipt['seconds'],'confirmed_matrix_sha256':sha(matrix),'confirmed_matrix_bytes':matrix.stat().st_size,
  'old_gcd_searches_repeated':False,'traced_Smith_repeated':False,'budget_authorized_seconds':600,'AS_limit_bytes':2200000000,
  'continuation_status':'interrupted_before_final_receipt; process absent after runtime restore',
  'total_controller_report_exists':False,'continued_algorithm_receipt_exists':False,
  'new_logged_pivot_witnesses':1,'new_pivot_matrix_update_confirmed':False,
  'full_coker_unit_proved':False,'actual_defect_curve_proved':False,'J0_defect_support_finite_proved':False,
  'source_of_interruption':'Cloud execution environment became unavailable; after restoration no Python/Singular process and no final receipts were present. Do not reinterpret as a completed600s timeout.',
  'next_operation':'Resume directly from archived exact state-after-12.sing. Verify the one new saved witness separately before using it; avoid first12 searches and unchanged tracedSmith. Consider polynomial coefficient clearing/native arithmetic if rational-function matrix updates dominate.',
  'archive_sha256':sha(archive),'archive_bytes':archive.stat().st_size,
  'members':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':sha(p)}for p in members},
  'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-t}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_saved_prefix_replay_checks_passed','archive_bytes','confirmed_saved_steps','full_coker_unit_proved','seconds']},indent=2))
if __name__=='__main__':main()
