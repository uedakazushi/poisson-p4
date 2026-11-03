#!/usr/bin/env python3
"""Freeze the exact J0 preparation and honest unfinished PID continuation."""
import argparse,gzip,hashlib,json,tarfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 work=ROOT/'work/nonsquare-J0-pid-v63';old=ROOT/'work/nonsquare-v62/cubic27-v6'
 source=old/'after21';assert sha(source)=='c65b3d80b2207e75ba794821bb41ebc3ffc864959dab27307af77cdeab9d8fd1'
 verification=json.loads((work/'substitution-verification.json').read_text());assert verification['all_passed']and verification['exact_entry_count']==2520
 for name,h in verification['input_sha256'].items():assert sha(ROOT/name)==h,name
 units=json.loads((work/'unit-span.receipt.json').read_text());assert units['status']=='completed'
 with (work/'J0-residual.txt').open()as f:assert list(map(int,f.readline().split()))==[23,57,13]
 attempts={name:json.loads((work/name/'receipt.json').read_text())for name in ['QQa-e-smith','QQe-a-smith','QQe-a-rowgcd','QQe-a-rowgcd-verified']}
 assert all(not r['exact_localized_coker_zero']for r in attempts.values())
 assert attempts['QQe-a-rowgcd']['status']=='arithmetic_error'and not attempts['QQe-a-rowgcd']['actual_nonunit_row_factor_found']
 assert attempts['QQe-a-rowgcd-verified']['status']=='timeout'and not attempts['QQe-a-rowgcd-verified']['arithmetic_error_detected']
 replay=json.loads((work/'QQe-a-rowgcd-verified/checkpoint-replay.json').read_text());assert not replay['exact_12_steps_replayed']
 common=['J0-matrix.txt','row-scalings.txt','prepare.log','units.txt','unit-span.log','unit-span.receipt.json','substitution-verification.json','substitution-verification.log','rowgcd-rational-smoke-matrix.sing']
 members=[old/n for n in ['after21','phase21-check.json','forced21.cpp','forced21','original-row-ids']]+[work/n for n in common]
 for case in ['QQa-e-smith','QQe-a-smith']:
  members += [work/case/n for n in ['matrix.sing','run.sing','receipt.json','solver.log']]
 for case in ['QQe-a-rowgcd','QQe-a-rowgcd-verified','rowgcd-rational-smoke-v3']:
  members += [work/case/n for n in ['run.sing','replay-library.sing','receipt.json','solver.log']]
 members += [work/'QQe-a-rowgcd-verified'/n for n in ['proof.sing','checkpoint-replay.sing','checkpoint-replay.log','checkpoint-replay.json']]
 members += [work/'rowgcd-rational-smoke-v3/proof.sing']
 assert all(f.is_file()and f.resolve().is_relative_to(ROOT)for f in members)
 archive=args.output/'continuation.tar.gz'
 with archive.open('wb')as stream:
  with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:
   with tarfile.open(fileobj=zipped,mode='w')as tar:
    for path in sorted(set(members)):
     info=tar.gettarinfo(str(path),arcname=str(path.relative_to(ROOT)));info.uid=info.gid=0;info.uname=info.gname='';info.mtime=0;info.mode=0o644
     with path.open('rb')as f:tar.addfile(info,f)
 memberdata={str(f.relative_to(ROOT)):{'sha256':sha(f),'bytes':f.stat().st_size}for f in sorted(set(members))}
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_J0_pid_prepare_20261103_v63.cpp',ROOT/'code/s61_nonsquare_J0_pid_verify_20261103_v63.py',ROOT/'code/s61_nonsquare_J0_pid_smith_20261103_v63.py',ROOT/'code/s61_nonsquare_J0_pid_rowgcd_20261103_v63.py',ROOT/'code/polynomial_unit_span_pivots_20261103_v62.cpp',ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json',ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json',ROOT/'research/handoff/20261103_v62/v6_cubic27_coker/manifest.json']
 report={'id':'S61_NONSQUARE_J0_PID_20261103_V63','frozen_continuation':True,'all_exact_preparation_checks_passed':True,
  'locus':'J=5ad-6e-4=0, a*F*B*C*d!=0','J_never_inverted':True,
  'safe_source_shape':[36,70],'safe_old_unit_pivots':21,'d_substitution':'(6e+4)/(5a)',
  'prepared_entry_count_exact_QQ':2520,'prepared_term_count':108530,'prepared_row_scaling_only_a_unit':True,
  'extra_bivariate_known_unit_pivots':13,'bivariate_residual_shape':[23,57],
  'field_unit_results':{'QQ(a)[e]':False,'QQ(e)[a]':False},'actual_defect_curve_proved':False,'J0_cubic_rank_defect_finite_proved':False,
  'traced_Smith_attempts':{name:{k:r[k]for k in ['field','status','CPU_limit_seconds','AS_limit_bytes','seconds','exact_localized_coker_zero']}for name,r in attempts.items()if'smith'in name},
  'invalid_builtin_backend_attempt_excluded_from_all_math_claims':True,
  'row_gcd_valid_attempt':{k:attempts['QQe-a-rowgcd-verified'][k]for k in ['field','status','CPU_limit_seconds','seconds','completed_unit_pivots','arithmetic_error_detected']},
  'fully_completed_forward_steps_conservative':11,'saved_step_witnesses':12,'saved_prefix_independent_replay_confirmed':False,
  'prefix_replay_attempt':replay,'resume_matrix_export_confirmed':False,
  'next_operation':'Replay saved cofactors with per-step progress, export exact24x70 state, then resume row-gcd within a new total<=600s/2.2GB budget; do not repeat first12 gcd searches or unchanged Smith runs.',
  'archive_sha256':sha(archive),'archive_bytes':archive.stat().st_size,'members':memberdata,
  'omitted_easy_to_regenerate_outputs':{'work/nonsquare-J0-pid-v63/J0-residual.txt':{'sha256':sha(work/'J0-residual.txt'),'bytes':(work/'J0-residual.txt').stat().st_size,'replay':'polynomial_unit_span_pivots_20261103_v62.cpp with J0-matrix.txt,units.txt,has-metadata=0,colclean; approximately10seconds'}},
  'input_sha256':{str(f.resolve().relative_to(ROOT)):sha(f)for f in inputs},'seconds':time.monotonic()-started}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_exact_preparation_checks_passed','archive_bytes','J0_cubic_rank_defect_finite_proved','seconds']},indent=2))

if __name__=='__main__':main()
