#!/usr/bin/env python3
"""One authorized120s/6GB final norm audit from the frozen primitive PRS prefix."""
import argparse,gzip,hashlib,json,math,resource,shlex,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_FINAL_NORM_IDENTITY_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main(out,work):
 started=time.monotonic();deadline=started+119;resource.setrlimit(resource.RLIMIT_AS,(6_000_000_000,6_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 prefix=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_PRS_CONTENT_PREFLIGHT_20261104_v63';pr=read(prefix/'report.json');assert pr['all_passed']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in pr['evidence_sha256'].items():assert sha(prefix/rel)==h,rel
 for name in ['C2','C1']:
  with gzip.open(prefix/(name+'-primitive.txt.gz'),'rb')as f,(work/(name+'.txt')).open('wb')as g:
   for b in iter(lambda:f.read(1<<20),b''):g.write(b)
 raw=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63/resultant.txt'
 cpp=ROOT/'code/q_new_family_two_shape_final_norm_review_20261104_v63.cpp';binary=work/'independent-final-norm';error=None;ret=None
 try:
  flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint','openssl'],text=True))
  subprocess.run(['c++','-std=c++17','-O2',str(cpp),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(.1,deadline-time.monotonic()-2))
  remaining=deadline-time.monotonic()-2;assert remaining>0
  def cap():
   resource.setrlimit(resource.RLIMIT_AS,(6_000_000_000,6_000_000_000));n=max(1,math.ceil(remaining));resource.setrlimit(resource.RLIMIT_CPU,(n,n))
  with(out/'native.log').open('w')as log:ret=subprocess.run([str(binary),str(work/'C2.txt'),str(work/'C1.txt'),str(raw),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=remaining,preexec_fn=cap).returncode
  if ret:raise RuntimeError('independent final-norm native exit '+str(ret))
 except Exception as e:error=repr(e)
 nr=read(out/'native-report.json')if(out/'native-report.json').exists()else None
 for p in list(out.glob('*.txt')):
  with p.open('rb')as f,gzip.open(p.with_suffix('.txt.gz'),'wb')as g:
   for b in iter(lambda:f.read(1<<20),b''):g.write(b)
  p.unlink()
 inputs=[Path(__file__),cpp,ROOT/'code/q_new_family_two_shape_prs_content_preflight_20261104_v63.cpp',ROOT/'code/q_new_family_two_shape_prs_review_20261104_v63.cpp',ROOT/'research/notes'/(STEM+'.md'),
         prefix/'report.json',prefix/'C2-primitive.txt.gz',prefix/'C1-primitive.txt.gz',raw,ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_SUBRESULTANT_PRS_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_COPRIME_DEGREES_20261104_v63_REVIEW/report.json']
 report={'id':STEM,'all_passed':error is None and bool(nr and nr['all_passed']),'error':error,'native_exit_code':ret,'native_result':nr,
         'one_authorized_DIFFERENT_final_norm_job_total_wall_seconds_limit':120,'current_execution_wall_seconds_limit':119,'prior_compile_only_failure_seconds':0.430552747,'prior_compile_only_failure_no_native_or_math_job_started':True,'address_space_bytes':6_000_000_000,'artifact_bytes_limit':32*1024*1024,
         'failed_original_two_GB_PRS_receipt_preserved':True,'no_full_PRS_resultant_gcd_point_prime_source26_or_automatic_retry':True,
         'raw_r_projection_identity_only_no_NEW_equation_or_target_field_index_one':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-started,'frozen':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','error','native_result','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.work.resolve())
