#!/usr/bin/env python3
"""One LIGHT60s known-content preflight on the saved, unchanged PRS prefix."""
import argparse,gzip,hashlib,json,math,resource,shlex,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_PRS_CONTENT_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main(out,work):
 started=time.monotonic();deadline=started+60;resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 old=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_SUBRESULTANT_PRS_20261104_v63_REVIEW';r=read(old/'report.json');assert not r['all_passed']and r['last_frontier']['stage']=='C1_saved'
 for name in ['C2','C1']:
  p=old/(name+'.txt.gz');assert sha(p)==r['evidence_sha256'][p.name]
  with gzip.open(p,'rb')as f,(work/(name+'.txt')).open('wb')as g:
   for b in iter(lambda:f.read(1<<20),b''):g.write(b)
 native=ROOT/'code/q_new_family_two_shape_prs_content_preflight_20261104_v63.cpp';bin=work/'known-content-only';error=None;ret=None
 try:
  flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint','openssl'],text=True))
  subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(bin)],check=True,timeout=max(.1,deadline-time.monotonic()-2))
  remaining=deadline-time.monotonic()-2;assert remaining>0
  def cap():
   resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));n=max(1,math.ceil(remaining));resource.setrlimit(resource.RLIMIT_CPU,(n,n))
  with(out/'native.log').open('w')as f:ret=subprocess.run([str(bin),str(work/'C2.txt'),str(work/'C1.txt'),str(out)],stdout=f,stderr=subprocess.STDOUT,timeout=remaining,preexec_fn=cap).returncode
  if ret:raise RuntimeError('known-content native exit'+str(ret))
 except Exception as e:error=repr(e)
 stats=read(out/'native-report.json')if(out/'native-report.json').exists()and error is None else None
 for p in list(out.glob('*.txt')):
  with p.open('rb')as f,gzip.open(p.with_suffix('.txt.gz'),'wb')as g:
   for b in iter(lambda:f.read(1<<20),b''):g.write(b)
  p.unlink()
 inputs=[Path(__file__),native,ROOT/'code/q_new_family_two_shape_prs_review_20261104_v63.cpp',ROOT/'research/notes'/(STEM+'.md'),old/'report.json',old/'C1.txt.gz',old/'C2.txt.gz']
 report={'id':STEM,'all_passed':error is None and bool(stats and stats['all_passed']),'error':error,'statistics':stats,'one_LIGHT_job_wall_seconds_limit':60,'address_space_bytes':2_000_000_000,'artifact_bytes_limit':32*1024*1024,
         'only_saved_prefix_common_known_factor_strip_no_norm_OR_PRS_rerun':True,'full_resultant_identity_remains_unverified':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-started,'frozen':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','error','statistics','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.work.resolve())
