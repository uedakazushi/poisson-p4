#!/usr/bin/env python3
"""One bounded structured L-frame split, heavy5 reduction and six remainder test."""
import argparse,hashlib,json,os,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def treebytes(p):return sum(p.stat().st_size for p in p.rglob('*')if p.is_file())
def main(out):
 start=time.monotonic();deadline=start+60;resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 source=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';pr=source/'report.json';r=json.loads(pr.read_text());assert r['frozen']and r['all_source_bound_coordinate_and_polynomial_input_checks_passed']
 for name in['quad10x26.txt','genuine16x31.txt','certificate.json','native.json']:assert sha(source/name)==r['artifact_sha256'][name]
 r10=ROOT/'reports/S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63/report.json';rr=json.loads(r10.read_text());assert rr['generic_actual_cubic11_rank_atleast10_proved']
 cpp=ROOT/'code/s61_nonsquare_ell_Lframe_borders_20261104_v63.cpp';exe=out/'Lframe-borders'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 stop=None;rss=0
 with(out/'process.log').open('w')as f:
  proc=subprocess.Popen([str(exe),str(source),str(out)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll()is None:
   if time.monotonic()>=deadline:stop='total60s_deadline';break
   if treebytes(out)>32*1024*1024:stop='artifact32MiB_guard';break
   try:
    for line in Path('/proc',str(proc.pid),'status').read_text().splitlines():
     if line.startswith('VmRSS:'):rss=max(rss,int(line.split()[1])*1024)
   except FileNotFoundError:pass
   time.sleep(.05)
  if stop:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
 exe.unlink(missing_ok=True)
 nativep=out/'native.json';n=json.loads(nativep.read_text())if nativep.exists()and proc.returncode==0 else{}
 status=n.get('status','blocked_no_identity_conclusion');complete=status in['all_six_zero_exact','at_least_one_border_nonzero_exact'];stagep=out/'current-stage.json';last=json.loads(stagep.read_text())if stagep.exists()else{}
 proof={'initial26_rows':['q0*L','scaled_q5*L','q0*x0','q0*y1','q0*y2','q0*y3','scaled_q5*x0','scaled_q5*y1','scaled_q5*y2','scaled_q5*y3']+['scaled_G'+str(i)for i in range(16)],
 'initial31_columns':'frozen Lframe preflight genuine31 order: noL20, L*y^2 six, H*(L,x0,y1,y2,y3)',
 'initial_only_known_units':'frozen frame a,C; source G row factors (-40C)^3/C^k',
 'quad_pivot_order':'Cpair two rows first, then A eight rows; no dense19 stage',
 'heavy5_fixed_genuine_relative_rows':[1,2,8,10,11],
 'remaining11_genuine_relative_rows':[0,3,4,5,6,7,9,12,13,14,15],
 'final10_genuine_relative_rows':[0,3,4,5,7,9,12,13,14,15],'last_tested_genuine_relative_row':6,
 'known_generic_lower':'the old accepted nonzero34minor uses all these15 genuine rows and the quad-product image; the invertible L-coordinate/monomial-quotient change therefore makes the first25 rows generically independent',
 'generic_elimination_record':'ELIMINATE t p c followed by polynomials u,v: row_t <- u*row_t-v*row_p. u=pivot/gcd(pivot,targetentry), v=targetentry/gcd. NORMALIZE t followed by g,c: row_t <- row_t/(g*c), all polynomial divisions exact.',
 'generic_factors_only':'all nonzero polynomial row factors/pivots inverted only over QQ(a,e,d); no new global unit claim',
 'checkpoint_scope':'current-matrix/current-stage store only the last FULLY completed pivot checkpoint; a killed operations tail must not be interpreted as a completed step',
 'R11_scope':'completed-R11x16, if present, is the actual canonical11 image after exactly10quad and5heavy pivots on this extra generic chart',
 'last_checkpoint':last,'native':n,'six_remainder_test_completed':complete,'native_returncode':proc.returncode,'wrapper_stop_reason':stop,
 'if_allzero':'generic canonical11 rank10. The final row is a rational combination of the other15 genuine rows and quad products; lift the restricted cubic relation through (Z,W) and multiply ell as in the accepted colon-source11 recipe. Independence from known37 requires its separate review.',
 'if_nonzero':'generic canonical11 rank11, no allpoint conclusion',
 'allpoint57_full239_and_uniform_kernel_count_unchanged':True,'no_newpoint_prime_GB_or_quartic_expansion':True}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 inputs=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_ell_Lframe_preflight_20261104_v63.cpp',ROOT/'research/notes'/(ID+'.md'),pr,source/'certificate.json',source/'quad10x26.txt',source/'genuine16x31.txt',r10,r10.parent/'certificate.json']
 report={'id':ID,'status':status,'six_remainder_test_completed':complete,'generic_rank_exactly10_proved':status=='all_six_zero_exact','generic_rank11_proved':status=='at_least_one_border_nonzero_exact',
 'no_allpoint_or_new_uniform_kernel_claim':True,'no_dense19_rerun_or_newpoint_prime_GB':True,'native_returncode':proc.returncode,'budget_stop_reason':stop,'last_complete_checkpoint':last,'native':n,
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'artifact_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
 'seconds':time.monotonic()-start,'observed_peak_native_RSS_bytes':rss,'combined_artifact_bytes':treebytes(out),'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','artifact_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
