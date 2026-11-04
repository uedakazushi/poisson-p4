#!/usr/bin/env python3
"""Protected Qdenominator fibers + all15norm filter onlyon unique-root chart."""
import hashlib,importlib.util,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_EXCEPTION_NORM_FILTER_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 def pin(stem,file='certificate.json'):
  b=ROOT/'reports'/stem;p=b/file;assert sha(p)==load(b/'report.json')['evidence_sha256'][file];return b,p,load(p)
 sb,sp,sc=pin('S61_NONSQUARE_RAW_H_B_EXCEPTION_SQUAREFREE_20261104_v63');S=sc['squarefree_Q'];qb,qp,qc=pin('S61_NONSQUARE_RAW_H_B_COMPONENT_EXACT_QUOTIENT_Q_20261104_v63');exceptions=qc['finite_e_exception_polynomials']+[[[3,'486'],[2,'243'],[1,'-297'],[0,'-54']]];ab,ap,ac=pin('S61_NONSQUARE_RAW_H_B_Q_ROOT_ALL_MINORS_20261104_v63','native-certificate.json');ids=ac['nonzero_remainder_indices'];norms=[[[m[1],v]for m,v in ac['all22_evaluations'][i]['remainder_quadratic_norm_e']]for i in ids];assert len(ids)==15
 assert len(S)==40 and sum(max(m for m,v in p)for p in exceptions)==323
 def line(p):return str(len(p))+''.join(' '+str(m)+' '+str(v)for m,v in p)+'\n'
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  f.write(line(S));f.write('20\n');f.writelines(map(line,exceptions));f.write('15\n');f.writelines(map(line,norms))
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_Q_exception_norm_filter_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-Q-exception-filter'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Qfilter-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');p=lambda P:{(0,m,0):Q(v)for m,v in P if Q(v)};one={(0,0,0):Q(1)};old=p(S);E=p(n['denominator_exception_product']);exc=p(n['protected_exception_gcd']);res=p(n['initial_generic_root_part']);ar.check([(exc,res)],old);ar.check([(exc,p(n['denominator_product_quotient_by_exc']))],E);ar.check([(p(n['protected_gcd_Bezout'][0]),old),(p(n['protected_gcd_Bezout'][1]),E)],exc)
 for source,rec in zip(norms,n['all15_generic_norm_gcd_steps']):
  assert p(rec['before'])==res and p(rec['norm'])==p(source);cg=p(rec['gcd_after']);ar.check([(cg,p(rec['before_quotient']))],res);ar.check([(cg,p(rec['norm_quotient']))],p(source));ar.check([(p(rec['Bezout'][0]),res),(p(rec['Bezout'][1]),p(source))],cg);res=cg
 new=p(n['new_safe_squarefree_cover']);assert res==p(n['final_generic_part']);ar.check([(exc,res)],new);ar.check([(new,p(n['old_cover_quotient_by_new']))],old)
 proof={'safe_union':'Protectexc=gcd(SF39,productE_Q*D); filterONLYSF39/exc byall15norms; newcover=exc*gcd(genericpart,norms).','protected_denominator_fibers_NEVER_discarded':True,'all_source_minor_indices':ids,'all_exact_GCD_Bezout_and_quotient_identities_independent_integer_checked':True,'no_fiber_solve_GB_resultant_point_prime':True,'new_Q_finite_QQ_cover_dimension_upper':8*n['new_degree'],'whole_B_geometric_bound_with_old1600_L36':1600+36+8*n['new_degree']};(out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),sb/'report.json',sp,qb/'report.json',qp,ab/'report.json',ap]
 report={'id':ID,'old_Q_degree':39,'protected_exception_degree':n['protected_degree'],'filtered_generic_part_degree':n['generic_degree'],'new_safe_Q_cover_degree':n['new_degree'],'new_Q_finite_QQ_algebra_dimension_upper':8*n['new_degree'],'new_combined_geometric_B_bound':1600+36+8*n['new_degree'],'native_seconds':n['native_seconds'],'all_independent_integer_checks':True,'fiber_tests_GB_resultant_point_prime_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
