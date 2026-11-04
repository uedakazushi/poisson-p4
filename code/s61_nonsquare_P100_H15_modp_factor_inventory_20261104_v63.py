#!/usr/bin/env python3
"""ONE factor inventory of the frozen monic H15 at the EXISTING prime32003."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';r=read(base/'report.json');n=read(base/'native.json')
 assert sha(base/'report.json')=='b44096ebdca0b36570e269aa6e171dd065f4c8b8dbff64e777f5db6b6afd890d'
 assert sha(base/'native.json')=='161888f3ad31a27f72e492ec7454c336c99f27817730588feea59693166580ed'
 assert r['complete'] and n['gcd_degree']==15 and n['all_passed']
 H=n['exact_gcd_monic'];assert max(m[0]for m,c in H)==15
 coeff=[0]*16;p=32003;den=[]
 for m,c in H:
  q=Q(c);assert q.denominator%p;den.append(dict(exponent=m[0],denominator=str(q.denominator),residue=q.denominator%p));coeff[m[0]]=q.numerator*pow(q.denominator,-1,p)%p
 assert coeff[15]==1
 ip=out/'native-input.txt';ip.write_text('16\n'+' '.join(map(str,coeff))+'\n')
 (out/'phase.json').write_text(json.dumps(dict(phase='input_bound_before_single_factor_job',source_H_monic=H,all_QQ_denominators_prime_units=den,prime=p,no_new_beta_or_zeta_point=True),indent=2)+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'factor-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,55-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,57-(time.monotonic()-start)))
 f=read(out/'native.json');assert f['all_passed'] and f['full_product_equals_Hbar'] and f['Hbar']==coeff
 degs=[(z['degree'],z['multiplicity'])for z in f['factors']];assert sum(d*m for d,m in degs)==15
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',base/'report.json',base/'native.json',base/'certificate.json']
 report=dict(id=ID,all_passed=True,frozen=True,independent_review_pending=True,prime=p,factor_degrees_and_multiplicities=degs,
   max_irreducible_degree=max(d for d,m in degs),full_product_and_Rabin_certificates_saved=True,
   no_new_prime_beta_or_zeta_point_or_coefficient_rank_job=True,no_char0_factor_or_full239_claim=True,
   input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
 try:main(a.output.resolve())
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=repr(e),no_auto_retry=True),indent=2)+'\n')
  raise
