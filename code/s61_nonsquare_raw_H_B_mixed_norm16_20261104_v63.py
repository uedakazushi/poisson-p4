#!/usr/bin/env python3
"""ONE mixed-m1 norm16 determinant/olda117 gcd, no full solve."""
import argparse,hashlib,json,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pd(p):return {int(k):int(v)for k,v in p}
def mul(a,b):
 z={}
 for i,x in a.items():
  for j,y in b.items():z[i+j]=z.get(i+j,0)+x*y
 return{k:v for k,v in z.items()if v}
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 ib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_INPUT_20261104_v63';ob=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63_REVIEW';assert sha(ib/'report.json')=='6a66674a009094d22012221f7ea881df38d1d848976a2c33184b8cf851193a34';assert read(ib/'report.json')['all_passed']
 c=read(ib/'certificate.json');cost=read(ib/'cost-preflight.json');D=int(cost['matrix_uniform_QQ_denominator']);matrix=[]
 for row in c['matrix_rows_source_order']:
  rr=[]
  for p in row:
   z={}
   for m,s in p:
    assert m[1:]==[0,0];q=Q(s)*D;assert q.denominator==1;z[m[0]]=q.numerator
   rr.append(z)
  matrix.append(rr)
 old={}
 for m,s in read(ob/'independent-certificate.json')['nonzero_primitive_a117']:
  assert m[1:]==[0,0];q=Q(s);assert q.denominator==1;old[m[0]]=q.numerator
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for p in[p for row in matrix for p in row]+[old]:
   f.write(str(len(p))+'\n')
   for k,v in sorted(p.items(),reverse=True):f.write(str(k)+' '+str(v)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_mixed_norm16_20261104_v63.cpp';exe=work/'mixed-norm16'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','mixed-norm16',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 n=read(out/'native-certificate.json')
 if n['determinant_zero']:classification='zero_m1_norm_STOP_no_next_mixed_coefficient';checks={'zero_determinant_stop':True}
 else:
  det=pd(n['determinant']);primitive=pd(n['primitive']);res=pd(n['new_primitive_without_a']);g=pd(n['gcd']);qn=pd(n['quotient_new']);qo=pd(n['quotient_old']);content=int(n['determinant_QQ_content']);ap=n['removed_known_a_power']
  assert det=={i:content*v for i,v in primitive.items()}and primitive=={i+ap:v for i,v in res.items()};assert mul(g,qn)==res and mul(g,qo)==old and pd(n['old_primitive_a117'])==old
  if n['Bezout_saved']:
   S,T=pd(n['Bezout_S']),pd(n['Bezout_T']);lhs=mul(S,res)
   for k,v in mul(T,old).items():lhs[k]=lhs.get(k,0)+v
   assert{k:v for k,v in lhs.items()if v}=={i:int(n['Bezout_constant'])*v for i,v in g.items()}
  classification='old_bad_on_DG_EMPTY_G0_RETAINED'if n['gcd_is_QQunit']else'nonconstant_gcd_necessary_a_cover_retained'
  checks=dict(whole_content_a_power_and_gcd_product_identities=True,Bezout_whole_integer_identity_checked=bool(n['Bezout_saved']),gcd_degree=n['gcd_degree'],classification=classification,G0_Q_L_cases_NOT_removed=True,nonzero_norm_degrees={'a':max(det),'min_a':min(det),'without_a':max(res)},no_actual_bad_existence_claim=True)
 (out/'independent-product-checks.json').write_text(json.dumps(checks,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',ib/'report.json',ib/'certificate.json',ib/'cost-preflight.json',ob/'report.json',ob/'independent-certificate.json']
 report=dict(id=ID,all_passed=True,classification=classification,determinant_zero=n['determinant_zero'],no_full16_solve_or_new_point_prime_GB=True,scope_only_old_DG_G0_retained=True,native_seconds=n.get('native_seconds'),input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},executable_SHA_not_published=sha(exe),seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(json.dumps(checks,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as error:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(error),'partial_checkpoints_preserved':[p.name for p in a.output.glob('*checkpoint.json')],'no_auto_retry_or_m2':True},indent=2)+'\n')
  raise
