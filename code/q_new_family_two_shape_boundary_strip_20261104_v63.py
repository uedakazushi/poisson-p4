#!/usr/bin/env python3
"""ONE60s exact named-boundary stripping of the frozen raw resultant."""
import argparse,hashlib,json,math,resource,shlex,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_BOUNDARY_STRIP_20261104_v63';BYTECAP=32*1024*1024
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def cap(seconds):
 resource.setrlimit(resource.RLIMIT_CPU,(seconds,seconds));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 started=time.monotonic();deadline=started+60;resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 rawdir=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63';rawpath=rawdir/'resultant.txt';rawreport=rawdir/'report.json';raw=read(rawreport);assert raw['all_passed']and raw['native_result']['FLINT_resultant_return_code']==1 and not raw['native_result']['zero']
 assert sha(rawpath)==raw['evidence_sha256']['resultant.txt']
 native=ROOT/'code/q_new_family_two_shape_boundary_strip_20261104_v63.cpp';binary=work/'named-strip-native';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint','openssl'],text=True));error=None;nr=None
 try:
  subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(.1,deadline-time.monotonic()-2))
  remaining=deadline-time.monotonic()-2;assert remaining>0
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(rawpath),str(out),str(BYTECAP-1_000_000)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=remaining,preexec_fn=lambda:cap(max(1,math.ceil(remaining))))
  nr=read(out/'native-report.json');assert nr['all_passed']and nr['resultant_calls']==0 and nr['general_factorization_calls']==0
  assert nr['input']['canonical_SHA256']==raw['native_result']['canonical_term_SHA256']
  assert nr['global_factorization_identity_verified']and nr['final_integer_content_one']
  assert nr['visible_r_power']==22 and nr['visible_Z0_power']==1
  for i,check in enumerate(nr['final_maximality_remainder_identities']):
   assert check['q_times_factor_plus_remainder_equals_final_quotient']and check['remainder_nonzero']
  files={'quotient':nr['final_quotient']}
  for label,check in zip(['rminus1','T'],nr['final_maximality_remainder_identities']):
   files['nondiv_quotient_'+label]=check['division_quotient'];files['nondiv_remainder_'+label]=check['division_remainder']
  for name,stat in files.items():
   h=hashlib.sha256();n=0
   with(out/(name+'.txt')).open('rb')as f:
    assert f.readline().split()==[b'6',str(stat['terms']).encode()]
    for line in f:h.update(line);n+=1
   assert n==stat['terms']and h.hexdigest()==stat['canonical_SHA256']
 except Exception as e:error=repr(e)
 model=dict(identity='raw=integer_content*r^22*Z0*(r-1)^a*T^b*quotient',variable_order=['s','r','R','Z0','Z1','Z2'],T='(R-2)*r-R',
  a=nr['rminus1_power']if nr else None,b=nr['T_power']if nr else None,integer_content=nr['integer_content']if nr else None,
  all_successful_intermediate_quotients='reconstruct by multiplying final quotient by the remaining suffix of the stored factor chain; exact step hashes and identities in native-report',
  chart_unit_reason='r,r-1,T,Z0 already units in the accepted Mobius/two-shape source chart; the nonzero integer is a QQ unit',
  retained_other_boundaries=['J','Fn','Hn','Nq','s','xi'],
  no_generic_degree_or_primitive_r_polynomial_claim=True,raw_resultant_identification_still_pending_independent_identity_audit=True)
 (out/'factor-chain.json').write_text(json.dumps(model,indent=2)+'\n')
 report=dict(id=STEM,all_passed=error is None,error=error,status='named_boundary_factorization_and_maximality_verified'if error is None else 'bounded_failure',
  one_shared_deadline_seconds=60,address_space_bytes=2_000_000_000,resultant_calls=0,general_factorization_calls=0,
  no_newpoints_primes_source26_or_r_N_elimination=True,raw_packet_unchanged=True,
  factor_powers=dict(r=22,Z0=1,rminus1=nr['rminus1_power']if nr else None,T=nr['T_power']if nr else None),
  final_quotient=nr['final_quotient']if nr else None,library_resultant_identity_not_independently_rechecked=True,
  generic_fiber_degree_actual_target_index_NEW_equation_unresolved=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),native,rawpath,rawreport]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<BYTECAP-10000
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','error','factor_powers','final_quotient','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
