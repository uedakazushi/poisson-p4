#!/usr/bin/env python3
"""Independent QQ proof and fresh one-pivot replay for the second J0 unit."""
import argparse,gzip,hashlib,json,math,resource,shlex,shutil,subprocess,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 source=ROOT/'work/nonsquare-J0-pid-v63/primitive-span-v1/residual.txt';work=ROOT/'work/nonsquare-J0-pid-v63/primitive-chain-v1/stage-1';previous=ROOT/'reports/S61_NONSQUARE_J0_PRIMITIVE_UNIT_20261103_v63/report.json'
 assert sha(source)==json.loads(previous.read_text())['new_residual_sha256'];cf={}
 with(work/'row13.witness').open()as f:
  assert int(f.readline())==56
  for j in range(56):
   z=f.readline().split()
   if int(z[0]):assert z[0]=='1'and z[2:]==['0','0','0'];cf[j]=F(z[1])
 assert list(cf)==[7,9,26,36,46]
 den=math.lcm(*(v.denominator for v in cf.values()));nums=[v.numerator*(den//v.denominator)for v in cf.values()];content=math.gcd(*nums);primitive=[v//content for v in nums]
 assert primitive==[5757,2400,12320,-22945,8225]
 R,a,e=ring('a,e',sp.QQ);polys={}
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[22,56,1];f.readline();f.readline()
  for i in range(22):
   for j in range(56):
    z=f.readline().split()
    if i==13 and j in cf:
     n=int(z[0]);assert len(z)==1+4*n and all(z[4+4*k]=='0'for k in range(n))
     polys[j]=R.from_dict({tuple(map(int,z[2+4*k:4+4*k])):sp.QQ(z[1+4*k])for k in range(n)})
 delta=sum((n*polys[j]for j,n in zip(cf,primitive)),R.zero);unit=a*a*(100*a+18*e**2-9*e-2)*(3*e+2)**2
 c=delta.exquo(unit);assert c and c.is_ground and delta==c*unit
 lifted=sum((R(sp.QQ(v.numerator,v.denominator))*polys[j]for j,v in cf.items()),R.zero);assert lifted==unit
 payload={'source_shape':[22,56],'target_row':13,'columns':list(cf),'primitive_integer_coefficients':primitive,'nonzero_QQ_scalar':str(c),
  'source_polynomials':{str(j):[[list(ex),str(v)]for ex,v in sorted(p.items())]for j,p in polys.items()},
  'unit_factorization':'a^2*C*(3e+2)^2','full_delta_terms':[[list(ex),str(v)]for ex,v in sorted(delta.items())]}
 with gzip.open(args.output/'identity.json.gz','wt')as f:json.dump(payload,f,separators=(',',':'))
 # Repeat only the cheap saved-witness pivot, never its167second search.
 native=ROOT/'code/polynomial_one_witness_pivot_20261103_v62.cpp';units=ROOT/'work/nonsquare-J0-pid-v63/units.txt'
 flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
 exe=args.output/'replay-pivot';subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-lgmp','-o',str(exe)],check=True)
 out=args.output/'replayed-residual.txt';t=time.monotonic()
 def limits():resource.setrlimit(resource.RLIMIT_AS,(1800000000,1800000000));resource.setrlimit(resource.RLIMIT_CPU,(30,31))
 with(args.output/'fresh-pivot.log').open('w')as f:
  p=subprocess.run([str(exe.resolve()),str(source),str(out.resolve()),str(units),'1',str(work/'row13.witness'),'13'],stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=35)
 assert p.returncode==0 and'PASS 1'in(args.output/'fresh-pivot.log').read_text();assert sha(out)==sha(work/'residual.txt')
 with out.open()as f:assert list(map(int,f.readline().split()))==[21,55,1]
 for name in ['row13.witness','row13.log','pivot.log']:shutil.copyfile(work/name,args.output/name)
 shutil.copyfile(units,args.output/'units.txt')
 # Store only compact evidence; the large output is reproducible in<1second.
 out_bytes=out.stat().st_size;out_sha=sha(out);out.unlink();exe.unlink()
 inputs=[Path(__file__),source,previous,work/'row13.witness',work/'row13.log',native,ROOT/'code/s61_nonsquare_J0_primitive_chain_20261103_v63.py',ROOT/'work/nonsquare-J0-pid-v63/primitive-chain-v1/producer-before-pathfix.py']
 report={'id':'S61_NONSQUARE_J0_SECOND_UNIT_20261103_V63','all_passed':True,'base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]',
  'J_not_inverted':True,'unknown_core_factor_inverted':False,'independent_full_QQ_identity_verified':True,
  'source_shape':[22,56],'target_row':13,'columns':list(cf),'primitive_integer_coefficients':primitive,'nonzero_QQ_scalar':str(c),'unit_factorization':'a^2*C*(3e+2)^2',
  'new_residual_shape':[21,55],'additional_known_unit_pivots':1,'total_safe_pivots_from57basis':36,
  'fresh_saved_witness_pivot_replay_seconds':time.monotonic()-t,'fresh_output_byte_identical_to_export_before_metadata_error':True,
  'controller_failure_was_metadata_only':'Path.relative_to after successful native export; fixed by resolving output path; no search was repeated.',
  'source_search_seconds_before_witness':167.0489888259999,'new_residual_sha256':out_sha,'new_residual_bytes':out_bytes,
  'new_large_residual_omitted_as_fast_reproducible':True,'full_coker_unit_or_rank57_allpoints_proved':False,'J0_defect_support_finite_proved':False,
  'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-started}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','new_residual_shape','unit_factorization','seconds']},indent=2))
if __name__=='__main__':main()
