#!/usr/bin/env python3
"""ONE native22-minor GCD, own Fraction reconstructions, independent reverse SymPy GCD."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MINOR_GCD_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(x+y for x,y in zip(m,n));v=z.get(k,Q())+c*d
   if v:z[k]=v
   else:z.pop(k,None)
 return z
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';ip=base/'certificate.json';pr=read(base/'report.json');assert sha(ip)==pr['evidence_sha256'][ip.name];c=read(ip);A=c['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];U=c['FULLU_product'];assert len(A)==22 and sum(map(len,A))==2963
 native=out/'native-input.txt';cpp=ROOT/'code/s61_nonsquare_raw_H_B_minor_gcd_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-minor-gcd'
 def line(p):return str(len(p))+''.join(' '+str(v)+' '+' '.join(map(str,m))for m,v in p)+'\n'
 with native.open('w')as f:
  f.write('22\n')
  for p in A+[U]:f.write(line(p))
 command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as f:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','B-gcd-native',str(exe),str(native),str(out/'native-certificate.json')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 nr=read(out/'native-certificate.json');g=poly(nr['gcd22']);h=poly(nr['remaining_h_coprime_U']);up=poly(U);assert g and h and len(nr['all22_literal_quotients'])==22
 for source,q in zip(A,nr['all22_literal_quotients']):assert mul(g,poly(q))==poly(source)
 current=g
 for step in nr['localizer_removal_steps']:
  t=poly(step['removed_gcd_with_U']);after=poly(step['remaining_h_after']);assert mul(t,after)==current and mul(t,poly(step['U_quotient_by_removed']))==up;current=after
 assert current==h and poly(nr['final_gcd_h_U'])=={(0,0,0):Q(1)}
 # A different exact implementation/order certifies greatest common divisor and final coprimality.
 R,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:R.from_dict({m:s.QQ(str(v))for m,v in poly(p).items()});monic=lambda p:p/p.LC
 reverse=monic(rp(A[-1]));checkchain=[]
 for p in reversed(A[:-1]):reverse=monic(reverse.gcd(rp(p)));checkchain.append(len(reverse))
 nativeg=rp(nr['gcd22']);assert reverse==nativeg
 independently_coprime=monic(rp(nr['remaining_h_coprime_U']).gcd(rp(U)));assert independently_coprime==R.one
 nonconstant=any(any(m)for m in h);classification='genuine_B_rankbad_divisor_on_FULLU'if nonconstant else'no_B_rankbad_hypersurface_on_FULLU_only_codimension_at_least2'
 proof={'coefficient_domain':'QQ[a,e,d] localized at U=adFBCJ','classification':classification,'native_exact_gcd22':nr['gcd22'],'remaining_h':nr['remaining_h_coprime_U'],
  'all22_common_factor_quotients_and_reconstructions_verified_by_own_Fraction':True,'every_removed_factor_supported_in_U_verified_by_literal_U_quotients':True,
  'independent_reverse_order_SymPy_exact_multivariate_gcd_matches_FLINT':True,'reverse_chain_support_counts':checkchain,'independent_fixed_gcd_h_U_is1':True,
  'no_irreducibility_assumption_on_named_factors_no_factorization':True,'logic_if_nonconstant':'h divides every actual localized maximal minor, is nonconstant and coprimeU, so at least one height1 component V(h) survives the openD(U) and is genuine B-rankbad. No irreducible factor needs to be selected.',
  'logic_if_constant':'After removing allU-supported factors, no height1 prime of the localized3variableUFD contains all maximalminors. Thus support has codimension>=2/dimension<=1, or is empty. It does NOT establish the unitideal, emptiness, or globalB3.',
  'no_GB_retry_point_prime_primarydec_or_full239_claim':True}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip]
 report={'id':ID,'classification':classification,'native_gcd_terms':len(g),'native_gcd_degree':max(sum(m)for m in g),'localizer_supported_removal_steps':len(nr['localizer_removal_steps']),'remaining_h_terms':len(h),'remaining_h_degree':max(sum(m)for m in h),
  'all22_quotients_and_independent_GCD_coprimality_verified':True,'native_seconds':nr['native_seconds'],'GB_jobs':0,'newpoint_prime_or_factorization':False,'globalB3_NOT_inferred_when_h_constant':True,
  'native_executable_sha256_not_committed':exesha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
