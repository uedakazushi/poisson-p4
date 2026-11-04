#!/usr/bin/env python3
"""ONE fixedsafeQ37 exact QQ factorization; finitecomponent inventory only."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_EXCEPTION_FACTOR_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True);helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('arith',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_Q_EXCEPTION_NORM_FILTER_20261104_v63';ip=base/'native-certificate.json';assert sha(ip)==load(base/'report.json')['evidence_sha256'][ip.name];S=load(ip)['new_safe_squarefree_cover'];P={(0,m,0):Q(v)for m,v in S};den=math.lcm(*(x.denominator for x in P.values()));integers=[int(x*den)for x in P.values()];cc=math.gcd(*integers);bits=max(abs(x//cc).bit_length()for x in integers);assert len(P)==38 and bits==198
 inp=out/'native-input.txt';inp.write_text(ar.line(P));cpp=ROOT/'code/s61_nonsquare_raw_H_B_common_curve_factor_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-Q-exception-factor'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Q37factor-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');e=s.symbols('e');cs,fs=s.factor_list(s.Poly.from_dict({(m[1],):s.Rational(v.numerator,v.denominator)for m,v in P.items()},(e,),domain=s.QQ));known={};cst=Q(str(cs))
 for p,ex in fs:
  cst*=Q(str(p.LC()))**ex;known[tuple(sorted(((0,m[0],0),str(v))for m,v in p.monic().terms()))]=ex
 actual={tuple(sorted((tuple(m),str(Q(v)))for m,v in r['factor'])):r['exponent']for r in n['factors']};assert known==actual and cst==Q(n['factor_constant'])
 prod={(0,0,0):Q(n['factor_constant'])};profiles=[]
 def multiply(a,b):
  aa,ad=ar.intp(a);bb,bd=ar.intp(b);r={}
  for m,c in aa.items():
   for z,d in bb.items():
    k=tuple(x+y for x,y in zip(m,z));r[k]=r.get(k,0)+c*d
  return{k:Q(v,ad*bd)for k,v in r.items()if v}
 for i,r in enumerate(n['factors']):
  p=ar.poly(r['factor']);assert r['exponent']==1
  prod=multiply(prod,p);de=max(m[1]for m in p);profiles.append({'factor_index':i,'degree_e':de,'terms':len(p),'dimension_QQ_before_a_clopen_upper':8*de,'full_ideal_matrix_entry_upper':22*(8*de)**2,'eligible_degree_at_most8':de<=8,'factor':r['factor']})
 assert prod==P;choices=sorted([r for r in profiles if r['eligible_degree_at_most8']],key=lambda r:(r['degree_e'],r['terms'],r['factor_index']));assert choices
 cert={'source_safe_SF37':S,'factor_constant':n['factor_constant'],'irreducible_factors':n['factors'],'all_irreducible_factors_independent_SymPy_and_expansion_verified':True,'factor_inventory':profiles,'selected_smallest_component_factor_index':choices[0]['factor_index'],'selection_rule':'minimumdegree,terms,nativefactorindex; ONLYonecomponent authorizednext','all_other_components_retained_not_tested':True,'no_membership_fiber_point_prime_GB_resultant_job_in_factor_stage':True};(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip]
 report={'id':ID,'factor_degrees':[r['degree_e']for r in profiles],'selected_smallest_factor_index':choices[0]['factor_index'],'selected_degree':choices[0]['degree_e'],'selected_factor':choices[0]['factor'],'native_seconds':n['native_seconds'],'all_independent_factorization_and_reconstruction_checked':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
