#!/usr/bin/env python3
"""ONE exact factor of fixedG49; independent QQ factorization/reconstruction."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar);base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63';np=base/'native-certificate.json';assert sha(np)==load(base/'report.json')['evidence_sha256'][np.name];G=ar.poly(load(np)['final_common_curve_gcd']);assert len(G)==49
 den=math.lcm(*(x.denominator for x in G.values()));contents=math.gcd(*(int(x*den)for x in G.values()));bits=max(abs(int(x*den)//contents).bit_length()for x in G.values());pre={'source_terms':49,'source_degree_a_e':[6,12],'source_box_upper':91,'primitive_integer_coefficient_bits':bits,'preflight_small_support':True};(out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');inp=out/'native-input.txt';inp.write_text(ar.line(G));cpp=ROOT/'code/s61_nonsquare_raw_H_B_common_curve_factor_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-common-curve-factor'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','curvefactor-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');const=Q(n['factor_constant']);factors=[(ar.poly(x['factor']),x['exponent'])for x in n['factors']];a,e=s.symbols('a e');GG=s.Poly.from_dict({(m[0],m[1]):s.Rational(v.numerator,v.denominator)for m,v in G.items()},(a,e),domain=s.QQ);cs,fs=s.factor_list(GG);sm={};con=Q(str(cs))
 for p,m in fs:
  lc=Q(str(p.LC()));con*=lc**m;pm=p.monic();key=tuple(sorted(((ex[0],ex[1],0),str(v))for ex,v in pm.terms()));sm[key]=m
 nm={tuple(sorted((m,str(v))for m,v in p.items())):ex for p,ex in factors};assert const==con and nm==sm
 # Independent denominator-cleared exact reconstruction of allfactors.
 prod={(0,0,0):const}
 for p,ex in factors:
  for j in range(ex):
   pi,di=ar.intp(prod);qi,ei=ar.intp(p);r={}
   for m,c in pi.items():
    for nn,d in qi.items():
     z=tuple(x+y for x,y in zip(m,nn));r[z]=r.get(z,0)+c*d
   prod={m:Q(v,di*ei)for m,v in r.items()if v}
 assert prod==G
 cert={'source_G':[[list(m),str(v)]for m,v in sorted(G.items(),reverse=True)],'factor_constant':str(const),'factors':n['factors'],'all_native_expansion_and_independent_integer_reconstruction_checked':True,'independent_SymPy_QQ_factor_list_matches_all_monic_factors_and_exponents':True,'factor_profiles':[{'exponent':ex,'terms':len(p),'degree_a_e_d':ar.stats(p)['degrees'],'totaldegree':ar.stats(p)['totaldegree']}for p,ex in factors],'no_pseudodivision_resultant_GB_point_prime_jobs':True,'factor_curves_only_necessary_not_actualbadexistence':True};(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),base/'report.json',np]
 report={'id':ID,'factors':cert['factor_profiles'],'native_seconds':n['native_seconds'],'independent_factorization_and_reconstruction':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps({k:v for k,v in report.items()},indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
