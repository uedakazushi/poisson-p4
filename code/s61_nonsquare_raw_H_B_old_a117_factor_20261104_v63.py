#!/usr/bin/env python3
"""ONE exactnative QQ factorinventory ofcompleted a117; productcheck only nofibers."""
import gzip,hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_OLD_A117_FACTOR_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True);helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_CANONICAL_20261104_v63';np=base/'native-certificate.json.gz';br=load(base/'report.json');assert sha(np)==br['evidence_sha256'][np.name];c=json.loads(gzip.decompress(np.read_bytes()));P=ar.poly(c['remaining_non_a_annihilator']);assert len(P)==118 and max(m[0]for m in P)==117;assert all(m[1]==m[2]==0 for m in P);inp=out/'native-input.txt';inp.write_text(ar.line(P));cpp=ROOT/'code/s61_nonsquare_raw_H_B_common_curve_factor_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-old-a117-factor'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','old117factor-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');prod={(0,0,0):Q(n['factor_constant'])};inventory=[]
 def multiply(a,b):
  aa,ad=ar.intp(a);bb,bd=ar.intp(b);r={}
  for m,c in aa.items():
   for z,d in bb.items():
    k=tuple(x+y for x,y in zip(m,z));r[k]=r.get(k,0)+c*d
  return{k:Q(v,ad*bd)for k,v in r.items()if v}
 for i,r in enumerate(n['factors']):
  p=ar.poly(r['factor']);assert all(m[1]==m[2]==0 for m in p)
  for j in range(r['exponent']):prod=multiply(prod,p)
  stats=ar.stats(p);inventory.append({'factor_index':i,'degree_a':stats['degrees'][0],'exponent':r['exponent'],'terms':len(p),'max_rational_numerator_bits':stats['max_coefficient_numerator_bits'],'factor':r['factor']})
 assert prod==P and sum(r['degree_a']*r['exponent']for r in inventory)==117
 cert={'source_a117':[[list(m),str(v)]for m,v in sorted(P.items(),reverse=True)],'factor_constant':n['factor_constant'],'factor_inventory':inventory,'native_exact_factor_algorithm_completed':True,'independent_integer_whole_product_reconstruction_verified':True,'irreducibility_source':'nativeFLINTQQunivariatefactoralgorithm; no secondlargefactorsearchrun','all_factors_preserved_no_unknown_unit_removed':True,'no_finitefiber_query_GB_point_prime_job':True};(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),base/'report.json',np]
 report={'id':ID,'factor_degrees_exponents':[[r['degree_a'],r['exponent']]for r in inventory],'sum_degrees_with_multiplicity':117,'native_seconds':n['native_seconds'],'independent_product_reconstruction':True,'no_finitefiber_query':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
