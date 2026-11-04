#!/usr/bin/env python3
"""Independent symbolic same-prime irreducibility reduction certificate."""
from pathlib import Path
from fractions import Fraction as Q
from functools import reduce
from math import gcd,lcm,comb
import hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1];STEM='S61_NONSQUARE_PREPOWER_R434_MODP_IRREDUCIBILITY_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 work=ROOT/'work/prePOWER-R434-symbolic-modp-factor-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,n='report.json'):return load('reports/'+st+'/'+n)
 r=get('S61_NONSQUARE_PREPOWER_R434_FACTOR_20261104_v63');assert r['all_passed'] and hashes['reports/S61_NONSQUARE_PREPOWER_R434_FACTOR_20261104_v63/report.json']=='d7b711aec63e1c3335b476078d8e44ae33ce38b80f7dc12ced97143d686de01d'
 for p,h in r['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for n,h in r['evidence_sha256'].items():p='reports/S61_NONSQUARE_PREPOWER_R434_FACTOR_20261104_v63/'+n;assert sha(p)==h;hashes[p]=h
 chart=get('S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW');assert chart['all_passed']
 assert hashes['reports/S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW/report.json']=='3a140e24f3cc2d2c97a9a3632ed15e92f87912f9ad5bbc19e7ac32e7f5c0e7dc'
 strip=get('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json');pol=strip['remaining_polynomial'];assert len(pol)==434
 coef=[Q(c) for m,c in pol];den=lcm(*(q.denominator for q in coef));ints=[int(q*den) for q in coef];g=reduce(gcd,map(abs,ints));content=Q(g,den)
 assert content==Q(131072,853553295) and content.numerator%32003 and content.denominator%32003
 primitive=[int(q/content) for q in coef];assert all(Q(q,1)*content==c for q,c in zip(primitive,coef))
 assert reduce(gcd,map(abs,primitive))==1 and max(abs(q).bit_length() for q in primitive)==98
 degree=max(sum(m) for m,c in pol if Q(c));assert degree==20
 reduced={tuple(m):q%32003 for (m,_),q in zip(pol,primitive) if q%32003};assert max(map(sum,reduced))==20
 dense=comb(23,3);visible=2*434+20*dense+4*dense+15000;assert visible<150000
 pre={'same_symbolic_prime':32003,'variables':['a','e','d'],'no_parameter_point_or_Saux_divisibility':True,'source434_terms':434,'QQ_content':str(content),'primitive_Z_content':1,'primitive_coefficient_bits':98,'QQ_and_modular_total_degree':20,'dense_degree20_envelope':dense,'visible_coefficients_upper':visible,'factor_internal_support_not_instrumentable':True,'external_limits':{'seconds':60,'memory_bytes':2<<30,'artifact_bytes':32<<20,'visible_coefficients':150000},'char0_scope_if_irreducible_modp':'QQ[a,e,d] only; no absolute/C irreducibility inference.'}
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 text=[str(content),str(len(pol))]+[' '.join(map(str,m))+' '+str(c) for m,c in pol]+[str(len(pol))]+[' '.join(map(str,m))+' '+str(q) for (m,_),q in zip(pol,primitive)]
 ip=work/'input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<32<<20
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-started)))
 with (work/'native.log').open('w') as log:subprocess.run([str(exe),str(ip),str(out/'native.json'),str(out/'stage.log')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,57-(time.monotonic()-started)))
 n=read('reports/'+STEM+'/native.json');assert n['all_factorization_checks_passed'] and n['factor_API_status']==1 and n['full_factor_product_matches_primitive']
 conclusion=bool(n['irreducible_mod32003'])
 proof={'Gauss_primitive_Z_model':'R434=content*Rprimitive, gcd of primitive integer coefficients1; content numerator/denominator32003-units.',
  'total_degree_preserved':20,'specialization_argument':'If primitive R factors into two nonconstant integer polynomials, reduction is a product of nonzero polynomials. Full total degree20 after reduction forces both factor degrees to be preserved, so an irreducible F32003 reduction contradicts such a factorization.',
  'QQ_irreducible_inference_enabled':conclusion,'irreducible_modp_is_not_absolute_irreducibility':True,'source_bound_actual_R434':True,'no_new_numeric_source_point_or_prime':True,'no_Saux_division_test':True}
 (out/'Gauss-degree-preservation-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in [str(Path(__file__).relative_to(ROOT)),str(cpp.relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':conclusion,'complete_symbolic_modular_factorization_and_product_checks_passed':True,'irreducible_mod32003':conclusion,'QQ_R434_irreducible_independently_proved':conclusion,'absolute_or_complex_irreducibility_claimed':False,'source_total_degree_preserved':20,'no_parameter_point_new_prime_Saux_division_or_QQ_retry':True,'factor_API_internal_support_not_instrumented':True,'input_sha256':hashes,'ignored_reproducible_input':{'path':str(ip.relative_to(ROOT)),'sha256':sha(str(ip.relative_to(ROOT)))},'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
