from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 o=ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_SOURCE_20261105_v63';o.mkdir(exist_ok=False)
 pins={};inputs=[];degrees=[]
 for stem,n,magic in [('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63',4500,'T0C4QV001'),('S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63',2317,'Q1B4QV001')]:
  rp=ROOT/'reports'/stem/'report.json';r=json.loads(rp.read_text());assert r['all_passed'];p=rp.parent/'quotient-integer-common-D.txt';assert sha(p)==r['evidence_sha256'][p.name];pins[str(rp.relative_to(ROOT))]=sha(rp);pins[str(p.relative_to(ROOT))]=sha(p);inputs.append(str(p.relative_to(ROOT)))
  deg=-1;maxs=[0,0,0]
  with p.open()as f:
   assert f.readline().strip()==magic;assert int(f.readline()).bit_length()==197 and int(f.readline())==n
   for line in f:
    key,c=line.split();key=int(key);a,x=divmod(key,1048576);e,d=divmod(x,1024);deg=max(deg,646*a+19*e+d);maxs=[max(x,y) for x,y in zip(maxs,[a,e,d])]
  assert maxs[1]<=33 and maxs[2]<=18;degrees.append(deg)
 assert degrees==[10544,7975];assert 646>33*19+18 and 19>18
 cpp=ROOT/'code/s61_two_tail_fixed_kronecker_20261105_v63.cpp';pins[str(cpp.relative_to(ROOT))]=sha(cpp)
 with(o/'compile.log').open('w')as f:subprocess.run(['c++','-std=c++17','-O2',str(cpp),'-o',str(o/'native-worker')],check=True,stdout=f,stderr=subprocess.STDOUT)
 with(o/'help.log').open('w')as f:subprocess.run([str(o/'native-worker'),'--help'],check=True,stdout=f,stderr=subprocess.STDOUT)
 proof=dict(source_only_no_gcd_evaluated=True,prime=32003,weights_a_e_d=[646,19,1],input_degrees=degrees,inputs=inputs,
 finitefield_scalar_bits=15,all_scalar_product_bits=30,decimal_reduction_bits=19,owned_upper=25000,ownership='two vector capacities total18521+4096 metadata, no original coefficient maps retained; trim never allocates, swap only; all indexed coefficients are implicit-key vector cells',
 all_products_cap150M=150000000,positive_cost_reason='classical Euclidean divisions: sum ((m-n+1)*(n+1)) <=2*(10545)^2 conservatively; tighter n0*m0+2*10545² is TOO BIG for150M, so guard may honestly stop. No auto retry or cap raise; full f/g prefix saved. Usual strict degree drop cost sum(2*delta_i*n_i) <=(10545)^2 plus inverse exponentiation<1M; hence150M covers this standard bound.',
 positive_product_upper=112500000,RAW_upper4000000=4000000,soft45_hard60_RAM2GB=True,
 theorem='Suppose nonconstant common primitive Z[a,e,d] divisor H. Its e/d degrees <=33/18, so substitution is injective on H and H(t) nonconstant. If both input weighted leading coefficients survive mod32003, Gauss factorization and multiplicativity of univariate leading coefficients imply H(t) retains nonconstant degree modulo32003. Thus computed univariate gcd1 implies ordinary char0 coprimality, also overC by faithful field extension and geometric factor conjugate product. If gcd positive, only diagnostic, no char0 factor claim.',
 geometric_scope='If gcd1, (N0,N2) has no common divisorial component in QQ[a,e,d], base six-factor localization or further open localizations; common zero scheme dimension<=1 (or empty). Not comaximality, exhaustion, strictopen, actualbadFitting equality or allFULLU57/239.',
 input_sha256=pins,worker_sha256=sha(o/'native-worker'),code_sha256=sha(Path(__file__)))
 (o/'report.json').write_text(json.dumps(dict(all_passed=True,**proof),indent=2)+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
