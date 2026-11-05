from pathlib import Path
import hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();o=ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63';o.mkdir(exist_ok=False);src=ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63';r=json.loads((src/'report.json').read_text());assert r['all_passed'] and r['result']['gcd_degree']==3899
 paths=[ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63/quotient-integer-common-D.txt',ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63/quotient-integer-common-D.txt',src/'monic-gcd-mod32003.txt'];cpp=ROOT/'code/s61_two_tail_gcd_quotient_verify_v2_20261105_v63.cpp';pins={str(p.relative_to(ROOT)):sha(p)for p in paths+[cpp,src/'report.json']}
 assert pins[str(paths[2].relative_to(ROOT))]==r['evidence_sha256'][paths[2].name]
 with(o/'compile.log').open('w')as f:subprocess.run(['c++','-std=c++17','-O2',str(cpp),'-o',str(o/'native-worker')],check=True,stdout=f,stderr=subprocess.STDOUT)
 with(o/'help.log').open('w')as f:subprocess.run([str(o/'native-worker'),'--help'],check=True,stdout=f,stderr=subprocess.STDOUT)
 proof=dict(fresh_independent_certificate_scope='GivenH divides bothoriginalKroneckerinputs by sparsemoniclongdivision; quotientscoprime implies givenH exactgcd. DoesnotrecomputeproducerrawEuclid.',prime32003=True,no_new_prime_or_char0factor=True,sparseH94=True,owned_upper30000=30000,whole_lifetime='FirstR10545+Q6646+Hcapacity<=256+4096; releaseR. SecondR7976+priorQa6646+Qb4077+H<=256+4096 <30000. ThenonlyquotientvectorsandH remain. Everyallocationchecked; nooriginalvectorcopy.',positive_product_upper=70000000,proof_products='sparse divisions <=94*(6646+4077)=1007962; quotientEuclid <=6646^2+1M <46M, so70M positivecovers all; hard150M unchanged.',modular_scalar_bits30=True,RAW_upper4000000=True,soft45_hard58driver60CPU_RAM2GB=True,input_sha256=pins,code_sha256=sha(Path(__file__)),binary_sha256=sha(o/'native-worker'))
 (o/'pre-arithmetic-source.json').write_text(json.dumps(proof,indent=2)+'\n')
 with(o/'native.log').open('w')as f:run=subprocess.run([str(o/'native-worker'),'--run',*[str(p)for p in paths],str(o)],stdout=f,stderr=subprocess.STDOUT,timeout=58)
 n=json.loads((o/'native-report.json').read_text());assert run.returncode==0 or n.get('prefix_saved')
 for p,h in pins.items():assert sha(ROOT/p)==h
 result=dict(all_passed=n['all_passed'],independent_exactdivision_and_quotient_coprimality=True,no_char0_commonfactor_or_coprime_inference=True,result=n,input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p)for p in o.iterdir()},seconds=time.monotonic()-start)
 assert sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in o.iterdir())+65536<4000000
 (o/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
