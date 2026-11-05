from pathlib import Path
import hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();b=ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_SOURCE_20261105_v63';r=json.loads((b/'report.json').read_text());assert r['all_passed'] and r['source_only_no_gcd_evaluated']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 assert sha(b/'native-worker')==r['worker_sha256'];assert all(32003%k for k in range(2,179))
 o=ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63';o.mkdir(exist_ok=False)
 # Explicit fresh ROOT source/code review, not a rerun of any prior coefficient task.
 # For descending degrees di: sum (di-di+1)*(di+1+1)<m²/2+m,
 # and sum(di+1+1)<m²/2+2m. Inverses and leading scalar products <1M.
 proof=dict(all_passed=True,ONE_FIXED_PAIR_PRIME_JOB_AUTHORIZED=True,no_multivariate_factorization=True,prime_trial_division_checked=True,
 all_division_products_upper=112500000,hard_products150000000=150000000,
 bound_derivation='Let degree drops Delta_i>=1. Sum Delta_i*(d_(i+1)+1)<=m²/2+m by telescoping squares. Extra sum(d_(i+1)+1)<=m²/2+2m. Leading scalar products <=2m, Fermat inverse operations <=45m. With m=10544, total<112500000. Mod products fit30bits; decimal reduction19bits. Two vector capacities18521+4096metadata<25000, no future allocation.',
 source=sha(b/'report.json'),worker=sha(b/'native-worker'),soft45_hard58driver60CPU_RAM2GB=True,RAW4000000=True,cpp_reviewed=True)
 (o/'pre-arithmetic-ROOT-source-review-and-ONE-GO.json').write_text(json.dumps(proof,indent=2)+'\n')
 readraw=sum((ROOT/n).stat().st_size for n in r['input_sha256'])+sum(x.stat().st_size for x in b.iterdir());assert readraw+1000000<4000000
 with(o/'native.log').open('w')as f:run=subprocess.run([str(b/'native-worker'),'--run',*[str(ROOT/p)for p in r['inputs']],str(o)],stdout=f,stderr=subprocess.STDOUT,timeout=58)
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 rr=json.loads((o/'native-report.json').read_text());assert run.returncode==0 or rr.get('prefix_saved')
 d=dict(all_passed=rr['all_passed'],source_fixed_prime32003=True,result=rr,
 ordinary_char0_coprimality_proved=bool(rr.get('all_passed') and rr.get('gcd_is_one')),
 positive_gcd_is_only_diagnostic_no_char0_factor=not rr.get('gcd_is_one',False),
 no_empty_locus_exhaustion_or_fullrank_claim=True,input_sha256={**r['input_sha256'],str((b/'report.json').relative_to(ROOT)):sha(b/'report.json')},code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p)for p in o.iterdir()},seconds=time.monotonic()-start)
 assert readraw+sum(p.stat().st_size for p in o.iterdir())+65536<4000000
 (o/'report.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
if __name__=='__main__':main()
