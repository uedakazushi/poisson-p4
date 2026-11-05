#!/usr/bin/env python3
"""Independent pre-execution proof of the fixed C^4 quotient operator and bounds."""
import hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 b=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63';rp=b/'report.json';assert sha(rp)=='2697636cf24fd85244f7f017f7f2cd3d940be9717c0a1748c1805f99417b51a0';r=read(rp);assert r['all_passed'] and r['value_worker_NOT_run']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(b/p)==h
 old=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63/native-input.txt';new=b/'native-input.txt';counts=[0]*21;L=0
 with old.open()as f,new.open()as g:
  assert f.readline().strip()=='T1CH001 7595 20 3';assert g.readline().strip()=='T1C4Q001 7595 20 3 16';D=int(f.readline());assert int(g.readline())==D
  for _ in range(7595):
   s=f.readline();assert s==g.readline();key,c=map(int,s.split());counts[key//1048576-3]+=1;L+=abs(c)
  assert not f.read().strip() and not g.read().strip()
 assert L.bit_length()<=363
 # Fixed expansion (100a-R)^(-4) at a=infinity, followed by coefficient
 # extraction. Divisibility is already proven by four exact zero Hasse layers.
 S=sum(math.comb(h+3,3)*29**h*100**(16-h)for h in range(17))
 Ql1=((1<<363)-1)*S;Cl1=129**4;rem_l1=((1<<363)-1)*100**20+Ql1*Cl1
 assert max(Ql1.bit_length(),rem_l1.bit_length())<=518
 qcap=sum(1 for a in range(17)for e in range(34)for d in range(19)if a+e+d<=50)
 rcap=sum(1 for a in range(21)for e in range(42)for d in range(19)if a+e+d<=58)
 assert qcap==10013 and rcap==14991
 pcap=sum(2*h+1 for h in range(17));ccap=sum((h+1)**2 for h in range(5));slab=74*19
 quotient_owned=16000+2*(7595+pcap+qcap+slab+1)
 back_owned=16000+2*(7595+ccap+qcap+rcap+1)
 assert max(quotient_owned,back_owned)<=81776
 inner=sum(counts[j]*(2*(j-i-4)+1)for i in range(17)for j in range(i+4,21))
 scaling=20+3*16**2+2*sum(counts[j]for i in range(17)for j in range(i+4,21))+4*sum(h*h for h in range(1,5))+7595+1
 allp=inner+qcap*25+scaling;assert allp<=2712113
 # Full simultaneous states fit without a second output/prefix. Decimal
 # coefficient width<=156, packedkey<=10digits, separators/line<=4bytes.
 prefix=(max(quotient_owned,back_owned)-16000)//2*170+50000
 raw=sum(p.stat().st_size for p in b.iterdir()if p.is_file())+prefix+200000;assert raw<10610960
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
 report=dict(id=out.name,all_passed=True,source_only=True,no_quotient_slab_rootpower_or_multiplyback_evaluated=True,all7595_input_lines_identical=True,known_C_valuation4_independent_peer='ae012d21cbef98ac9554efb634412cd0062014e630a1aed7572582386ccdd406',operator='D*100^20 Q_i=sum_(j>=i+4) N_j binom(j-i-1,3) R^(j-i-4)100^(20-j+i)',quotient_total_l1_bound_bits=Ql1.bit_length(),full_remainder_l1_bound_bits=rem_l1.bit_length(),hard_bit_bound=518,quotient_support=qcap,partial_slab_support=slab,remainder_support=rcap,source_owned_peak=max(quotient_owned,back_owned),contract_owned_limit=81776,all_product_upper=allp,RAW_upper=raw,manual_code_review='Retain whole high-e partial slabs; validate final degree bounds after complete; release all R powers before C cache; verify every coefficient of C^4Q-100^20 U1, including missing negative-a consistency. Before-commit guards/full cursors, no I/O second state.',no_other_tail_factor_GCD_or_bad_locus_claim=True,input_sha256={str(rp.relative_to(ROOT)):sha(rp),str(old.relative_to(ROOT)):sha(old),str(new.relative_to(ROOT)):sha(new),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},seconds=time.monotonic()-t)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k}));print(sha(out/'report.json'))
if __name__=='__main__':main()
