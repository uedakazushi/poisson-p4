#!/usr/bin/env python3
"""Independent pre-emission source/budget review of the T1-only C valuation."""
import hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 b=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63';rp=b/'report.json';assert sha(rp)=='355a68985f9c50262933864ab1d32da8fce477e3e2d71ff4fc809197273baf5f';r=read(rp);assert r['all_passed'] and r['value_worker_NOT_run']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 for n,h in r['evidence_sha256'].items():assert sha(b/n)==h
 old=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63/integer-column4-and-result.txt';new=b/'native-input.txt';counts=[0]*21;l1=0;maxbits=0
 with old.open()as f,new.open()as g:
  assert f.readline().strip()=='J1KV001';assert g.readline().strip()=='T1CH001 7595 20 3'
  for expected in [66,94,122,150]:
   node,D=map(int,f.readline().split());assert node==expected;m=int(f.readline())
   for _ in range(m):f.readline()
  node,D=map(int,f.readline().split());assert node==155 and int(g.readline())==D and int(f.readline())==7595
  for _ in range(7595):
   line=f.readline();assert g.readline()==line
   key,c=map(int,line.split());a=key//1048576;rem=key-a*1048576;e,d=divmod(rem,1024)
   assert 3<=a<=23 and 0<=e<=41 and 0<=d<=18 and a+e+d<=61 and c
   counts[a-3]+=1;l1+=abs(c);maxbits=max(maxbits,abs(c).bit_length())
  assert not f.read().strip() and not g.read().strip()
 assert l1.bit_length()<=363 and D.bit_length()<=395
 # H_k(e,d)*D*100^20 = sum N_ied*binom(i,k)*100^(20-i+k)
 #                              *(-18e^2+9e+2)^(i-k)*e^e*d^d.
 # Root-power l1<=29^h, binom(i,k)<=2^20 and29^h<=100^h.
 # This bounds EACH term and EVERY partial sum, including reverse passes.
 worst=((1<<363)-1)*(1<<20)*(100**20);assert worst.bit_length()<=518
 loops=[sum(counts[i]*(2*(i-k)+1)for i in range(k,21))for k in range(21)]
 worst_products=max(sum(loops[:k])+3*loops[k]for k in range(21));assert worst_products<=6528123
 input_cells=7595;root_cells=sum(2*h+1 for h in range(21));layer_cells=82*19;assert root_cells==441 and layer_cells==1558
 owned=2*(input_cells+root_cells+layer_cells)+16000;assert owned==35188
 output_bound=2200000;raw=sum(p.stat().st_size for p in b.iterdir()if p.is_file())+output_bound+100000;assert raw<19442896
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
 report=dict(id=out.name,all_passed=True,source_only=True,no_root_power_or_Hasse_layer_evaluated=True,all7595_input_integer_lines_byte_identical=True,actual_input_l1_bits=l1.bit_length(),actual_input_max_coefficient_bits=maxbits,positive_integer_bound=518,whole_owned_bound=owned,root_power_slots=441,max_output_terms=1558,max_integer_products=worst_products,RAW_upper=raw,formula='D*100^20 H_k=sum N_ied binom(i,k) 100^(20-i+k) (-18e^2+9e+2)^(i-k)e^e d^d; C=100(a-a0).',exact_C_valuation_iff_first_nonzero_Hasse_layer=True,manual_code_review='Only fixed T1; firstnonzero order0..20; before-commit bits/support/ownership, full prefix cursors; reverse-to-zero/reproduce only first nonzero layer; I/O failure cannot emit second full state.',no_quotient_GCD_other_tail_or_Fitting_claim=True,input_sha256={str(rp.relative_to(ROOT)):sha(rp),str(old.relative_to(ROOT)):sha(old),str(new.relative_to(ROOT)):sha(new),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k}));print(sha(out/'report.json'))
if __name__=='__main__':main()
