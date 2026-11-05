#!/usr/bin/env python3
"""Independent fixed-divisor integer long division and forced scalar normalization."""
from pathlib import Path
from itertools import product
import hashlib,json,math,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
MA=1048576;ME=1024
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def mul(a,b):
 assert abs(a).bit_length()+abs(b).bit_length()<=1556
 z=a*b;assert abs(z).bit_length()<=518;return z
def update(p,k,c):
 assert max(abs(c).bit_length(),abs(p.get(k,0)).bit_length())+1<=1556
 z=p.get(k,0)+c;assert abs(z).bit_length()<=518
 if z:p[k]=z
 elif k in p:del p[k]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 p=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_ENTRY_20261105_v63';rp=p/'report.json';assert sha(rp)=='be5fca9748f1aabdee550ea05b4a3584ddff22b5f43af33b0bb7ebe74ab7736b';r=read(rp);assert r['all_passed']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 for n,h in r['evidence_sha256'].items():assert sha(p/n)==h
 b=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63';inp=b/'native-input.txt';value=p/'quotient-integer-common-D.txt';SCALE=100**20
 # First verify the fixed scalar normalization independently. C is primitive
 # over ZZ, since9-8=1, so Gauss also proves this divisibility abstractly.
 with value.open()as f:
  assert f.readline().strip()=='T1C4QV001';Dout=int(f.readline());nq=int(f.readline());assert nq==4811
  for _ in range(nq):k,c=map(int,f.readline().split());assert c%SCALE==0
  assert not f.read().strip()
 R={};L=0
 with inp.open()as f:
  assert f.readline().strip()=='T1C4Q001 7595 20 3 16';D=int(f.readline());assert Dout==D*SCALE
  for line in f:
   k,c=map(int,line.split());assert k-3*MA not in R;R[k-3*MA]=c;L+=abs(c)
 assert len(R)==7595 and L.bit_length()<=363
 # Independent multinomial enumeration of the fixed four-term C^4 divisor.
 terms=[(MA,100),(2*ME,18),(ME,-9),(0,-2)];div={}
 for xs in product(terms,repeat=4):
  k=sum(x[0]for x in xs);c=math.prod(x[1]for x in xs);update(div,k,c)
 assert div[4*MA]==100**4 and len(div)==25
 S=sum(math.comb(h+3,3)*29**h*100**(16-h)for h in range(17));Qbound=(((1<<363)-1)*S+SCALE-1)//SCALE
 assert (((1<<363)-1)+Qbound*129**4).bit_length()<=518
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63';out.mkdir(exist_ok=False);Q={};steps=0;peak=0;maxbits=0
 while R:
  k=max(R)
  if k//MA<4:break
  c=R[k];assert c%(100**4)==0;qc=c//(100**4);qk=k-4*MA;assert qk not in Q;Q[qk]=qc
  for dk,dc in div.items():update(R,qk+dk,-mul(qc,dc))
  steps+=1;assert steps<=10013;peak=max(peak,18000+2*(len(R)+len(Q)+len(div))+2*len(Q));assert peak<100000
  maxbits=max(maxbits,max([abs(c).bit_length()for c in R.values()]+[0]),abs(qc).bit_length())
  if steps%128==0:assert time.monotonic()-start<50
 assert not R and len(Q)==4811
 with value.open()as f:
  f.readline();f.readline();assert int(f.readline())==len(Q)
  for _ in range(len(Q)):
   k,c=map(int,f.readline().split());assert c//SCALE==Q[k] and c%SCALE==0
  assert not f.read().strip()
 # Publish only the fixed proven scalar normalization, no content search.
 norm=out/'quotient-original-common-D.txt'
 with norm.open('w')as f:
  f.write('T1C4QN001\n'+str(D)+'\n'+str(len(Q))+'\n')
  for k in sorted(Q):f.write(str(k)+' '+str(Q[k])+'\n')
 raw=sum(x.stat().st_size for x in [inp,value,Path(__file__),norm,rp])+200000;assert raw<8*1024*1024
 report=dict(id=out.name,all_passed=True,independent_fixed_C4_long_division_exact=True,full_remainder_zero=True,all4811_quotient_coefficients_match=True,fixed_scalar_100_power20_divides_all_native_coefficients=True,no_generic_polynomial_or_content_GCD=True,original_common_denominator_bits=D.bit_length(),normalized_integer_l1_bits=sum(abs(c)for c in Q.values()).bit_length(),normalized_terms=4811,steps=steps,actual_max_integer_bits=maxbits,source_stored_bound=518,source_temporary_bound=1556,source_owned_peak=peak,hard_owned_limit=100000,RAW_upper=raw,input_sha256={str(x.relative_to(ROOT)):sha(x)for x in [rp,inp,value,Path(__file__)]},evidence_sha256={norm.name:sha(norm)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k}));print(sha(out/'report.json'))
if __name__=='__main__':main()
