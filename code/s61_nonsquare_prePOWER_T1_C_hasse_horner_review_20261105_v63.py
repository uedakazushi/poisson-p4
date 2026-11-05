#!/usr/bin/env python3
"""Independent truncated-Horner verification of T1 C-valuation4, no binomial engine."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def plus(a,b):
 assert max(abs(a).bit_length(),abs(b).bit_length())+1<=1580
 c=a+b;assert abs(c).bit_length()<=526;return c
def times(a,b):
 assert abs(a).bit_length()+abs(b).bit_length()<=1580
 c=a*b;assert abs(c).bit_length()<=526;return c
def put(p,k,v):
 z=plus(p.get(k,0),v)
 if z:p[k]=z
 elif k in p:del p[k]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 b=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63';p=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_ENTRY_20261105_v63';rp=p/'report.json'
 assert sha(rp)=='c0272684a8393e1a4fdb800fca785ec8ff86774376d7c63a93b94ef10bd5c746';r=read(rp);assert r['all_passed'] and r['result']['C_valuation']==4
 inputs={rp:sha(rp),b/'native-input.txt':sha(b/'native-input.txt'),Path(__file__):sha(Path(__file__))}
 for name,h in r['input_sha256'].items():assert sha(ROOT/name)==h;inputs[ROOT/name]=h
 for name,h in r['evidence_sha256'].items():assert sha(p/name)==h;inputs[p/name]=h
 bucket=[{}for _ in range(21)];L=0
 with (b/'native-input.txt').open()as f:
  assert f.readline().strip()=='T1CH001 7595 20 3';D=int(f.readline())
  for line in f:
   key,c=map(int,line.split());aa=key//1048576;rest=key-aa*1048576;ee,dd=divmod(rest,1024);k=64*ee+dd
   assert 3<=aa<=23 and k not in bucket[aa-3];bucket[aa-3][k]=c;L+=abs(c)
 assert sum(map(len,bucket))==7595 and L.bit_length()<=363
 # J[n,j]=D*100^n*[y^j] partial Horner polynomial at a0+y.
 # Jnew_j = E*Jold_j +100*Jold_(j-1) +100^n*N_degree*[j=0].
 # n<=21; coefficient l1<=input_l1*100^21*2^21 gives <=526bits.
 assert (((1<<363)-1)*100**21*2**21).bit_length()<=526
 old=[{}for _ in range(5)];scale=1;owned=0;maxbits=0
 for n,degree in enumerate(range(20,-1,-1),1):
  scale=times(scale,100);new=[]
  for j in range(5):
   q={}
   for key,c in old[j].items():
    for delta,w in [(0,2),(64,9),(128,-18)]:put(q,key+delta,times(c,w))
   if j:
    for key,c in old[j-1].items():put(q,key,times(c,100))
   else:
    for key,c in bucket[degree].items():put(q,key,times(c,scale))
   new.append(q);live=16000+2*7595+2*sum(map(len,old))+2*sum(map(len,new))+2*len(q)+2048;owned=max(owned,live);assert live<80000
  old=new;maxbits=max(maxbits,max([abs(c).bit_length()for q in old for c in q.values()]+[0]))
 assert all(not old[j]for j in range(4)) and old[4]
 matched=0
 with (p/'first-nonzero-Hasse.txt').open()as f:
  assert f.readline().strip()=='T1CHV001';k,Dclear=map(int,f.readline().split());assert k==4 and Dclear==D*100**20;n=int(f.readline());native={}
  for _ in range(n):
   ee,dd,c=map(int,f.readline().split());native[64*ee+dd]=c*100
  assert not f.read().strip()
 assert old[4]==native and len(native)==665
 raw=sum(x.stat().st_size for x in inputs)+200000;assert raw<8*1024*1024
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_HORNER_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
 result=dict(id=out.name,all_passed=True,C_valuation_T1_exact=4,all_four_lower_Hasse_layers_independently_zero=True,all665_first_nonzero_coefficients_match=True,independent_method='Five truncated Taylor jets by descending-a Horner, multiplying only E(e)+100*y; no binomial coefficient/power expansion engine imported.',scale_difference='Horner uses D*100^21, exactly100 times native D*100^20.',source_preoperation_committed_bound=526,source_temporary_bound=1580,actual_stored_bits=maxbits,source_owned_peak=owned,hard_owned_limit=80000,RAW_upper=raw,no_quotient_other_factor_tail_point_or_gcd=True,input_sha256={str(x.relative_to(ROOT)):h for x,h in inputs.items()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k}));print(sha(out/'report.json'))
if __name__=='__main__':main()
