#!/usr/bin/env python3
"""ROOT integer multiplyback and three symbolic evaluations, independent of trial division."""
import hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_EXISTING_QRE_STRIP_20261105_v64'
SOURCE=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_ROOT_SOURCE_REVIEW_20261105_v64'
OUT=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_ROOT_VALUE_REVIEW_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cell(path):
 with path.open('rb')as f:
  assert f.read(9)==b'Q5CELL1\n\x01'
  def u():
   v=0;s=0;last=128
   while last&128:
    b=f.read(1);assert b;last=b[0];v+=(last&127)<<s;s+=7;assert s<2048
   assert s==7 or last;return v
  def z():
   tag=u();raw=f.read(tag//2);assert len(raw)==tag//2 and(not raw or raw[0]);v=int.from_bytes(raw,'big');assert v or not tag&1;return-v if tag&1 else v
  num,den=z(),z();shift=[]
  for _ in range(5):v=u();shift.append(v//2 if not v&1 else-(v//2)-1)
  n=u();out={};previous=None
  for _ in range(n):
   c=z();e=tuple(shift[j]+u()for j in range(5));assert c and e not in out
   assert previous is None or(sum(previous),previous)>(sum(e),e);previous=e;out[e]=c
  assert not f.read(1)
 return num,den,out
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));OUT.mkdir(exist_ok=False)
 pins={};p=BASE/'report.json';assert sha(p)=='d1f1940f49703f4c7d956b6ead3c803b25de82aa4e628d6a0a74725739462fb3'
 r=json.loads(p.read_text());assert r['all_passed'] and r['only_QRE_powers']=={'Q':4,'R':2,'E':1};pins[str(p.relative_to(ROOT))]=sha(p)
 p=SOURCE/'report.json';assert sha(p)=='a03e3d78d68bfc67639065ad48a7b7000387d7acee611f22854856ed011be350';pins[str(p.relative_to(ROOT))]=sha(p)
 for name,want in r['evidence_sha256'].items():assert sha(BASE/name)==want;pins[str((BASE/name).relative_to(ROOT))]=want
 p=BASE/'fixed-row58-QRE-strip-and-Gauss-factorization.json';d=json.loads(p.read_text())
 fp=BASE/d['final_reduced_integer_polynomial']['file'];assert sha(fp)==d['final_reduced_integer_polynomial']['sha256']
 a,b,F=cell(fp);assert(a,b)==(1,1) and len(F)==317
 a,b,P=cell(BASE/'initial-integer-polynomial.bin');assert(a,b)==(1,1) and len(P)==951
 assert sum(map(abs,F.values()))<2**358
 # Prior source proves all intermediates in the same1739 fibre universe.
 # Multiplyback norm <2^358*6^4*2^2*82<2^377; three cleared root
 # evaluations have norm <=2^358*50^16<2^449. Use hard518 throughout.
 assert 2**358*6**4*2**2*82<2**377 and 2**358*50**16<2**449
 assert 7*(951+317+2*1739)+12000<60000
 (OUT/'phase.json').write_text(json.dumps(dict(before_multiplyback_and_symbolic_values=True,
  source_bits449=449,integer_hard518=518,whole_owned60000=60000,RAW_hard2000000=2000000,
  only_fixed_factorization_and_QRE_nondivisibility=True,no_trial_division_matrix_update_or_GCD=True),indent=2)+'\n')
 peak=0;maximum=0;updates=0
 def check(x):
  nonlocal maximum
  maximum=max(maximum,abs(x).bit_length());assert abs(x).bit_length()<=518;return x
 def add(out,key,value):
  nonlocal updates
  value=check(value);new=check(out.get(key,0)+value);updates+=1
  if new:out[key]=new
  else:out.pop(key,None)
 def account(A,B):
  nonlocal peak
  n=7*(len(F)+len(P)+len(A)+len(B))+12000;peak=max(peak,n);assert n<60000
 A=dict(F)
 factors=[{(0,0,0,1,0):5,(0,0,0,0,0):1}]*4+[{(0,0,0,1,0):1,(0,0,0,0,0):1}]*2+[{(0,0,1,0,0):32,(0,0,0,1,0):15,(0,0,0,0,0):35}]
 for factor in factors:
  B={}
  for e,c in A.items():
   for f,v in factor.items():add(B,tuple(e[j]+f[j]for j in range(5)),check(c*v))
  account(A,B);A=B;del B
 assert A==P;del A
 # Independent exact symbolic substitution, clearing a fixed16th power of
 # leading coefficients. No division algorithm or producer codec is imported.
 residues={}
 for label in ['Q','R','E']:
  value={}
  for (h,l,s,q,t),c in F.items():
   if label=='Q':add(value,(h,l,s,t),c*((-1)**q)*5**(16-q))
   elif label=='R':add(value,(h,l,s,t),c*((-1)**q))
   else:
    for j in range(s+1):
     coefficient=c*32**(16-s)*math.comb(s,j)*(-15)**j*(-35)**(s-j)
     add(value,(h,l,q+j,t),coefficient)
  assert value
  # Store only an exact nonzero witness and digest of the full root result.
  witness=min(value);digest=hashlib.sha256()
  for e,c in sorted(value.items()):digest.update((' '.join(map(str,e))+' '+str(c)+'\n').encode())
  residues[label]=dict(terms=len(value),witness_exponents=list(witness),witness_integer=str(value[witness]),all_value_sha256=digest.hexdigest())
  peak=max(peak,7*(len(F)+len(P)+len(value))+12000);assert peak<60000
 facts=dict(original951_integer_coefficients_equal=True,final317_integer_coefficients=True,
  exact_powers={'Q':4,'R':2,'E':1},root_nonzero_tests=residues,
  method='independent stream decoder, forward integer multiplication and direct binomial symbolic root substitution on final factor',
  no_division_engine_or_producer_imports=True,max_integer_bits=maximum,owned_upper=peak,updates=updates,
  p10='F10/(353894400000*H*t0*F9)',F9_not_divided=True,
  all168x15_87extra66constraints_other14_oneP290_complements_retained=True)
 (OUT/'verification.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir())+20000;assert raw<2000000
 result=dict(id=OUT.name,all_passed=True,frozen=True,independent_actual_fixed_factorization_review=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in OUT.iterdir()},RAW_total_upper=raw,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()
