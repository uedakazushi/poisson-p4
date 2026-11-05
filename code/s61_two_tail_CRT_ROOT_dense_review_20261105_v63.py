"""Independent dense tensor convolution of the actual two integer products.

No CRT, gcd, division or producer binary-reader imports. A fresh positive
height/ownership gate precedes every signed product. Existing verified modp
gcd then proves maximality of the independently verified char0 commonfactor.
"""
from pathlib import Path
import hashlib,json,math,resource,signal,struct,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def loadmap(p,magic,which,cap):
 with p.open('rb')as f:
  assert f.read(8)==magic and int.from_bytes(f.read(4),'little')==which
  n=int.from_bytes(f.read(4),'little');assert n<=cap;rows=[];prev=-1
  for _ in range(n):
   key=int.from_bytes(f.read(8),'little');sign=f.read(1);assert sign in [b'\0',b'\1'];l=int.from_bytes(f.read(4),'little');assert l<=1024
   b=f.read(l);assert len(b)==l and (not l or b[-1]);v=int.from_bytes(b,'little');assert v and key>prev;prev=key
   if sign==b'\1':v=-v
   legacy=key%(1<<26);a=legacy>>20;e=(legacy%(1<<20))>>10;d=legacy%1024;assert key//(1<<26)==2*a+e
   rows.append((legacy,v))
  assert not f.read(1)
 return rows
def original(p,magic,n):
 with p.open()as f:
  assert f.readline().strip()==magic;den=int(f.readline());assert den.bit_length()==197;assert int(f.readline())==n
  rows=[tuple(map(int,line.split()))for line in f];assert len(rows)==n
 assert all(v and abs(v).bit_length()<=301 and (i==0 or rows[i-1][0]<k)for i,(k,v)in enumerate(rows))
 return den,rows
def coords(k):return k>>20,(k%(1<<20))>>10,k%1024
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};dense=[];cursor={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def dump(n,x):(out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 base=ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63';pr=pin(base/'report.json');assert sha(pr)=='f2ff75ad00931757c3f0eade5d3250d6d265199f12d86e69e10fd1327e4c32b6';r=json.loads(pr.read_text());assert r['all_passed']
 for n in ['factor.bin','cofactor-N0.bin','cofactor-N2.bin']:assert sha(pin(base/n))==r['evidence_sha256'][n]
 H=loadmap(base/'factor.bin',b'TTFAC001',2,207);Q=[loadmap(base/'cofactor-N0.bin',b'TTCOF001',0,3300),loadmap(base/'cofactor-N2.bin',b'TTCOF001',1,1309)]
 assert len(H)<=207 and sum(map(len,Q))<=2800
 O=[];D=[]
 for st,magic,n in [('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63','T0C4QV001',4500),('S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63','Q1B4QV001',2317)]:
  p=pin(ROOT/'reports'/st/'quotient-integer-common-D.txt');assert sha(p)==r['input_sha256'][str(p.relative_to(ROOT))];d,rows=original(p,magic,n);D.append(d);O.append(rows)
 assert D[0]==D[1]
 assert math.gcd(*(v for k,v in H))==1
 assert max(k for k,v in H)==6*(1<<20)+(1<<10)+4
 assert all(a<=6 and e<=9 and d<=4 and a+e+d<=11 for a,e,d in map(lambda t:coords(t[0]),H))
 assert [max(coords(k)[j]for k,v in H)for j in range(3)]==[6,9,4]
 shapes=[(16,32,18),(12,25,14)];qshapes=[(10,23,14),(6,16,10)]
 for rows,shape in zip(Q,qshapes):assert all(all(x<=b for x,b in zip(coords(k),shape))for k,v in rows)
 normH=sum(abs(v)for k,v in H);normQ=[sum(abs(v)for k,v in q)for q in Q]
 heights=[max(normH.bit_length()+v.bit_length(),max(abs(c).bit_length()for k,c in o))+1 for v,o in zip(normQ,O)]
 assert max(heights)<=1000  # Also bounds a complete dense failure JSON below4MiB.
 # Charge even a full duplicate dense output when saving a failure frontier.
 owned=16000+2*(6817+2800+207)+2*math.prod(s+1 for s in shapes[0]);assert owned==56966 and owned<60000
 raw_read=sum((ROOT/p).stat().st_size for p in pins)+Path(__file__).stat().st_size;raw_prior=raw_read+4*1024*1024+65536;assert raw_prior<8*1024*1024
 dump('pre-arithmetic-dense-tensor-SOURCE.json',dict(input_sha256=pins,positive_height_bounds=heights,positive_normH_bits=normH.bit_length(),positive_normQ_bits=[v.bit_length()for v in normQ],owned_including_failure_copy=owned,cap=60000,RAW_prior=raw_prior,coefficient_cap=8192,bigint_product_cap=600000,hard_seconds=60,soft_seconds=45,dense_shapes=shapes,actual_read_records=[len(H),*[len(q)for q in Q]],no_division_CRT_or_new_gcd=True))
 maximum=0;products=0;checks=[]
 try:
  for which,((A,E,Z),q,originals)in enumerate(zip(shapes,Q,O)):
   dense=[0]*((A+1)*(E+1)*(Z+1));cursor={'which':which,'q':0,'h':0}
   def ix(key):a,e,d=coords(key);return(a*(E+1)+e)*(Z+1)+d
   # Reverse both sparse input orders; use a dense full-box result, unlike
   # producer weighted residual maps/division and sparse multiplyback.
   for qi,(qkey,qv)in enumerate(reversed(q)):
    qi0=ix(qkey)
    for hi,(hk,hv)in enumerate(reversed(H)):
     cursor={'which':which,'q':qi,'h':hi};products+=1;assert products<=600000
     if products%4096==0:assert time.monotonic()-started<45
     assert abs(qv).bit_length()+abs(hv).bit_length()<=heights[which]
     ix0=qi0+ix(hk);v=dense[ix0]+qv*hv;maximum=max(maximum,abs(v).bit_length());assert abs(v).bit_length()<=heights[which];dense[ix0]=v
   pos=0;stream=hashlib.sha256()
   for i,v in enumerate(dense):
    expected=originals[pos][1]if pos<len(originals)and ix(originals[pos][0])==i else 0
    assert v==expected,(which,i)
    if expected:pos+=1
    stream.update(f'{i} {v}\n'.encode())
   assert pos==len(originals);checks.append(dict(which=which,all_dense_coefficients=len(dense),matched_nonzero_original_terms=pos,full_coefficient_stream_sha256=stream.hexdigest()));dense=[]
  # Maximality uses the independently VERIFIED old modp gcd, not a fresh gcd.
  peer=pin(ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63/report.json');assert sha(peer)=='9041a2b792e537070184a5c2b3ef33c9c35b96401a3e9b81199de2286b83652f';assert json.loads(peer.read_text())['all_passed']
  seedpath=pin(ROOT/'reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63/monic-gcd-mod32003.txt')
  with seedpath.open()as f:assert int(f.readline())==3899;seed=dict(tuple(map(int,line.split()))for line in f)
  p=32003;toK=lambda k:sum(c*w for c,w in zip(coords(k),[646,19,1]));hLC=max(H,key=lambda t:toK(t[0]))[1];assert hLC%p
  inv=pow(hLC%p,p-2,p);modH={toK(k):(v*inv)%p for k,v in H if v%p};assert modH==seed and max(modH)==3899
  qdegs=[]
  for q in Q:
   lead=max(q,key=lambda t:toK(t[0]));assert lead[1]%p;qdegs.append(toK(lead[0]))
  assert qdegs==[6645,4076]
  rep=dict(id=ID,all_passed=True,independent_dense_tensor_integer_products=True,original_multivariate_identities=checks,primitive_factor_terms=len(H),factor_degree=[6,9,4],factor_total_degree=11,actual_max_partial_bits=maximum,source_height_bounds=heights,products=products,owned_with_failure_copy=owned,RAW_source_prior=raw_prior,
   exact_char0_gcd_proved_up_to_QQ_scalar=True,maximality_proof='H divides both originalZZinputs by ALLdense coefficients. Hmonicmod32003 is the old independently verifiedgcd3899; bothcofactorKroneckerLCs survivep withdegrees6645/4076. Oldmodularquotients coprime. Any primitivecommonfactor ofchar0cofactors has e<=16,d<=10 soKronleading survives/noalias, wouldreduce nonconstantdivisorofmodgcd1, contradiction.',
   irreducibility_proved=False,Fitting_or_global_rank_or_exhaustion_claimed=False,input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={'pre-arithmetic-dense-tensor-SOURCE.json':sha(out/'pre-arithmetic-dense-tensor-SOURCE.json')},seconds=time.monotonic()-started)
  dump('report.json',rep);print(json.dumps(rep))
 except BaseException as e:
  dump('failure-full-dense-prefix.json',dict(reason=repr(e),cursor=cursor,dense=dense,completed=checks,input_sha256=pins,no_retry=True));raise
if __name__=='__main__':main()
