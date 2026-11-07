#!/usr/bin/env python3
"""ONE R7 absolute test: verify cached4 jets,16 subsets and14 complete products; NO root construction."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];P=32003;CAP=100000;OWN=12000;RAW=2097152;Z=(0,0);ONE=(1,0)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,int):return 1
 if isinstance(x,(list,tuple)):return sum(numeric(v) for v in x)
 if isinstance(x,dict):return len(x)+sum(numeric(v) for v in x.values())
 return 0
def main(custody,custody_sha,out):
 resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 pins={};readbytes=0;peak=0;bits=0
 s=dict(version='R7ABS401',stage='pins',products=0,pending=0,cursor={},field_work={},poly_work={},subsets={},failures=[])
 def pin(p,h):
  nonlocal readbytes
  assert sha(p)==h;readbytes+=p.stat().st_size;pins[str(p.relative_to(ROOT))]=h
 def guard():
  nonlocal peak
  n=6000+numeric(s);peak=max(peak,n);assert n<=OWN,('owned',n);assert time.monotonic()-start<45,'soft_clock'
 def bm(a,b):
  nonlocal bits
  assert 0<=a<P and 0<=b<P;s['pending']=a*b;s['products']+=1;bits=max(bits,s['pending'].bit_length())
  assert s['products']<=CAP and s['pending'].bit_length()<=30
  if s['products']%512==0:guard()
  v=s['pending']%P;s['pending']=0;return v
 def bp(a,n):
  r=1;b=a
  while n:
   s['power_cursor']=dict(a=a,r=r,b=b,n=n,step='multiply')
   if n&1:r=bm(r,b)
   s['power_cursor'].update(r=r,step='square');b=bm(b,b);n>>=1
  s['power_cursor']={};return r
 def ka(x,y):return((x[0]+y[0])%P,(x[1]+y[1])%P)
 def ks(x,y):return((x[0]-y[0])%P,(x[1]-y[1])%P)
 def kn(x):return((-x[0])%P,(-x[1])%P)
 def km(x,y):
  w=dict(x=x,y=y,step='ac');s['field_work']=w;ac=bm(x[0],y[0]);w.update(ac=ac,step='bd');bd=bm(x[1],y[1]);w.update(bd=bd,step='cross')
  xx=(x[0]+x[1])%P;yy=(y[0]+y[1])%P;w.update(sumx=xx,sumy=yy);cross=bm(xx,yy);w.update(cross=cross,step='constant')
  real=(ac-bm(s['q'][0],bd))%P;w.update(real=real,step='linear');imag=(cross-ac-bd-bm(s['q'][1],bd))%P
  s['field_work']={};return(real,imag)
 def smul(a,b,n):
  r=[Z]*n;s['poly_work']=dict(kind='series',left=a,right=b,result=r,cut=n)
  for i in range(n):
   for j in range(n-i):
    s['poly_work'].update(i=i,j=j)
    if a[i]!=Z and b[j]!=Z:r[i+j]=ka(r[i+j],km(a[i],b[j]))
  s['poly_work']={};return r
 def evaluate(root):
  r=[Z]*7;s['eval_work']=dict(root=root,current=r,a=4)
  for a in range(4,-1,-1):
   s['eval_work'].update(a=a,current=r);r=smul(r,root,7)
   for e in range(7):r[e]=ka(r[e],s['Rbar'][a][e])
  s['eval_work']={};return r
 def pmul(a,da,b,db,cut):
  r=[Z]*((da+db+1)*cut);s['poly_work']=dict(kind='bivar',left=a,right=b,result=r,cut=cut,da=da,db=db)
  for i in range(da+1):
   for e in range(7):
    x=a[7*i+e]
    if x==Z:continue
    for j in range(db+1):
     for f in range(min(7,cut-e)):
      y=b[7*j+f];s['poly_work'].update(i=i,e=e,j=j,f=f)
      if y!=Z:k=(i+j)*cut+e+f;r[k]=ka(r[k],km(x,y))
  s['poly_work']={};return r
 def allocate(a,degree,b):
  if not b:return a[:]
  r=[Z]*len(a);s['allocation_work']=dict(input=a,result=r)
  for i in range(degree+1):
   for e in range(7):
    x=a[7*i+e];s['allocation_work'].update(i=i,e=e,x=x,step='55real');v=bm(55,x[0]);s['allocation_work'].update(real=v,step='55imag')
    r[7*i+e]=ka(r[7*i+e],(v,bm(55,x[1])))
    if e<6:
     s['allocation_work']['step']='66real';v=bm(66,x[0]);s['allocation_work'].update(real=v,step='66imag');r[7*i+e+1]=ka(r[7*i+e+1],(v,bm(66,x[1])))
  s['allocation_work']={};return r
 def dump(n,x):
  other=readbytes+sum(p.stat().st_size for p in out.iterdir() if p.is_file());written=0
  with (out/n).open('w') as f:
   for chunk in json.JSONEncoder(separators=(',',':')).iterencode(x):
    written+=len(chunk.encode());assert other+written+65536<=RAW;f.write(chunk)
   f.write('\n')
 def finish(matched):
  s['stage']='complete_fixed_slice_factor_found' if matched else 'complete_R7_absolute_PRIMARY';guard();dump('four-cached-jets-and14-products-certificate.json',s)
  r=dict(id=out.name,all_passed=True,computations_completed=True,R7_geometrically_irreducible_PRIMARY=not matched,
   fixed_modular_R7_factor_only=matched,no_char0_reducibility_claim=matched,failed_candidate_pairs=len(s['failures']),required_if_positive=14,
   no_seed_or_Hensel_root_construction=True,no_new_point_prime_or_factor_family=True,
   products=s['products'],peak_owned=peak,max_field_temporary_bits=bits,input_sha256=pins,
   evidence_sha256={'four-cached-jets-and14-products-certificate.json':sha(out/'four-cached-jets-and14-products-certificate.json')},
   independent_finite_identity_and_theorem_peer_REQUIRED=True,no_Fitting_rank_or_exhaustion_claim=True,seconds=time.monotonic()-start)
  dump('report.json',r);print(json.dumps(r))
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text());assert c['prepared_only'];pin(Path(__file__),c['worker_sha256'])
  for p,h in c['runtime_input_sha256'].items():pin(ROOT/p,h)
  s['prime_trial_remainders']=[[d,P%d] for d in range(2,179)];assert all(v for d,v in s['prime_trial_remainders']) and 178*178<P<179*179
  s['q']=c['quadratic_field'];disc=(bm(s['q'][1],s['q'][1])-bm(4,s['q'][0]))%P;leg=bp(disc,(P-1)//2);assert disc and leg==P-1
  s['field_certificate']=dict(q=s['q'],discriminant=disc,Legendre=leg)
  alpha=140625%P;ia=bp(alpha,P-2);assert bm(alpha,ia)==1;s['alpha_inverse']=ia
  s['Rbar']=[[Z]*7 for _ in range(5)]
  for (a,e),v in c['R7_d1_integer_coefficients']:s['Rbar'][a][e]=(bm(int(v)%P,ia),0)
  assert s['Rbar'][4]==[(55,0),(66,0)]+[Z]*5
  assert [x for row in s['Rbar'] for x in row]==[tuple(x) for x in c['old_Q_e_less7']]
  s['roots']=[[tuple(x) for x in r] for r in c['four_cached_root_jets']];assert len(s['roots'])==4 and all(len(r)==7 for r in s['roots'])
  assert len(set(r[0] for r in s['roots']))==4
  s['stage']='verify_cached_four_R7_root_jets';s['root_residuals']=[evaluate(r) for r in s['roots']]
  assert s['root_residuals']==[[Z]*7 for _ in range(4)]
  s['linears']=[[kn(x) for x in r]+[ONE]+[Z]*6 for r in s['roots']];s['subsets'][0]=[ONE]+[Z]*6;s['stage']='sixteen_subset_products'
  for mask in range(1,16):
   bit=(mask&-mask).bit_length()-1;old=mask^(1<<bit);s['cursor']=dict(mask=mask,old=old,factor=bit);s['subsets'][mask]=pmul(s['subsets'][old],old.bit_count(),s['linears'][bit],1,7)
  assert allocate(s['subsets'][15],4,1)==[x for row in s['Rbar'] for x in row];s['complete_root_product_verified']=True;s['stage']='fourteen_FULL65_products'
  for mask in range(1,15,2):
   other=15^mask;da=mask.bit_count();db=4-da
   for b in (0,1):
    G=allocate(s['subsets'][mask],da,b);Q=allocate(s['subsets'][other],db,1-b);s['candidate_work']=dict(mask=mask,other=other,LC_allocation=b,G=G,Q=Q)
    product=pmul(G,da,Q,db,13);s['candidate_work']['full_product']=product;failure=None
    for a in range(5):
     for e in range(13):
      diff=ks(product[13*a+e],s['Rbar'][a][e] if e<7 else Z)
      if diff!=Z:failure=dict(mask=mask,LC_allocation=b,a=a,e=e,difference=diff);break
     if failure is not None:break
    if failure is None:s['matched_pair']=s['candidate_work'];finish(True);return
    s['failures'].append(failure);s['candidate_work']={};guard()
  assert len(s['failures'])==14;finish(False)
 except BaseException as e:
  s['failure']=repr(e);dump('failure-complete-prefix.json',dict(state=s,input_sha256=pins,no_retry=True))
  dump('report.json',dict(id=out.name,all_passed=False,prefix_only=True,no_irreducibility_claim=True,reason=repr(e),input_sha256=pins));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True);p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())
