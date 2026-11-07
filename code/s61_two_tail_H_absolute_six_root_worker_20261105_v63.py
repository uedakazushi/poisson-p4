#!/usr/bin/env python3
"""ONE fixed absolute-H test via six root jets and62 complete recombinations; help is inert."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];P=32003;CAP=2000000;OWN=24000;RAW=4194304
Z=(0,0);ONE=(1,0)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,int):return 1
 if isinstance(x,(list,tuple)):return sum(numeric(y) for y in x)
 if isinstance(x,dict):return len(x)+sum(numeric(y) for y in x.values())
 return 0
def main(custody,custody_sha,out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60)
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={};readbytes=0;peak=0;bits=0
 s=dict(version='HABS6001',stage='pins',p=P,products=0,pending=0,cursor={},field_work={},poly_work={},eval_work={},roots=[],subsets={},updates=[],failures=[])
 def pin(p,h):
  nonlocal readbytes
  assert sha(p)==h;readbytes+=p.stat().st_size;pins[str(p.relative_to(ROOT))]=h
 def guard():
  nonlocal peak
  n=6000+numeric(s);peak=max(peak,n);assert n<=OWN,('owned',n)
  assert time.monotonic()-start<45,'soft_clock'
 def bm(a,b):
  nonlocal bits
  assert 0<=a<P and 0<=b<P
  s['pending']=a*b;s['products']+=1;bits=max(bits,s['pending'].bit_length())
  assert s['products']<=CAP and s['pending'].bit_length()<=30
  if s['products']%2048==0:guard()
  x=s['pending']%P;s['pending']=0;return x
 def bp(a,n):
  r=1;b=a
  while n:
   s['base_power_cursor']=dict(a=a,r=r,b=b,n=n,step='multiply')
   if n&1:r=bm(r,b)
   s['base_power_cursor'].update(r=r,step='square');b=bm(b,b);n>>=1
  s['base_power_cursor']={};return r
 def bi(a):
  assert a;r=bp(a,P-2);assert bm(a,r)==1;return r
 def ka(x,y):return((x[0]+y[0])%P,(x[1]+y[1])%P)
 def ks(x,y):return((x[0]-y[0])%P,(x[1]-y[1])%P)
 def kn(x):return((-x[0])%P,(-x[1])%P)
 def km(x,y):
  w=dict(x=x,y=y,step='ac');s['field_work']=w
  ac=bm(x[0],y[0]);w.update(ac=ac,step='bd');bd=bm(x[1],y[1]);w.update(bd=bd,step='cross')
  # BOTH sums reduce before multiplication, preserving the30-bit premise.
  xx=(x[0]+x[1])%P;yy=(y[0]+y[1])%P;w.update(sumx=xx,sumy=yy)
  cross=bm(xx,yy);w.update(cross=cross,step='quadratic_constant');real=(ac-bm(s['q'][0],bd))%P
  w.update(real=real,step='quadratic_linear');imag=(cross-ac-bd-bm(s['q'][1],bd))%P
  s['field_work']={};return(real,imag)
 def ki(x):
  s['inverse_work']=dict(value=x,step='conjugate')
  conjugate=((x[0]-bm(s['q'][1],x[1]))%P,(-x[1])%P);s['inverse_work'].update(conjugate=conjugate,step='norm');norm=km(x,conjugate)
  assert norm[1]==0 and norm[0];s['inverse_work'].update(norm=norm,step='base_inverse');v=bi(norm[0])
  s['inverse_work'].update(norm_inverse=v,step='real_scale');r0=bm(conjugate[0],v);s['inverse_work'].update(real=r0,step='imag_scale')
  r=(r0,bm(conjugate[1],v));s['inverse_work'].update(result=r,step='verify');assert km(x,r)==ONE;s['inverse_work']={}
  return r,dict(value=x,conjugate=conjugate,norm=norm[0],norm_inverse=v,inverse=r)
 def kval(c,x):
  r=Z
  for j in range(len(c)-1,-1,-1):r=ka(km(r,x),(c[j],0))
  return r
 def bmul(a,b):
  r=[0]*(len(a)+len(b)-1);s['base_poly_work']=dict(kind='product',left=a,right=b,result=r)
  for i,x in enumerate(a):
   for j,y in enumerate(b):s['base_poly_work'].update(i=i,j=j);r[i+j]=(r[i+j]+bm(x,y))%P
  s['base_poly_work']={};return r
 def bdiv(a,b):
  assert b[-1]==1;r=a[:];q=[0]*(len(a)-len(b)+1);s['base_poly_work']=dict(kind='monic_division',left=a,right=b,current=r,quotient=q)
  for k in range(len(a)-len(b),-1,-1):
   t=r[k+len(b)-1];q[k]=t
   for j,v in enumerate(b):s['base_poly_work'].update(k=k,j=j,t=t);r[k+j]=(r[k+j]-bm(t,v))%P
  assert all(x==0 for x in r[len(b)-1:]);s['base_poly_work']={};return q,r[:len(b)-1]
 def smul(a,b,n):
  r=[Z]*n;s['poly_work']=dict(kind='series',left=a,right=b,result=r,cut=n)
  for i in range(n):
   for j in range(n-i):
    s['poly_work'].update(i=i,j=j)
    if a[i]!=Z and b[j]!=Z:r[i+j]=ka(r[i+j],km(a[i],b[j]))
  s['poly_work']={};return r
 def evaluate(root,n):
  r=[Z]*n;s['eval_work']=dict(root=root,cut=n,current=r,a=6)
  for a in range(6,-1,-1):
   s['eval_work'].update(a=a,current=r);r=smul(r,root,n)
   for e in range(n):r[e]=ka(r[e],s['Fbar'][a][e])
  s['eval_work']={};return r
 def pmul(a,da,b,db,cut):
  r=[Z]*((da+db+1)*cut);s['poly_work']=dict(kind='bivar',left=a,left_degree=da,right=b,right_degree=db,result=r,cut=cut)
  for i in range(da+1):
   for e in range(10):
    x=a[10*i+e]
    if x==Z:continue
    for j in range(db+1):
     for f in range(min(10,cut-e)):
      y=b[10*j+f];s['poly_work'].update(i=i,e=e,j=j,f=f)
      if y!=Z:
       k=(i+j)*cut+e+f;r[k]=ka(r[k],km(x,y))
  s['poly_work']={};return r
 def allocate(a,degree,b):
  if b==0:return a[:]
  r=[Z]*len(a);s['allocation_work']=dict(input=a,result=r,degree=degree)
  for i in range(degree+1):
   for e in range(10):
    s['allocation_work'].update(i=i,e=e);x=a[10*i+e]
    s['allocation_work'].update(step='55real',value=x);real=bm(55,x[0]);s['allocation_work'].update(step='55imag',scaled_real=real)
    r[10*i+e]=ka(r[10*i+e],(real,bm(55,x[1])))
    if e<9:
     s['allocation_work'].update(step='66real');real=bm(66,x[0]);s['allocation_work'].update(step='66imag',scaled_real=real)
     r[10*i+e+1]=ka(r[10*i+e+1],(real,bm(66,x[1])))
  s['allocation_work']={};return r
 def dump(name,obj):
  p=out/name;other=readbytes+sum(x.stat().st_size for x in out.iterdir() if x.is_file())
  written=0
  with p.open('w') as f:
   for chunk in json.JSONEncoder(separators=(',',':')).iterencode(obj):
    written+=len(chunk.encode());assert other+written+65536<=RAW,('RAW',other+written)
    f.write(chunk)
   f.write('\n')
 def finish(matched):
  s['stage']='completed_fixed_image_reducible' if matched else 'completed_absolute_irreducible_PRIMARY';guard()
  dump('complete-six-root-certificate.json',s)
  r=dict(id=out.name,all_passed=True,computations_completed=True,
   absolute_irreducible_mod32003_PRIMARY=not matched,geometric_H_irreducible_PRIMARY=not matched,
   fixed_image_factor_found_only=matched,no_char0_reducibility_claim=matched,
   failed_complete_candidate_pairs=len(s['failures']),required_if_positive=62,
   one_fixed_d1_prime32003_only=True,no_other_sample_prime_factor_family=True,
   products=s['products'],peak_owned=peak,max_field_temporary_bits=bits,input_sha256=pins,
   evidence_sha256={'complete-six-root-certificate.json':sha(out/'complete-six-root-certificate.json')},
   independent_identity_and_theorem_peer_REQUIRED=True,no_Fitting_fullrank_or_exhaustion_claim=True,seconds=time.monotonic()-start)
  dump('report.json',r);print(json.dumps(r))
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text());assert c['prepared_only']
  pin(Path(__file__),c['worker_sha256'])
  for p,h in c['runtime_input_sha256'].items():pin(ROOT/p,h)
  d=json.loads((ROOT/c['source_input']).read_text());assert d['primitive_in_a_over_Fp_e_and_closure'] and d['remainder_mod32003']
  s['prime_trial_certificate']=[[div,P%div] for div in range(2,179)]
  assert 178*178<=P<179*179 and all(rem for div,rem in s['prime_trial_certificate'])
  assert d['remainder_mod32003']==24204 and d['remainder']%P==24204
  s['stage']='fixed_Fbar_input';s['Fbar']=[[Z for e in range(10)] for a in range(7)]
  ia=bi(17578125%P);s['alpha_inverse']=ia
  for (a,e),v in d['canonical_integer_coefficients']:s['Fbar'][a][e]=(bm(int(v)%P,ia),0)
  assert s['Fbar'][6]==[(55,0),(66,0)]+[Z]*8
  i55=bi(55);s['L_constant_inverse']=i55;f0=[bm(s['Fbar'][a][0][0],i55) for a in range(7)]
  assert f0==[5457,11750,555,10369,14038,10960,1];g4=[2668,18535,17561,16979,1]
  q,rem=bdiv(f0,g4);assert rem==[0]*4 and len(q)==3 and q[-1]==1 and bmul(g4,q)==f0
  s['seed_f0']=f0;s['seed_quartic']=g4;s['q']=q[:2];s['stage']='fixed_field_quartic_factor_algorithm';fp_roots=[];s['seed_Fp_roots']=fp_roots
  for x in range(P):
   v=1
   for j in range(3,-1,-1):
    s['cursor']=dict(stage='seed_quartic_root_enumeration',field_element=x,coefficient=j,accumulator=v);v=(bm(v,x)+g4[j])%P
   if v==0:fp_roots.append(x);assert len(fp_roots)<=4
  assert len(fp_roots)==4;prod=[1]
  for x in fp_roots:prod=bmul(prod,[(-x)%P,1])
  assert prod==g4
  disc=(bm(q[1],q[1])-bm(4,q[0]))%P;leg=bp(disc,(P-1)//2);assert disc and leg==P-1
  s['quadratic_seed_certificate']=dict(q=q,discriminant=disc,Legendre=leg,complete_quartic_product=prod)
  r0=[(x,0) for x in fp_roots]+[(0,1),((-q[1])%P,P-1)];assert len(set(r0))==6
  s['seed_six_root_residuals']=[kval(f0,r) for r in r0];assert s['seed_six_root_residuals']==[Z]*6;s['seed_six_roots']=r0
  s['roots']=[[r]+[Z]*9 for r in r0];derivative=[bm(a,s['Fbar'][a][0][0]) for a in range(1,7)];inverses=[];s['derivative_inverse_certificates']=[]
  for r in r0:
   v=kval(derivative,r);iv,cert=ki(v);inverses.append(iv);s['derivative_inverse_certificates'].append(cert)
  s['derivative_inverses']=inverses;s['stage']='six_root_jets_e10'
  for n in range(1,10):
   for i in range(6):
    s['cursor']=dict(stage='root_Newton_coefficient',order=n,root=i);error=evaluate(s['roots'][i],n+1)[n]
    value=kn(km(error,inverses[i]));s['roots'][i][n]=value;s['updates'].append(dict(order=n,root=i,error=error,value=value))
  s['root_full_residuals']=[evaluate(r,10) for r in s['roots']];assert all(all(x==Z for x in row) for row in s['root_full_residuals'])
  s['linears']=[[kn(x) for x in r]+[ONE]+[Z]*9 for r in s['roots']]
  s['subsets'][0]=[ONE]+[Z]*9;s['stage']='all64_subset_products'
  for mask in range(1,64):
   bit=(mask&-mask).bit_length()-1;old=mask^(1<<bit);s['cursor']=dict(mask=mask,previous=old,factor=bit)
   s['subsets'][mask]=pmul(s['subsets'][old],old.bit_count(),s['linears'][bit],1,10)
  whole=allocate(s['subsets'][63],6,1);assert whole==[x for row in s['Fbar'] for x in row]
  s['complete_six_factor_product_times_L_verified']=True;s['stage']='all62_complete_recombinations'
  for mask in range(1,63,2):
   complement=63^mask;da=mask.bit_count();db=6-da
   for b in (0,1):
    G=allocate(s['subsets'][mask],da,b);Q=allocate(s['subsets'][complement],db,1-b)
    s['candidate_work']=dict(mask=mask,complement=complement,LC_allocation=b,G=G,Q=Q)
    product=pmul(G,da,Q,db,19);s['candidate_work']['full_product']=product
    failure=None
    for a in range(7):
     for e in range(19):
      target=s['Fbar'][a][e] if e<10 else Z;diff=ks(product[19*a+e],target)
      if diff!=Z:failure=dict(mask=mask,LC_allocation=b,a=a,e=e,difference=diff);break
     if failure is not None:break
    if failure is None:
     s['matched_polynomial_factor_pair']=s['candidate_work'];finish(True);return
    s['failures'].append(failure);s['candidate_work']={};guard()
  assert len(s['failures'])==62;finish(False)
 except BaseException as e:
  s['failure']=repr(e);dump('failure-complete-prefix.json',dict(state=s,input_sha256=pins,no_retry=True))
  dump('report.json',dict(id=out.name,all_passed=False,prefix_only=True,no_irreducibility_claim=True,
   reason=repr(e),products=s['products'],peak_owned=peak,max_field_temporary_bits=bits,input_sha256=pins,
   evidence_sha256={'failure-complete-prefix.json':sha(out/'failure-complete-prefix.json')},no_retry=True))
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True);p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())
