#!/usr/bin/env python3
"""ONE named-root U-saturation, with every parameter/evaluation unit certified."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_NAMED_SATURATION_20261104_v63';sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(a):int(v)for a,v in p}
def pa(a,b,c=1):
 z=dict(a)
 for i,v in b.items():
  z[i]=z.get(i,0)+c*v
  if not z[i]:del z[i]
 return z
def pm(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():z[i+j]=z.get(i+j,0)+c*d
 return{k:v for k,v in z.items()if v}
class Stop(Exception):
 def __init__(self,status,data):self.status=status;self.data=data
class Ring:
 def __init__(self,g):self.g=g;self.L=g[13];self.operations=0
 def make(self,p,d=1):
  if not p:return ({},1)
  if d<0:p={i:-v for i,v in p.items()};d=-d
  while p and max(p)>=13:
   k=max(p)-13;c=p[max(p)];p=pa({i:v*self.L for i,v in p.items()},{i+k:c*v for i,v in self.g.items()},-1);d*=self.L
  if not p:return({},1)
  gcd=math.gcd(d,*(abs(v)for v in p.values()));return({i:v//gcd for i,v in p.items()},d//gcd)
 def c(self,n):return self.make({0:n})if n else({},1)
 def add(self,x,y):
  self.operations+=1;D=math.lcm(x[1],y[1]);return self.make(pa({i:v*(D//x[1])for i,v in x[0].items()},{i:v*(D//y[1])for i,v in y[0].items()}),D)
 def neg(self,x):return({i:-v for i,v in x[0].items()},x[1])
 def sub(self,x,y):return self.add(x,self.neg(y))
 def mul(self,x,y):self.operations+=1;return self.make(pm(x[0],y[0]),x[1]*y[1])
 def scale(self,x,n,d=1):return self.make({i:v*n for i,v in x[0].items()},x[1]*d)
 def enc(self,x):return dict(numerator=[[i,str(v)]for i,v in sorted(x[0].items(),reverse=True)],denominator=str(x[1]))
 def dec(self,r):return self.make(ip(r['numerator']),int(r['denominator']))
 def eval(self,p,r):
  x=self.c(0)
  for c in reversed(p):x=self.add(self.mul(x,r),c)
  return x
 def synth(self,p,r):
  n=len(p)-1;q=[self.c(0)for _ in range(n)];q[-1]=p[-1]
  for i in range(n-2,-1,-1):q[i]=self.add(p[i+1],self.mul(r,q[i+1]))
  assert self.add(p[0],self.mul(r,q[0]))==self.c(0)
  product=[self.c(0)for _ in range(n+1)]
  for i,c in enumerate(q):product[i]=self.sub(product[i],self.mul(r,c));product[i+1]=self.add(product[i+1],c)
  assert product==p;return q
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 ub=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_U_INPUT_20261104_v63';hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63';bridge=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW';assert sha(ub/'report.json')=='ff55a3ca80fb7dee461faef458c9d866de8db8773de5e66c0278d30270d2a424';assert read(bridge/'report.json')['all_passed']
 h=read(hb/'H4-monic-field-coefficients.json');u=read(ub/'U-field-coefficients.json');g={a:int(v)for a,v in h['g13']};ring=Ring(g)
 from fractions import Fraction as Q
 def load(p):
  z={a:Q(v)for a,v in p};D=math.lcm(*(v.denominator for v in z.values()));return ring.make({a:int(v*D)for a,v in z.items()},D)
 H=list(map(load,h['H4_coefficients_d0_to4']));U=list(map(load,u['U_coefficients_d0_to3']));t=ring.make(ip(h['e_numerator']),int(h['e_denominator']));a=ring.make({1:1});one=ring.c(1)
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp';exe=work/'unit-check'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 units=[];cache={};events=[];removals=[];parameters={};remaining=H
 def saveprogress(status,extra=None):
  d=dict(status=status,parameters={k:ring.enc(v)for k,v in parameters.items()},units=units,events=events,removals=removals,remaining_H=[ring.enc(c)for c in remaining],extra=extra,elapsed_seconds=time.monotonic()-start)
  payload=json.dumps(d,indent=2)+'\n';assert len(payload.encode())<28*1024*1024;(out/'progress.json').write_text(payload)
 def inv(tag,x):
  if not x[0]:raise Stop('zero_parameter_or_evaluation',dict(tag=tag,input=ring.enc(x)))
  key=(tuple(sorted(x[0].items())),x[1])
  if key in cache:
   j,v=cache[key];events.append(dict(reuse_unit=j,tag=tag));return v
  content=math.gcd(*(abs(v)for v in x[0].values()));p={i:v//content for i,v in x[0].items()};m=max(p);B=max(abs(v).bit_length()for v in p.values());bg=max(abs(v).bit_length()for v in g.values());height=13*(B+math.ceil(math.log2(m+1)))+m*(bg+4)+math.factorial(m+13).bit_length()+10;bound=(m+14)*(150+(height*30103+99999)//100000)
  existing=sum(p.stat().st_size for p in out.iterdir()if p.is_file());assert bound<28*1024*1024
  if existing+bound+2*1024*1024>32*1024*1024:raise Stop('inverse_output_size_guard',dict(tag=tag,predicted_unit_bytes=bound,existing_bytes=existing))
  j=len(units);inp=work/f'unit-{j}.txt';result=work/f'unit-{j}.json'
  with inp.open('w')as f:
   for z in[g,p]:
    f.write(str(len(z))+'\n')
    for i,v in sorted(z.items(),reverse=True):f.write(str(i)+' '+str(v)+'\n')
  subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','named-unit',str(exe),str(inp),str(result)],check=True,timeout=max(1,56-(time.monotonic()-start)))
  n=read(result)
  if not n['unit']:raise Stop('proper_base_component_retained',dict(tag=tag,input=ring.enc(x),gcd=n['gcd']))
  S,T=ip(n['S']),ip(n['T']);R=int(n['R']);assert pa(pm(S,p),pm(T,g))=={0:R};v=ring.make({i:c*x[1]for i,c in S.items()},content*R);assert ring.mul(x,v)==one
  rec=dict(tag=tag,input=ring.enc(x),primitive_numerator=[[i,str(c)]for i,c in sorted(p.items(),reverse=True)],numerator_content=str(content),S=n['S'],T=n['T'],R=n['R'],inverse=ring.enc(v),unit_certificate_bound_bytes=bound);units.append(rec);cache[key]=(j,v);saveprogress('running');return v
 status='pending';extra={}
 try:
  t2=ring.mul(t,t);F=ring.add(ring.add(ring.scale(a,500),ring.scale(t2,90)),ring.add(ring.scale(t,3),ring.c(-24)));C=ring.add(ring.add(ring.scale(a,100),ring.scale(t2,18)),ring.add(ring.scale(t,-9),ring.c(-2)));b1=ring.add(ring.add(ring.scale(a,50),ring.scale(t2,18)),ring.add(ring.scale(t,-9),ring.c(-2)));b0=ring.add(ring.scale(t,60),ring.c(-40));parameters.update(a=a,t=t,F=F,C=C,b1=b1,b0=b0)
  ia=inv('a',a)
  if not F[0]or not C[0]:
   assert all(not x[0]for x in U);status='FULLU_closed_cover_EMPTY_parameter_factor_zero';extra=dict(zero_factor='F'if not F[0]else'C')
  else:
   iF=inv('F',F);iC=inv('C',C);ib1=inv('b1',b1);rB=ring.neg(ring.mul(b0,ib1));rJ=ring.scale(ring.mul(ring.add(ring.scale(t,6),ring.c(4)),ia),1,5);lead=ring.scale(ring.mul(ring.mul(ring.mul(ring.mul(a,a),F),C),b1),5);parameters.update(rB=rB,rJ=rJ,U_leading_unit=lead)
   factor=[ring.c(0),ring.mul(lead,ring.mul(rB,rJ)),ring.neg(ring.mul(lead,ring.add(rB,rJ))),lead];assert factor==U;events.append(dict(U_factorization_all4_coefficients_exact=True,factor_coefficients=[ring.enc(v)for v in factor]))
   for tag,root in [('0',ring.c(0)),('B',rB),('J',rJ)]:
    while len(remaining)>1:
     value=ring.eval(remaining,root)
     if not value[0]:
      previous=remaining;remaining=ring.synth(previous,root);removals.append(dict(root_tag=tag,root=ring.enc(root),previous_H=[ring.enc(v)for v in previous],quotient_H=[ring.enc(v)for v in remaining],all_factor_coefficients_exact=True));assert len(removals)<=4;saveprogress('running')
     else:
      inv('remaining_eval_at_'+tag,value);events.append(dict(root_tag=tag,remaining_degree=len(remaining)-1,nonroot_eval=ring.enc(value)));break
   if len(remaining)==1:
    assert remaining==[one];status='FULLU_closed_cover_EMPTY_all_factors_removed'
   else:
    for tag,root in [('0',ring.c(0)),('B',rB),('J',rJ)]:inv('FINAL_remaining_eval_at_'+tag,ring.eval(remaining,root))
    status='clopen_U_unit_cover_retained';extra=dict(remaining_degree=len(remaining)-1,QQ_length=13*(len(remaining)-1),no_actual_Bbad_existence=True)
 except Stop as e:status=e.status;extra=e.data
 saveprogress(status,extra);complete=status.startswith('FULLU_closed_cover_EMPTY')or status=='clopen_U_unit_cover_retained';empty=status.startswith('FULLU_closed_cover_EMPTY')
 proof=dict(status=status,all_parameter_and_nonzero_evaluation_inverses_Bezout_certified=complete,complete_source_U_factorization=complete,all_monic_synthetic_divisions_exact=complete,remaining_H_degree=len(remaining)-1,removed_multiplicities=len(removals),scheme_scope='Ka13[d]/H4 localizedatFULLU; namedfactorremovalpreserveslocalizedideal overANYcoefficientring. NoKa13fieldoretaleguess.',empty_closed_cover_proved=empty,remaining_cover_U_is_unit=complete and not empty,no_global_Bbad_claim=True,no_dEuclid_12x28_Gaussian_Bminor_or_newnorm=True)
 proof['explicit_inverse_recipe_when_retained']='For eachnamedrootr, H_U(d)=(d-r)*Q_r(d)+H_U(r) byuniversalmonicdivision. ModH_U, inverse(d-r)=-Q_r*certifiedinverse(H_U(r)). Multiplythree andcertifiedinverse(Uleadingunit). No expandedinversepolynomial computed.'
 (out/'scheme-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',ub/'report.json',ub/'U-field-coefficients.json',ub/'known-factor-saturation-preflight.json',hb/'report.json',hb/'H4-monic-field-coefficients.json',bridge/'report.json']
 report=dict(id=ID,all_passed=complete,status=status,FULLU_closed_cover_empty=empty,remaining_degree=len(remaining)-1,removed_multiplicities=len(removals),certified_unit_queries=len(units),field_integer_operations=ring.operations,no_field_guess_dEuclid_Gaussian_minor_norm_point_prime=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(e),'progress_preserved':(a.output/'progress.json').exists(),'no_auto_retry':True},indent=2)+'\n')
  raise
