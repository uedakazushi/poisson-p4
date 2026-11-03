#!/usr/bin/env python3
"""One bounded two-variable scalar-row GB/lift on entry primitives.

Known-unit valuations and QQ contents are retained, so every selected primitive
lift can be restored and checked against the original polynomial row.
"""
import argparse,gzip,json,math,pathlib,resource,subprocess,time
from fractions import Fraction as F
import sympy as sp
from sympy.polys.rings import ring
from s61_nonsquare_J0_degree1_identity_20261103_v63 import poly,sha,encoded,identity
ROOT=pathlib.Path(__file__).resolve().parents[1]
R,a,e=ring('a,e',sp.QQ);units=[a,500*a+90*e**2+3*e-24,300*a*e+54*e**3+9*e**2-24*e-4,100*a+18*e**2-9*e-2,3*e+2]
def express(p):
 if not p:return '0'
 return '('+'+'.join('('+str(c)+')'+('a^'+str(ex[0])if ex[0]else '')+('*e^'+str(ex[1])if ex[1]else '')for ex,c in sorted(p.items()))+')'
def m2expr(p):
 if not p:return '0'
 terms=[]
 for ex,c in sorted(p.items()):
  s='('+str(c)+')'
  for name,k in zip(['a','e'],ex):
   if k:s+='*'+name+'^'+str(k)
  terms.append(s)
 return '('+'+'.join(terms)+')'
def read_sparse(p):
 s=p.read_text().split();at=0;n=int(s[at]);at+=1;polys=[]
 for j in range(n):
  nt=int(s[at]);at+=1;q={}
  for k in range(nt):
   c=F(s[at]);ex=tuple(map(int,s[at+1:at+4]));at+=4;assert ex[2]==0;q[ex]=c
  polys.append(q)
 assert at==len(s);return polys
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--row',type=pathlib.Path,required=True);ap.add_argument('--output',type=pathlib.Path,required=True);ap.add_argument('--resume-saved-gb',action='store_true');args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=args.resume_saved_gb);start=time.monotonic()
 with args.row.open()as f:
  n,m,k=map(int,f.readline().split());assert n==1 and k==0;raw=[poly(f.readline())for _ in range(m)];assert not f.read().strip()
 base=[];valuations=[];scalars=[]
 for p in raw:
  q=R.from_dict({ex[:2]:sp.QQ(c.numerator,c.denominator)for ex,c in p.items()});v=[]
  if not q:base.append(q);valuations.append([0]*5);scalars.append(F(1));continue
  for unit in units:
   power=0
   while True:
    quotient,remainder=divmod(q,unit)
    if remainder:break
    q=quotient;power+=1
   v.append(power)
  den=math.lcm(*(int(c.denominator)for c in q.values()));num=math.gcd(*(int(c.numerator)*(den//int(c.denominator))for c in q.values()));scalar=F(num,den)
  b=q.exquo(R(sp.QQ(scalar.numerator,scalar.denominator)));assert all(c.denominator==1 for c in b.values())
  reconstructed=b*R(sp.QQ(scalar.numerator,scalar.denominator))
  for u,power in zip(units,v):reconstructed*=u**power
  assert reconstructed==R.from_dict({ex[:2]:sp.QQ(c.numerator,c.denominator)for ex,c in p.items()});base.append(b);valuations.append(v);scalars.append(scalar)
 metadata={'row_sha256':sha(args.row),'columns':m,'known_units':[str(u)for u in units],'valuations':valuations,'QQ_scalar_contents':list(map(str,scalars)),'all_original_reconstruction_identities_verified':True}
 (out/'primitive-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
 export='exportPolynomial=(out,p)->(tt:=terms p;out<<#tt;scan(tt,t->(ee:=first exponents t;out<<" "<<toString leadCoefficient t<<" "<<ee#0<<" "<<ee#1<<" 0";));out<<endl;);\n'
 script='R=QQ[a,e,MonomialOrder=>GRevLex];\nM=matrix{{'+','.join(map(m2expr,base))+'}};\nprint "PRIMITIVE_SCALAR_INPUT_READY";<<flush;\nG=gens gb M;print("PRIMITIVE_GB_DONE",numColumns G,dim(R/ideal G));<<flush;\n'+export
 script+='outG=openOut '+json.dumps(str(out/'GB.txt'))+';outG<<numColumns G<<endl;scan(flatten entries G,p->exportPolynomial(outG,p));close outG;\n'
 script+='assert(numColumns G<=20);L=G//M;assert(M*L==G);print("EXACT_ALL_GB_LIFT_PASS",sum(flatten entries L,p->#terms p));<<flush;\n'
 script+='outL=openOut '+json.dumps(str(out/'lifts.txt'))+';outL<<numRows L*numColumns L<<endl;scan(toList(0..numColumns L-1),j->scan(toList(0..numRows L-1),i->exportPolynomial(outL,L_(i,j))));close outL;\n'
 (out/'run.m2').write_text(script)
 def limits():resource.setrlimit(resource.RLIMIT_AS,(1800000000,1800000000));resource.setrlimit(resource.RLIMIT_CPU,(60,61))
 t=time.monotonic()
 if args.resume_saved_gb:
  assert'EXACT_ALL_GB_LIFT_PASS'in(out/'M2.log').read_text()and(out/'GB.txt').exists()and(out/'lifts.txt').exists();status='completed'
 else:
  with(out/'M2.log').open('w')as f:
   try:p=subprocess.run(['M2','--script',str(out/'run.m2')],stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=65);status='completed'if p.returncode==0 else f'exit_{p.returncode}'
   except subprocess.TimeoutExpired:status='wall_timeout'
 report={'status':status,'solver_seconds':None if args.resume_saved_gb else time.monotonic()-t,'saved_gb_reused_without_solver_rerun':args.resume_saved_gb,'preprocessing_seconds':t-start,'AS_limit_bytes':1800000000,'CPU_limit_seconds':60,'original_row_sha256':sha(args.row)}
 if status=='completed'and'EXACT_ALL_GB_LIFT_PASS'in(out/'M2.log').read_text():
  gs=read_sparse(out/'GB.txt');lifts=read_sparse(out/'lifts.txt');assert len(lifts)==m*len(gs);results=[]
  for j,g in enumerate(gs):
   q=lifts[j*m:(j+1)*m];bb=[{(*ex,0):F(str(c))for ex,c in b.items()}for b in base];assert identity(bb,q)==g
   powers=[max((v[u]for v,cf in zip(valuations,q)if cf),default=0)for u in range(5)]
   lifted=[];common=R.one
   for u,power in zip(units,powers):common*=u**power
   for cf,v,c in zip(q,valuations,scalars):
    if not cf:lifted.append({});continue
    h=R.from_dict({ex[:2]:sp.QQ(x.numerator,x.denominator)for ex,x in cf.items()})/R(sp.QQ(c.numerator,c.denominator))
    for u,top,old in zip(units,powers,v):h*=u**(top-old)
    lifted.append({(*ex,0):F(str(x))for ex,x in h.items()})
   target=identity(raw,lifted);expected=R.from_dict({ex[:2]:sp.QQ(c.numerator,c.denominator)for ex,c in g.items()})*common
   assert target=={(*ex,0):F(str(c))for ex,c in expected.items()};remaining=expected;targetpowers=[]
   for u in units:
    power=0
    while True:
     quotient,remainder=divmod(remaining,u)
     if remainder:break
     remaining=quotient;power+=1
    targetpowers.append(power)
   unit=bool(remaining and remaining.is_ground)
   entry={'GB_index':j,'global_known_unit':unit,'pure_e_primitive':all(ex[0]==0 for ex in g),'active':sum(bool(cf)for cf in lifted),'terms':sum(map(len,lifted)),'degree':max((sum(ex)for cf in lifted for ex in cf),default=0),'unit_factor_exponents':targetpowers if unit else None}
   results.append(entry)
   with gzip.open(out/f'original-lift-{j}.json.gz','wt')as f:json.dump({'raw_row':list(map(encoded,raw)),'coefficients':list(map(encoded,lifted)),'target':encoded(target),'result':entry},f,separators=(',',':'))
  report.update(all_original_QQ_lifts_verified=True,GB_generators=len(gs),results=results)
 report['seconds']=time.monotonic()-start;(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
