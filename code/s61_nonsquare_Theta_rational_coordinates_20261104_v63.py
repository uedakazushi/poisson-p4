#!/usr/bin/env python3
"""Fixed rational Theta coordinates, scheme inverse and exact unit factors."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def add(p,q,c=Q(1)):
 r=p.copy()
 for m,v in q.items():
  r[m]=r.get(m,Q())+c*v
  if not r[m]:del r[m]
 return r
def sm(p,c):return{m:c*v for m,v in p.items()if c*v}
def mul(p,q):
 r={}
 for m,v in p.items():
  for n,w in q.items():
   k=tuple(x+y for x,y in zip(m,n));r[k]=r.get(k,Q())+v*w
 assert len(r)<50000
 return{k:v for k,v in r.items()if v}
def shift(p,s,c=Q(1)):return{tuple(x+y for x,y in zip(m,s)):c*v for m,v in p.items()if c*v}
def power(p,k,n):
 r={(0,)*n:Q(1)}
 for _ in range(k):r=mul(r,p)
 return r
def enc(p):return[[list(m),str(v)]for m,v in sorted(p.items())]
def dec(p):return{tuple(m):Q(v)for m,v in p if Q(v)}
def subs(p,images,n):
 r={};cache={}
 for m,c in p.items():
  t={(0,)*n:c}
  for i,k in enumerate(m):
   assert k>=0
   if(i,k)not in cache:cache[i,k]=power(images[i],k,n)
   t=mul(t,cache[i,k])
  r=add(r,t)
 return r
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 prior=ROOT/'reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63';paths=[Path(__file__),prior/'report.json',prior/'certificate.json']
 rd=lambda p:json.loads(p.read_text());pr,pc=[rd(p)for p in paths[1:]];assert sha(paths[2])==pr['certificate_sha256']
 a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={(0,0,0):Q(1)}
 k=add(sm(e,3),sm(one,-2));g1=dec(pc['g1']);g2=dec(pc['g2']);Theta=dec(pc['Theta'])
 C=add(sm(a,100),mul(k,add(sm(e,6),one)))
 B=add(mul(add(sm(a,50),mul(k,add(sm(e,6),one))),d),sm(k,20))
 F=add(add(add(sm(a,500),sm(power(e,2,3),90)),sm(e,3)),sm(one,-24));J=add(add(sm(mul(a,d),5),sm(e,-6)),sm(one,-4))
 N=add(mul(d,B),mul(g1,g2),Q(-1))
 assert g2==add(sm(g1,2),sm(d,5))
 assert add(power(N,2,3),sm(mul(g1,power(g2,2,3)),20),Q(-1))==mul(power(d,2,3),Theta)
 # Target ring QQ[z,d,d^-1], with precisely translated FULLU localizers.
 z={(1,0):Q(1)};D={(0,1):Q(1)};O={(0,0):Q(1)}
 zp20=add(z,sm(O,20));L=add(power(z,2,2),sm(D,50));K=add(mul(zp20,add(z,sm(O,-10))),sm(D,50))
 A=shift(mul(zp20,K),(0,-2),Q(1,500));E=add(sm(O,Q(2,3)),shift(add(power(z,2,2),sm(O,-400)),(0,-1),Q(1,60)))
 images=[A,E,D]
 G1=sm(power(z,2,2),Q(1,20));G2=sm(L,Q(1,10))
 BB=shift(mul(mul(L,z),zp20),(0,-1),Q(1,200));CC=shift(mul(L,power(zp20,2,2)),(0,-2),Q(1,200))
 mapped={name:subs(p,images,2)for name,p in [('a',a),('e',e),('d',d),('g1',g1),('g2',g2),('B',B),('C',C),('F',F),('J',J),('Theta',Theta),('N',N)]}
 for name,value in [('g1',G1),('g2',G2),('B',BB),('C',CC)]:assert mapped[name]==value
 assert not mapped['Theta'];assert mapped['N']==mul(z,G2)
 Fnum=shift(mapped['F'],(0,2),Q(40));Jnum=shift(mapped['J'],(0,1),Q(100))
 Fshort=add(add(mul(power(z,2,2),power(zp20,2,2)),mul(D,add(add(sm(power(z,2,2),82),sm(z,2000)),sm(O,7200)))),sm(power(D,2,2),720))
 Jshort=add(mul(power(z,2,2),zp20),sm(mul(D,add(z,sm(O,4))),50))
 assert Fnum==Fshort and Jnum==Jshort
 # Reverse composition is tested by polynomial numerator identities moduloTheta.
 monic=shift(Theta,(0,0,-2),Q(1,2500))
 def divtheta(f):
  r=f.copy();q={};steps=0
  while r and max(m[0]for m in r)>=2:
   h=max(m[0]for m in r);b={(m[0]-2,m[1],m[2]):v for m,v in r.items()if m[0]==h};q=add(q,b);r=add(r,mul(b,monic),Q(-1));steps+=1;assert steps<20
  q=shift(q,(0,0,-2),Q(1,2500));assert add(mul(q,Theta),r)==f;assert not r;return q
 reverse=[]
 def record(name,numerator,denominator,target):
  f=add(numerator,mul(denominator,target),Q(-1));q=divtheta(f)
  reverse.append({'coordinate':name,'substituted_numerator':enc(numerator),'substituted_denominator':enc(denominator),'difference_quotient_by_Theta':enc(q)})
 Gsq=power(g2,2,3);N2=power(N,2,3);Np20=add(N,sm(g2,20));Lnum=add(N2,sm(mul(d,Gsq),50))
 record('e',add(N2,mul(add(sm(d,40),sm(one,-400)),Gsq)),sm(mul(d,Gsq),60),e)
 Knum=add(add(N2,sm(mul(N,g2),10)),mul(add(sm(d,50),sm(one,-200)),Gsq))
 record('a',mul(Np20,Knum),sm(mul(power(d,2,3),power(g2,3,3)),500),a)
 record('B',mul(mul(Lnum,N),Np20),sm(mul(d,power(g2,4,3)),200),B)
 record('C',mul(Lnum,power(Np20,2,3)),sm(mul(power(d,2,3),power(g2,4,3)),200),C)
 record('g1',N2,sm(Gsq,20),g1)
 record('g2',Lnum,sm(Gsq,10),g2)
 # Regular source representative z=N*C*g1/B^2 equals N/g2 moduloTheta.
 sourceZnum=mul(mul(N,C),g1);sourceZden=power(B,2,3)
 assert add(mul(sourceZnum,g2),mul(N,sourceZden),Q(-1))==sm(mul(N,Theta),Q(-1))
 factor_record={'a':[enc(zp20),enc(K),'1/(500d^2)'],'B':[enc(L),enc(z),enc(zp20),'1/(200d)'],
                'C':[enc(L),enc(zp20),enc(zp20),'1/(200d^2)'],'F':[enc(Fnum),'1/(40d^2)'],'J':[enc(Jnum),'1/(100d)']}
 certificate={'source_ring':'R=QQ[a,e,d,1/(adFBCJ)]/(Theta)',
  'target_ring':'S=QQ[z,d,d^-1,1/(z*(z+20)*L*K*Fnum*Jnum)]',
  'source_Z_definition':'z=(dB-g1*g2)/g2','source_Z_regular_representative':{'numerator':enc(sourceZnum),'denominator':enc(sourceZden)},
  'squarecompletion':{'N':enc(N),'identity':'N^2-20*g1*g2^2=d^2*Theta'},
  'target_coordinates':{'L':enc(L),'K':enc(K),'a':enc(A),'e':enc(E),'B':enc(BB),'C':enc(CC),'g1':enc(G1),'g2':enc(G2),'Fnum':enc(Fnum),'Jnum':enc(Jnum)},
  'forward_all_images':{name:enc(p)for name,p in mapped.items()},'forward_source_Theta_zero':True,'forward_Z_composition_cleared_zero':True,
  'all_six_reverse_composition_quotients':reverse,'unit_localizer_factorizations':factor_record,
  'source_unit_proofs':{'g1_inverse':'C*g2/B^2','g2_inverse':'C*g1/B^2','z_inverse':'z*C*g2/(20B^2)',
   '(z+20)_inverse':'g2*(z+20)/(20d^2*C)','L_inverse':'C*g1/(10B^2)','K_inverse':'(z+20)/(500d^2*a)',
   'Fnum_inverse':'1/(40d^2*F)','Jnum_inverse':'1/(100d*J)'},
  'scheme_isomorphism_scope':'Mutual regular maps after exactly translated FULLU localizers, with literal polynomial identities moduloTheta. No point, irreducibility, reducedness, B-rank, minor, GB or physicalextension premise.',
  'other_pending_Theta_defects_not_computed':True}
 (out/'certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n')
 result={'id':ID,'literal_squarecompletion_verified':True,'all_forward_coordinate_images_verified':True,
         'mutual_regular_compositions_verified':True,'all_six_reverse_Theta_quotient_identities_verified':True,
         'FULLU_localizer_translation_explicit':True,'Fnum_terms':len(Fnum),'Jnum_terms':len(Jnum),
         'no_B_solver_minor_GB_point_prime_or_other_defect_test':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
