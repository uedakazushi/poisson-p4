#!/usr/bin/env python3
"""Independent Laurent arithmetic, actual quartic coefficients, and genus-bound proof audit."""
import argparse,gzip,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from ambient_B18_relative_weighted_model_review_20261104_v63 import parse
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_RECIPROCAL_GENUS_20261104_v63_REVIEW'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def add(*aa):
 d={}
 for a in aa:
  for e,c in a.items():d[e]=d.get(e,Q(0))+c
 return {e:c for e,c in d.items()if c}
def scale(a,c):return {e:v*c for e,v in a.items()if v*c}
def mul(a,b):
 d={}
 for e,c in a.items():
  for f,v in b.items():
   g=tuple(x+y for x,y in zip(e,f));d[g]=d.get(g,Q(0))+c*v
 return {e:c for e,c in d.items()if c}
def shift(a,e,c=1):return {tuple(x+y for x,y in zip(f,e)):v*c for f,v in a.items()if v*c}
def power(a,n):
 b={(0,)*len(next(iter(a))):Q(1)}
 for _ in range(n):b=mul(b,a)
 return b
def variables(names):return [{tuple(int(i==j)for i in range(len(names))):Q(1)}for j in range(len(names))]
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 primary=ROOT/'reports/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63';genus=ROOT/'reports/Q_FIRSTJET_RESIDUAL_GENUS_BOUND_20261104_v63/report.json'
 author=read(primary/'report.json');ac=read(primary/'certificate.json');gr=read(genus);assert author['all_passed']and gr['all_passed']
 for folder,report in [(primary,author),(genus.parent,gr)]:
  for rel,h in report['input_sha256'].items():assert sha(ROOT/rel)==h,rel
  for rel,h in report.get('evidence_sha256',{}).items():assert sha(folder/rel)==h,rel
 old=ROOT/'reports/Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63/certificate.json';oc=read(old)
 names=['a','c','d','f','u','v','R','T','K','B','x','y','r','t','z','D']
 a,c,d,f,u,v,R,T,K,B,x,y,r,t,z,D=variables(names);n=len(names)
 def mon(**ee):return {tuple(ee.get(name,0)for name in names):Q(1)}
 uv=mul(u,v);mm=add(u,scale(v,-1))
 rr=add(scale(mul(v,R),-1),mul(mul(f,u),mon(R=-1)))
 tt=add(scale(mul(c,T),-1),mul(uv,mon(T=-1)))
 hh=add(mul(v,R),mul(mul(f,u),mon(R=-1)))
 jj=add(mul(c,T),mul(uv,mon(T=-1)));zz=mul(hh,jj)
 L1=add(mul(x,rr),mul(y,tt),scale(K,-1))
 L2=add(B,mul(c,mul(R,T)),mul(f,mul(power(u,2),mon(R=-1,T=-1))),
        scale(mul(mul(c,f),mul(T,mon(R=-1))),-1),scale(mul(power(v,2),mul(R,mon(T=-1))),-1),
        scale(mul(mul(c,d),T),-2),scale(mul(mul(a,mul(f,u)),mon(R=-1)),2),scale(mul(mul(a,v),R),-2),scale(mul(mul(d,uv),mon(T=-1)),2))
 assert mul(mul(R,T),L1)==parse(oc['auxiliary_polynomials']['Q1'],names)
 assert mul(mul(R,T),L2)==parse(oc['auxiliary_polynomials']['Q2'],names)
 def iota(p):
  result={}
  for e,co in p.items():
   er,et=e[6],e[7];ex=list(e);ex[6]=-er;ex[7]=-et
   # R -> -fu/(vR), T -> -uv/(cT).
   ex[3]+=er;ex[4]+=er+et;ex[5]+=-er+et;ex[1]-=et
   ee=tuple(ex);result[ee]=result.get(ee,Q(0))+co*(-1 if(er+et)%2 else 1)
  return {e:c for e,c in result.items()if c}
 assert iota(L1)==L1 and iota(L2)==L2
 assert iota(rr)==rr and iota(tt)==tt and iota(zz)==zz and iota(hh)==scale(hh,-1)
 qr=add(power(rr,2),scale(mul(f,uv),4));qt=add(power(tt,2),scale(mul(c,uv),4))
 assert power(hh,2)==qr and power(jj,2)==qt and power(zz,2)==mul(qr,qt)
 assert add(scale(mul(uv,add(L2,scale(B,-1),scale(mul(a,rr),-2),scale(mul(d,tt),-2))),2),
            scale(mul(mul(add(u,v),rr),tt),-1),scale(mul(mm,zz),-1))=={}
 # Reconstruct coordinates without dividing by a non-monomial h.
 assert add(scale(mul(v,R),2),scale(hh,-1),rr)=={}
 assert add(scale(mul(mul(c,T),hh),2),scale(zz,-1),mul(tt,hh))=={}
 W=add(scale(mul(uv,D),2),mul(mul(add(u,v),r),t))
 psi=add(mul(power(r,2),power(t,2)),mul(mul(mul(add(u,v),r),t),D),mul(uv,power(D,2)),
         scale(mul(power(mm,2),add(mul(c,power(r,2)),mul(f,power(t,2)),scale(mul(mul(c,f),uv),4))),-1))
 assert add(power(W,2),scale(mul(power(mm,2),mul(add(power(r,2),scale(mul(f,uv),4)),add(power(t,2),scale(mul(c,uv),4)))),-1),scale(mul(uv,psi),-4))=={}
 # Independent DIRECT Laurent substitution of all actual coefficients,
 # unlike the author's hand-cleared bar-polynomial computation.
 anames=['x','y','a','k1','k2','r','m'];x,y,a,k1,k2,r,m=variables(anames)
 def amon(**ee):return {tuple(ee.get(name,0)for name in anames):Q(1)}
 S=add(power(x,2),power(y,2));delta=add(power(x,2),scale(power(y,2),-1))
 ell=add(power(add(a,x),2),scale(power(y,2),-1))
 d=scale(mul(add(mul(a,x),S),amon(y=-1)),-1);c=add(m,scale(ell,-1));f=add(m,scale(mul(mul(power(x,2),ell),amon(y=-2)),-1))
 b2=add(scale(mul(power(a,2),power(x,2)),-4),scale(mul(a,power(x,3)),-2),scale(mul(mul(a,x),power(y,2)),-6),scale(power(x,4),3),power(y,4))
 u=scale(mul(add(k2,scale(b2,-1),scale(mul(power(add(x,scale(y,-1)),2),m),-1)),amon(x=-1,y=-1)),Q(1,4))
 v=scale(mul(add(k2,scale(b2,-1),scale(mul(power(add(x,y),2),m),-1)),amon(x=-1,y=-1)),Q(1,4))
 assert add(u,scale(v,-1))==m
 ay_dx=add(mul(a,y),mul(d,x));uv=mul(u,v)
 A0=add(mul(mul(a,d),ay_dx),mul(ay_dx,add(u,v)),mul(mul(c,d),y),mul(mul(f,a),x))
 K1=add(scale(k1,Q(1,2)),scale(A0,-1))
 B0=add(mul(add(power(a,2),c),add(power(d,2),f)),power(m,2),scale(mul(mul(a,d),add(u,v)),2),scale(power(delta,2),-1))
 tt=mul(add(K1,scale(mul(x,r),-1)),amon(y=-1))
 D0=add(B0,scale(mul(a,r),2),scale(mul(d,tt),2))
 actual=add(mul(power(r,2),power(tt,2)),mul(mul(mul(add(u,v),r),tt),D0),mul(uv,power(D0,2)),
            scale(mul(power(m,2),add(mul(c,power(r,2)),mul(f,power(tt,2)),scale(mul(mul(c,f),uv),4))),-1))
 plane=shift(actual,(2,10,0,0,0,0,0),16)
 assert len(plane)==1592 and all(min(e)>=0 and co.denominator==1 for e,co in plane.items())
 grouped={}
 for key,terms in ac['actual_primitive_plane_coefficients_r_m'].items():
  er,em=map(int,key.split(','))
  for ex,c0 in terms:
   e=tuple(ex)+(er,em);assert e not in grouped;grouped[e]=Q(c0)
 assert grouped==plane
 assert max(e[5]+e[6]for e in plane)==4 and len({e[5:]for e in plane})==15
 assert {e[:5]:co for e,co in plane.items()if e[5]==4}=={(4,8,0,0,0):Q(16)}
 norm=add(power(B0,2),scale(mul(power(a,2),mul(f,uv)),16))
 assert max(e[6]for e in norm)==4
 assert {e[:6]:co for e,co in norm.items()if e[6]==4}=={(0,0,0,0,0,0):Q(4)}
 q=add(power(r,2),scale(mul(f,uv),4));assert max(e[5]+e[6]for e in q)==3
 # Matrix interpretation of the reciprocal map; square roots are formal
 # scalars and disappear from the resulting R,T rational formulas.
 aa,dd,cc,ff,gg,hh,ii,jj,xx,yy,ss,ts=s.symbols('a d c f g h i j x y s t',nonzero=True)
 C=s.Matrix([[aa,1,gg,hh],[cc,-aa,ii,jj],[-jj,hh,dd,1],[ii,-gg,ff,-dd]])
 D=s.diag(xx,-xx,yy,-yy);J2=s.Matrix([[0,1],[-1,0]]);J4=s.diag(J2,J2);tor=s.diag(ss,1/ss,ts,1/ts)
 Cnew=tor*J4*(-C)*J4.inv()*tor.inv();Dnew=tor*J4*(-D)*J4.inv()*tor.inv()
 assert Dnew==D
 images=[-ss/ts*jj,ss*ts*ii,hh/(ss*ts),-ts/ss*gg]
 assert [s.cancel(Cnew[i,j]-v)for(i,j),v in zip([(0,2),(0,3),(1,2),(1,3)],images)]==[0]*4
 assert(Cnew.subs({cc:1/ss**2,ff:1/ts**2})-s.Matrix([[aa,1,*images[:2]],[1/ss**2,-aa,*images[2:]],[-images[3],images[1],dd,1],[images[2],-images[0],1/ts**2,-dd]])).applyfunc(s.cancel)==s.zeros(4)
 table=[]
 for degree in range(1,5):
  genus0=(degree-1)*(degree-2)//2;poles=3*degree;table.append({'plane_component_degree':degree,'plane_genus_upper':genus0,'pole_divisor_degree_upper':poles,'ramification_degree_upper':2*poles,'connected_quadratic_genus_upper':2*genus0-1+poles})
 assert table==gr['degree_table']and max(v['connected_quadratic_genus_upper']for v in table)==17
 proof={'all1592_actual_plane_coefficients_independently_verified':True,'independent_actual_clearing':'F=16*x^2*y^10*Psi',
        'actual_quartic_total_degree':4,'r4_coefficient':'16*x^4*y^8','q_total_degree_r_m':3,
        'same_original_ideal_binding':'R*T*L1=Q1, R*T*L2=Q2 in the frozen octic presentation',
        'involution':'R->-fu/(vR), T->-uv/(cT), m fixed; r,t,z fixed and h changes sign',
        'matrix_orientation':'(C,D)->(-C,-D), then diag(J2,J2) conjugation and diagonal normalization; C signs/offdiagonal positions independently checked',
        'generic_base_units':['x','y','x^2-y^2','m','c','f','u','v'],
        'h_unit_proof':'If Psi=q=0, r is algebraic over QQ(x,y,a,k2,m) while t is transcendental because K1=k1/2-A0. The linear W=0 forces B0+2ar=0, hence B0^2+16a^2fuv=0, impossible since its m^4 coefficient is4.',
        'scheme_inverse':'t=(K1-xr)/y, z=-(2uvD0+(u+v)rt)/m, h^2=q, R=(h-r)/(2v), T=(z/h-t)/(2c)',
        'inverse_units':'1/R=(h+r)/(2fu), 1/T=(z/h+t)/(2uv); h^-1=h*q^-1; q^-1 may be computed by adjugate of its4x4 multiplication matrix',
        'all_algebras_not_only_points':'quartic base is free4 (leading x^2/y^2 unit), monic quadratic is free2 and h unit gives etale; h-sign fixed algebra is exactly the quartic base even if nonreduced',
        'genus_proof':'For every m-dominant geometric component, normalize its plane component of degree d<=4. Genus<=3, q poles<=3d<=12, odd zero/pole valuations at most6d<=24. RH gives connected genus<=17; split genus<=3.',
        'genus_degree_table':table,
        'scope':'ordered eigenvalue base and existing relative-sign quotient; no irreducibility/smoothness/exact genus/selected degree8 or Weyl descent assumption'}
 ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'code/ambient_B18_relative_weighted_model_review_20261104_v63.py',primary/'report.json',primary/'certificate.json',genus,old,
         ROOT/'research/notes/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63.md',ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW/report.json',ROOT/'reports/Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63_REVIEW/report.json']
 report={'id':STEM,'all_passed':True,'independent_Laurent_involution_and_all_universal_quotient_identities':True,
         'literal_original_Q1_Q2_ideal_generators_bound':True,'independent_direct_Laurent_actual_coefficient_substitution_matches_all1592_terms':True,
         'quartic_total_degree4_quadratic_cover_q_degree3':True,'generic_h_unit_and_full_scheme_free4_times_free2_inverse_proof_accepted':True,
         'independent_matrix_negative_pair_and_reciprocal_orientation_checked':True,'every_m_dominant_geometric_component_normalization_genus_upper':17,
         'no_irreducibility_smoothness_exact_genus_selected_degree8_or_Weyl_descent_claim':True,
         'no_author_rerun_beta_point_new_invariant_GB_octic_factorization_or_specialization':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','every_m_dominant_geometric_component_normalization_genus_upper','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
