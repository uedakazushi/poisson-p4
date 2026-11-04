#!/usr/bin/env python3
"""Exact generic full-PGL stabilizer: Ga x mu2; no points or GB."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(M):return [[str(s.factor(x))for x in row]for row in M.tolist()]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(55,55));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(55);out.mkdir(parents=True,exist_ok=False)
 geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
 sig=ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63/certificate.json'
 review=ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63_REVIEW/report.json'
 G=json.loads(geo.read_text());B=json.loads(sig.read_text());assert G['all_passed']and json.loads(review.read_text())['all_passed']
 u,v,a,c,x,y,z,w,H,K,Z,T=s.symbols('u v a c x y z w H K Z T');p=32*u+v;coords=(x,y,z,w);n0,n1,n2=s.symbols('n0 n1 n2')
 P=s.Matrix(G['affine_P']).subs({n0:u,n1:1,n2:v});N=s.Matrix(G['homogeneous_modular_N']).subs({n0:u,n1:1,n2:v})
 def Ga(t):return s.eye(5)+t*N+t*t*N*N/2
 Ms=[s.Matrix(b['homogeneous_matrix'])for b in B['branches']]
 R=(Ms[1].inv()*Ms[0]).applyfunc(s.cancel);t0=-1/(128*u*v)
 assert(R*N-N*R).applyfunc(s.cancel)==s.zeros(5)
 assert(R*R-Ga(t0)).applyfunc(s.cancel)==s.zeros(5)
 E=(Ga(-t0/2)*R).applyfunc(s.cancel)
 expected=s.Matrix([[1,0,0,0,0],[-4/(3*v),-1,0,0,0],[4/(3*v*v),2/v,1,0,0],[-2/(3*v**3),-1/v**2,0,1,0],[s.Rational(1,1536),0,0,0,-1]])
 assert E==expected
 assert(E*E-s.eye(5)).applyfunc(s.cancel)==s.zeros(5)
 assert(E*N-N*E).applyfunc(s.cancel)==s.zeros(5)
 D=a*a-64*c*c;kap=1/D;j=D*D;e=(D-a)/v;g=-64*c/v;F=x+v*y
 Pyz=u*v*y*y/2+u*x*y+w*w-x*x/64+64*w*y
 xp=a*x+64*c*w;yp=e*x+D*y+g*w;wp=c*x+a*w
 qdiff=s.Poly(s.expand(s.cancel(kap*j*(e*(-64*w*F)+D*Pyz+g*(-x*F))-(u*v*yp*yp/2+u*xp*yp+wp*wp-xp*xp/64+64*wp*yp))),x,y,w)
 Q1=-p*D*D+p*a*a+4096*a*c-64*v*c*c
 Q2=-D*D+a*a+u*a*c+64*c*c
 Q3=v*D*D-v*a*a+4096*a*c+64*p*c*c
 for mon,expr in{(2,0,0):Q1/(64*v),(1,0,1):64*Q2/v,(0,0,2):Q3/v}.items():assert s.cancel(qdiff.coeff_monomial(s.prod(q**b for q,b in zip((x,y,w),mon)))-expr)==0
 assert len(qdiff.terms())==3
 assert s.expand(Q1+Q3+(v-p)*Q2+32*a*c*(u*u-256))==0
 assert s.expand(Q3+v*Q2-a*c*(4096+u*v)-64*c*c*(p+v))==0
 assert s.expand(Q2.subs(c,0)-a*a*(1-a*a))==0
 branches=[]
 for aa in[1,-1]:
  tx=2*(aa-1)/(3*v);tw=(1-aa)/s.Integer(3072);ii=s.Rational(3,4)*(tx+v*T)
  A=s.Matrix([[aa,0,0,0],[(1-aa)/v,1,0,0],[H,ii,1,K],[0,0,0,aa]])
  tr=s.Matrix([tx,T,Z,tw]);Y=A*s.Matrix(coords)+tr
  diff=(A*P*A.T-P.subs(dict(zip(coords,Y)),simultaneous=True)).applyfunc(s.expand)
  equations=[];labels=[]
  for i in range(4):
   for k in range(i+1,4):
    for mon,co in s.Poly(diff[i,k],*coords).terms():
     co=s.cancel(co)
     if co:equations.append(co);labels.append([i,k,*mon])
  LM,rhs=s.linear_eq_to_matrix(equations,[H,K,Z]);rr=list(LM.T.rref()[1]);minor=s.factor(LM[rr[:3],:].det());assert minor!=0
  ans=[3*T/4,0,3*v*T*T/8+T/(64*u)]if aa==1 else[(3*T*v*v-8)/(4*v*v),0,(72*T*T*u*v**3-192*T*u*v+3*T*v*v-4)/(192*u*v*v)]
  sub=dict(zip([H,K,Z],ans));assert all(s.cancel(q.subs(sub))==0 for q in equations)
  MM=s.eye(5);MM[1:,0]=tr;MM[1:,1:]=A;MM=MM.subs(sub)
  timepar=3*T/(4*u)if aa==1 else s.Rational(3,4)/u*(T-4/(3*v*v))
  desired=Ga(timepar)if aa==1 else Ga(timepar)*E
  assert(MM-desired).applyfunc(s.cancel)==s.zeros(5)
  branches.append(dict(a=aa,c=0,full_linear_equations=list(map(str,equations)),coefficient_labels=labels,
    linear_matrix=rows(LM),rhs=list(map(str,rhs)),rank3_minor_rows=rr[:3],rank3_minor=str(minor),
    unique_H_K_Z=list(map(str,ans)),all_full_coefficient_identities_zero=True,
    homogeneous_matrix=rows(MM),Ga_parameter=str(timepar),component='Ga'if aa==1 else'Ga*E'))
 cert=dict(generic_parameter_field='k(u,v), k=QQ or C; geometric generic exhaustion',
  intrinsic_flags='unique tripleH, vertex of Q|H, unique compatible constant bivector line C=dy wedge dz',
  initial_general_affine_flag='xprime=a*x+b*w+tx; wprime=c*x+d*w+tw; yprime=e*x+f*y+g*w+ty; zprime=h*x+i*y+j*z+k*w+tz',
  necessary_bracket_equations=['kappa*(a*d-b*c)=1','eta=kappa*f','eta*b=64*c','eta*a=d','eta*d=a','eta*c=b/64'],
  modular_N_y_equation='kappa*f*u=u hence eta=1 for u!=0',
  resulting_flag=dict(b='64*c',d='a',D='a^2-64*c^2',kappa='1/D',f='D',j='D^2',e='(D-a)/v',g='-64*c/v',
   tx='-c/48+2*(a-1)/(3*v)',tw='(1-a)/3072+2*c/(3*v)',i='3*D*(tx+v*T)/4+D*(1-D)/(64*u)'),
  quadratic_coefficients=dict(x2=str(Q1/(64*v)),xw=str(64*Q2/v),w2=str(Q3/v)),
  quadratic_combination='Q1+Q3+(v-p)*Q2=-32*a*c*(u^2-256)',
  second_combination='Q3+v*Q2=a*c*(4096+u*v)+64*c^2*(p+v)',
  generic_exhaustion='u^2-256 and p+v units => ac=0 => c=0; Dunit => a!=0; Q2=a^2(1-a^2) => a=+1 or-1',
  R=rows(R),R_square_Ga_time=str(t0),R_commutes_N=True,N=rows(N),
  E=rows(E),E_squared_identity=True,E_commutes_N=True,E_outside_Ga_by_x_w_linear_block=True,
  all_full_flag_branches=branches,
  exact_generic_projective_Poisson_stabilizer='Ga x mu2',
  explicit_group_map='(q,epsilon) -> exp(q*N)*E^epsilon, epsilon=0,1',
  no_all_special_parameter_stabilizer_classification=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),geo,sig,review]
 report=dict(id=ID,all_passed=True,actual_involution_and_all_generic_stabilizer_cases_verified=True,
  exact_generic_full_PGL5_stabilizer='Ga x mu2 over k(u,v), geometrically too',
  no_new_point_prime_GB_or_unbounded_solver=True,
  no_special_parameter_group_exhaustion=True,independent_review_pending=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','exact_generic_full_PGL5_stabilizer','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
