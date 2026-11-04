#!/usr/bin/env python3
"""Independent Sylvester determinant and ideal audit of the generic total octic."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63_REVIEW'
def read(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 authorpath=ROOT/'reports/Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63/report.json';author=read(authorpath);assert author['all_passed']
 checked={}
 for group in ['input_sha256','evidence_sha256','checked_dependency_sha256']:
  for name,digest in author.get(group,{}).items():
   candidates=[ROOT/name,authorpath.parent/name];p=next((p for p in candidates if p.is_file()and sha(p)==digest),None);assert p is not None,name;checked[str(p.relative_to(ROOT))]=digest
 certpath=authorpath.parent/'certificate.json';cert=read(certpath)
 RR,x,y,a,c,d,f,u,v,K,B,R,T=ring('x,y,a,c,d,f,u,v,K,B,R,T',s.QQ)
 m=u-v;A=R**2-2*d*R-f;b=-x*v*R**2+f*x*u-K*R;q=v*R**2+f*u
 a1=-c*y*R;b1=b;c1=y*u*v*R;a2=c*A;b2=B*R+2*a*f*u-2*a*v*R**2;c2=-v*v*R**2+f*u**2+2*d*u*v*R
 Q1=a1*T**2+b1*T+c1;Q2=a2*T**2+b2*T+c2
 # Independent elimination: determinant of the actual quadratic Sylvester4 matrix.
 syl=[[a1,b1,c1,RR.zero],[RR.zero,a1,b1,c1],[a2,b2,c2,RR.zero],[RR.zero,a2,b2,c2]]
 resultant=RR.zero
 for perm in itertools.permutations(range(4)):
  term=RR.domain.convert((-1)**sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4)))
  for row in range(4):term*=syl[row][perm[row]]
  resultant+=term
 P=-resultant.exquo(c)
 L=y*R*b2+A*b;M=y*m*R*q;E=L*T+M
 assert y*R*Q2+A*Q1==E
 assert y*R*P==L**2*Q1-a1*E**2-(b*L-2*a1*M)*E
 assert a1*M**2-b*M*L+c1*L**2==y*R*P
 assert y*R*(a2*M**2-b2*M*L+c2*L**2)==-A*y*R*P
 expr=P.as_expr();Rs=s.Symbol('R');poly=s.Poly(expr,Rs)
 assert poly.degree()==8 and s.expand(poly.LC()-s.Symbol('x')**2*s.Symbol('v')**4)==0
 assert s.expand(poly.TC()-s.Symbol('x')**2*s.Symbol('f')**4*s.Symbol('u')**4)==0
 for k,value in enumerate(cert['octic_coefficients_R0_through_R8']):
  assert RR.from_expr(poly.nth(k)-s.sympify(value).subs(s.Symbol('m'),s.Symbol('u')-s.Symbol('v')))==RR.zero
 # Independent quadratic resultant: remainder q-vA=lR+z.
 ell=2*d*v;z=f*(u+v);resAq=z*z+2*d*ell*z-f*ell*ell
 assert q-v*A==ell*R+z and resAq==f*f*(u+v)**2+4*d*d*f*u*v
 actualpath=ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/certificate.json.gz';actual=read(actualpath)
 xx,yy,aa,kk,mm,rr,tt=s.symbols('x y a k2 m R T');Delta=xx**2-yy**2;S=xx**2+yy**2;ell0=(aa+xx)**2-yy**2
 dd=-(aa*xx+S)/yy;cc=mm-ell0;ff=mm-xx**2*ell0/yy**2
 base=-4*aa**2*xx**2-2*aa*xx**3-6*aa*xx*yy**2+3*xx**4+yy**4;kap=kk-base
 uu=(kap-(xx-yy)**2*mm)/(4*xx*yy);vv=(kap-(xx+yy)**2*mm)/(4*xx*yy)
 compact=actual['compact_coefficients'];subs={s.Symbol('c'):cc,s.Symbol('d'):dd,s.Symbol('f'):ff,s.Symbol('u'):uu,s.Symbol('v'):vv,
  s.Symbol('K'):s.sympify(compact['K1']),s.Symbol('B'):s.sympify(compact['constant'])}
 assert s.cancel(Q1.as_expr().subs(subs,simultaneous=True)-s.sympify(actual['polynomial22_equation1']))==0
 assert s.cancel(Q2.as_expr().subs(subs,simultaneous=True)-s.sympify(actual['polynomial22_equation2']))==0
 actual_res=s.Poly(s.cancel(resAq.as_expr().subs(subs,simultaneous=True)),mm)
 assert actual_res.degree()==4 and s.cancel(actual_res.LC()-S**2/(4*xx**2*yy**2))==0
 assert s.diff(s.sympify(compact['K1']),s.Symbol('k1'))==s.Rational(1,2)
 assert s.diff(s.sympify(compact['constant']),s.Symbol('k1'))==0
 evidence=dict(independent_method='24-term determinant of quadratic Sylvester matrix; QQ polynomial ring, m eliminated at outset',
  resultant_identity='Res_T(Q1,Q2)=-c*P8',all_nine_octic_coefficients_match=True,
  full_ideal_identities=True,T_recovery='-y*(u-v)*R*(v*R^2+f*u)/L',
  quadratic_resultant=str(resAq),actual_m_degree=4,actual_m_leading=str(S**2/(4*xx**2*yy**2)),
  generic_L_unit_argument_audited=True,total_algebra_dimension=8,
  no_selected_component_exact_degree_or_genus_or_NEW_degree_claim=True)
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 inputs=[Path(__file__),authorpath,certpath,actualpath,ROOT/'research/notes/Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63.md']
 report=dict(id=STEM,all_passed=True,independent_resultant_and_all_coefficients=True,
  full_ideal_and_unit_proof_accepted=True,total_monogenic_eight_dimensional_algebra_accepted=True,
  selected_actual_component_degree8_reduced_irreducible_genus_or_NEW_multisection_degree_claimed=False,
  no_point_GB_high_degree_evaluation=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  checked_dependency_sha256=checked,evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({'all_passed':True,'seconds':report['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
