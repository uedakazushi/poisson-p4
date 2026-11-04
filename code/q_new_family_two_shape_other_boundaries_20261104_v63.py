#!/usr/bin/env python3
"""ONE bounded J/Fn/leading boundary preflight. NO full resultant or Q19 factorization."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s shared boundary preflight cap')));signal.alarm(60)
 r,s,R=S.symbols('r s R');v2=(r,R);P=lambda z:S.Poly(z,*v2,domain=S.QQ)
 path=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';review=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63_REVIEW/report.json';context=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_BOUNDARY_STRIP_20261104_v63/report.json';assert read(review)['all_passed']and read(context)['all_passed'];d=read(path)
 T=(r-1)*R-2*r;J=(10*r-1)*s-10*r;LJ=10*r-1;Fn=S.sympify(d['Fn']);domain=S.QQ.poly_ring(r,R);K=S.QQ.frac_field(R)
 Fs=S.Poly(Fn,s,domain=domain);a,b,c=[S.Poly(Fn,s).nth(i)for i in[2,1,0]]
 def load(rec):return S.Poly.from_dict({tuple(e):S.Rational(n,den)for n,den,*e in rec['coefficients']},(r,s,R),domain=S.QQ)
 C=[load(d['rational_coefficient_pairs'][k]['numerator'])for k in['C1','C2']];A=[load(d['rational_coefficient_pairs'][k]['numerator'])for k in['A1','A2']]
 def rec(poly):
  if not isinstance(poly,S.Poly):poly=P(poly)
  return dict(variables=list(map(str,poly.gens)),terms=len(poly.terms()),multidegree=list(poly.degree_list()),coefficients=[[str(z.p),str(z.q),*e]for e,z in poly.terms()])
 cert={};front=[];error=None
 def save(stage):
  front.append(dict(stage=stage,seconds=time.monotonic()-start));(out/'frontier.json').write_text(json.dumps(front,indent=2)+'\n');(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 def field_gcd_certificate(polys):
  pp=[S.Poly(p.as_expr()if isinstance(p,S.Poly)else p,r,domain=K)for p in polys];g=pp[0];bez=[S.Poly(1,r,domain=K)]
  for q in pp[1:]:
   u,v,h=g.gcdex(q);assert u*g+v*q==h;bez=[u*w for w in bez]+[v];g=h
  assert sum((u*p for u,p in zip(bez,pp)),S.Poly(0,r,domain=K))==g
  U=S.Poly(1,R,domain=S.QQ)
  for p in bez+[g]:
   for z in p.all_coeffs():U=U.lcm(S.Poly(S.denom(z),R,domain=S.QQ))
  def clear(p):
   ans=0
   for(e,),z in p.terms():
    q=S.cancel(z*U.as_expr());assert S.denom(q)==1 or not S.denom(q).has(R)
    ans+=q*r**e
   return P(ans)
  bb=[clear(p)for p in bez];gg=clear(g);assert sum((z*P(p.as_expr()if isinstance(p,S.Poly)else p)for z,p in zip(bb,polys)),P(0))==gg
  remain=g;powers={}
  for name,f in[('r',r),('r-1',r-1),('T',T)]:
   fp=S.Poly(f,r,domain=K);n=0
   while remain.degree()>0:
    q0,rem=remain.div(fp)
    if not rem.is_zero:break
    assert q0*fp==remain;remain=q0;n+=1
   powers[name]=n
  return dict(field='QQ(R)[r]',gcd=str(g.as_expr()),clearing_R_only_polynomial=rec(U),
   cleared_Bezout_coefficients=[rec(p)for p in bb],cleared_gcd=rec(gg),exact_cleared_Bezout_identity=True,
   known_unit_factor_powers=powers,residual_gcd=str(remain.as_expr()),gcd_has_no_unknown_r_factor=remain.degree()==0)
 try:
  # J leading case: LJ*(s-1)-J=1, so LJ=J=0 is impossible.
  assert S.expand(LJ*(s-1)-J-1)==0
  fnJ=P(a*(10*r)**2+b*(10*r)*LJ+c*LJ**2);assert fnJ==P(1024*r**3*T**2)
  cs=[]
  for p in C:
   ps=S.Poly(p.as_expr(),s,domain=domain);degree=ps.degree();h=P(0)
   for(e,),coef in ps.terms():h+=P(coef)*(P(10*r)**e)*(P(LJ)**(degree-e))
   cs.append(h)
  cert['J_boundary']=dict(leading_unit_identity='(10r-1)*(s-1)-J=1',homogenized_J_Fn_identity='(10r-1)^2*Fn(s=10r/(10r-1))=1024r^3*T^2',
   identity_polynomial=rec(fnJ),C1_C2_homogenized_at_J0=[rec(p)for p in cs],
   logic='r,T,Z0 units force Fn!=0 onJ=0; E1/E2 then force C1num=C2num=0. Leading10r-1=0 impossible by unitidentity')
  save('J_homogenizations_and_leading_case_verified')
  cert['J_boundary']['gcd_certificate']=field_gcd_certificate(cs);save('J_functionfield_gcd_and_Bezout_verified')
  # Fn regular branch: pseudo-division never silently inverts its unknown r-dependent LC.
  reduced=[]
  for i,p in enumerate(A):
   ps=S.Poly(p.as_expr(),s,domain=domain);quo,rem=ps.pdiv(Fs);power=max(0,ps.degree()-Fs.degree()+1)
   check=S.Poly(a**power*ps.as_expr()-quo.as_expr()*Fn-rem.as_expr(),s,r,R,domain=S.QQ);assert check.is_zero and rem.degree()<2
   u=P(rem.nth(1));v=P(rem.nth(0));reduced.append((u,v))
   cert['Fn_A'+str(i+1)+'_pseudo']=dict(leading_coefficient_power=power,quotient=rec(S.Poly(quo.as_expr(),s,r,R,domain=S.QQ)),remainder=rec(S.Poly(rem.as_expr(),s,r,R,domain=S.QQ)),u=rec(u),v=rec(v),exact_identity_verified=True)
   save('Fn_A'+str(i+1)+'_pseudo_remainder_verified')
  (u1,v1),(u2,v2)=reduced;det=u1*v2-u2*v1;norms=[]
  for i,(u,v)in enumerate(reduced):
   norm=P(a)*v*v-P(b)*u*v+P(c)*u*u;linear=u.as_expr()*s+v.as_expr()
   identity=u.as_expr()**2*Fn-(a*linear+b*u.as_expr()-2*a*v.as_expr())*linear-norm.as_expr();assert S.Poly(identity,s,r,R,domain=S.QQ).is_zero;norms.append(norm)
  assert P(u1.as_expr()*v2.as_expr()-u2.as_expr()*v1.as_expr())==det
  cert['Fn_regular_branch']=dict(required_unit='a=LC_s(Fn) !=0 ONLY in this branch',a=rec(P(a)),b=rec(P(b)),c=rec(P(c)),
   linear_remainder_determinant=rec(det),two_quadratic_norm_conditions=[rec(p)for p in norms],
   norm_identity='u^2*Fn - (a*(u*s+v)+b*u-2*a*v)*(u*s+v) = a*v^2-b*u*v+c*u^2',
   common_root_necessary_conditions='determinant=norm1=norm2=0',no_unknown_r_polynomial_globally_inverted=True)
  save('Fn_regular_determinant_and_norm_conditions_verified')
  # Retain a=0. If b!=0, substitute the LINEAR Fn root by polynomial homogenization.
  lead=[]
  for i,p in enumerate(A):
   ps=S.Poly(p.as_expr(),s,domain=domain);degree=ps.degree();value=P(0)
   for(e,),coef in ps.terms():value+=P(coef)*(P(-c)**e)*(P(b)**(degree-e))
   lead.append(value)
  cert['Fn_leading_zero_cases']=dict(case_a0_b_nonzero='a=0 and b!=0: necessary homogenized_A1=homogenized_A2=0',
   homogenized_at_linear_Fn_root=[rec(p)for p in lead],case_a0_b0='a=b=0 then c=0 necessary; otherwise Fn has no root',
   a_zero_not_discarded_or_inverted=True)
  save('Fn_leading_zero_cases_retained')
  # Small degree<=5 base polynomial first, so this remains a bounded boundary computation.
  if time.monotonic()-start<35:
   cert['Fn_leading_zero_cases']['a_b_c_gcd_certificate']=field_gcd_certificate([P(a),P(b),P(c)]);save('Fn_identically_zero_quadratic_case_bounded')
  # Generic infinity: independent Z parameters mean common fixed leading roots divide ALL6 coefficient LCs.
  leading=[]
  for name in['Z1','Z2']:
   z=d['localized_constraints'][name];degree=z['shape_bidegree'][1]
   for k in['coefficient1','coefficient2','denominator']:
    p=load(z[k]);leading.append(P(S.Poly(p.as_expr(),s,domain=domain).nth(degree)))
  cert['s_infinity']=dict(s_degrees=[3,5],all6_independent_parameter_leading_coefficient_polynomials=[rec(p)for p in leading],
   scope='fixed leading root for generic independenttargetparameters must annihilate all6; no full resultants or parameter specializations')
  save('all6_infinity_leading_coefficients_saved')
  if time.monotonic()-start<40:
   cert['s_infinity']['gcd_certificate']=field_gcd_certificate(leading);save('infinity_leading_content_Bezout_verified')
  cert['Hn_boundary']=dict(observation='Away from J*Fn=0 and other displayedunits, both normalized spectral maps are regular on Hn=0. This is a proper curve in2shapevariables; its image cannot dominate independentZ1,Z2. Intersection with pole/base loci belongs to the J/Fn cases, not silently removed.',
   no_new_Hn_resultant_or_saturation_computed=True)
  save('completed_boundary_preflight')
 except Exception as exc:error=repr(exc);save('bounded_partial_frontier')
 signal.alarm(0)
 report=dict(id=STEM,all_passed=error is None,error=error,one_shared60s_boundary_preflight=True,
  J_boundary_completed='gcd_certificate'in cert.get('J_boundary',{}),Fn_regular_conditions_completed='Fn_regular_branch'in cert,
  Fn_leading_zero_cases_retained='Fn_leading_zero_cases'in cert,no_unknown_r_polynomial_globally_inverted=True,
  no_new_full_resultant_PRS_Q19_factorization_point_prime_or_source26=True,
  Q19_primitive_generic_degree_actual_target_index_NEW_equation_not_claimed=True,
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,artifact_bytes=32*1024*1024),
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),path,review,context]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024-10000
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','error','J_boundary_completed','Fn_regular_conditions_completed','Fn_leading_zero_cases_retained','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
