#!/usr/bin/env python3
"""Exact finite vertex-proportional / scroll-degeneration intersection.

The companion note proves integrality of the scroll base on E5!=0 by an explicit
P1 parameterization and records the reviewed quad78 finite-defect consequence.
"""
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import subprocess,json,sys,time
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/vertex-curve-E5-v85'
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/vertex-scroll-finite-replay-v85'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();a,e,P,t,H,v,z=s.symbols('a e P t H v z');load=lambda ts,vs:s.Poly(sum(s.Rational(c)*s.prod(x**k for x,k in zip(vs,ex))for ex,c in ts),*vs,domain=s.QQ)
 for stem,folder in [('s61_nonsquare_four_quadric_binary_reduction_20261103_v63','binary'),('s61_nonsquare_pair_empty_QQ_20261103_v63','pair')]:
  with(OUT/(folder+'.log')).open('w')as log:subprocess.run([sys.executable,str(ROOT/'code'/(stem+'.py')),str(OUT/folder)],check=True,stdout=log)
 binary=json.loads((OUT/'binary/certificate.json').read_text());pair=json.loads((OUT/'pair/finite-source/certificate.json').read_text());profile=json.loads((DATA/'profile.json').read_text());intersection=json.loads((DATA/'intersection-profile.json').read_text());cert=json.loads((DATA/'intersection-certificate.json').read_text());J=5*P-6*e-4
 cc=[[s.sympify(q)for q in row]for row in binary['C1_C2_four_coordinate_coefficients']];alpha,beta,gamma,delta=cc[0];eps,zeta,eta,theta=cc[1];assert s.cancel(beta-32*a**3*J**2/243)==0
 r1=s.Poly(s.cancel(P**2*(eps*beta-(zeta-40*z/6075)*(alpha-z/6075))),z);r2=s.Poly(s.cancel(P**2*(beta*(eta-theta*z)-(zeta-40*z/6075)*(gamma-delta*z))),z)
 for c in r1.all_coeffs()+r2.all_coeffs():s.Poly(c,a,e,P,domain=s.QQ)
 aa,bb=r1.all_coeffs(),r2.all_coeffs();assert len(aa)==len(bb)==3;M=s.Matrix([aa+[0],[0]+aa,bb+[0],[0]+bb]);DM=DomainMatrix.from_Matrix(M);res=DM.domain.to_sympy(DM.det());Eold=load(profile['E5_original'],(a,e,P))
 expected=-s.Rational(67108864,4749088069944112812890625)*a**9*J**6*P**2*Eold.as_expr();assert s.Poly(res-expected,a,e,P,domain=s.QQ).is_zero
 natural=s.Poly(s.cancel(Eold.as_expr().subs(a,t*P).subs(e,(v-20*t)/3)/P**3).subs(P,H-8*t),t,H,v,domain=s.QQ);E5=load(profile['E5_natural'],(t,H,v));assert natural==E5
 D=55*H-66*v+18;N=2*(v+2)*(7*v-11)**2;assert E5==s.Poly(N-t*D**2,t,H,v,domain=s.QQ)
 # Use two actual binary-vertex minors. No division by f(v) is made.
 chosen=[0,2]
 # The actual coefficient columns 2 and 4, together with column 5, are
 # precisely S^2, S*U, U^2.  Column 3 is NOT a vertex coefficient.
 assert [pair['linear_equation_normalizations'][i]['pair_coefficient_index']for i in chosen]==[2,4]
 equations=[load(pair['linear_t_equations'][i],(t,v,H)).reorder(t,H,v)for i in chosen];F=[]
 for eq,record,encoded in zip(equations,intersection['normalizations'],intersection['generic_projection_polynomials']):
  pp=s.Poly(eq.as_expr(),t);raw=s.Poly(s.expand(D**2*pp.coeff_monomial(1)+N*pp.coeff_monomial(t)),H,v,domain=s.QQ);assert record['removed_D_J_N_f']==[0,0,0,0];f=load(encoded,(H,v));assert raw==f*s.Rational(record['QQcontent']);F.append(f)
 assert F==[load(q,(H,v))for q in cert['projection_polynomials']]
 coefficients=[]
 for index,degree in cert['coefficient_labels']:coefficients.append(s.Poly(s.Poly(F[index].as_expr(),H).coeff_monomial(H**degree),v,domain=s.QQ))
 weights=[load(q,(v,))for q in cert['coefficient_Bezout']];assert sum((c*q for c,q in zip(weights,coefficients)),s.Poly(0,v,domain=s.QQ))==s.Poly(1,v,domain=s.QQ)
 special=[s.Poly(q.as_expr().subs(v,cert['good_special_v']),H,domain=s.QQ)for q in F];assert[s.degree()for s in special]==cert['special_H_degrees']==[q.degree(H)for q in F];u,w=[load(q,(H,))for q in cert['special_Bezout']];assert u*special[0]+w*special[1]==s.Poly(1,H,domain=s.QQ)
 for branch in cert['D_zero_branches']:
  vv=s.Rational(branch['v']);hh=s.Rational(branch['H']);assert D.subs({v:vv,H:hh})==0 and N.subs(v,vv)==0
  actual=[s.Poly(q.as_expr().subs({v:vv,H:hh}),t,domain=s.QQ)for q in equations];assert actual==[load(q,(t,))for q in branch['polynomials']];u,w=[load(q,(t,))for q in branch['Bezout']];assert u*actual[0]+w*actual[1]==s.Poly(1,t,domain=s.QQ)
 result=dict(all_exact_QQ_actual_source_checks_passed=True,scroll_rankdrop_resultant_exact=True,E5_natural='2*(v+2)*(7*v-11)^2-t*(55H-66v+18)^2',scroll_base_geometrically_integral_on_a_P_J_E5_open=True,actual_vertex_coefficient_columns=[2,4,5],vertex_binary_rankdrop_intersection_E5_dimension_at_most_zero=True,projection_finite_on_D_nonzero=True,D_zero_branches_empty_by_QQ_Bezout=True,no_f_chart_or_unknown_unit_division=True,quad78_finite_defect_corollary_uses_reviewed_G1_G2_and_corrected_criterion=True,quad78_corollary_open='a*P*J!=0',original_light_quadratic_quartic_rank146_outside_finite_corollary=True,original_rank146_corollary_open='a*P*F*B0*C*J!=0',all_point_quad78_or_rank146_lower_not_proved=True,no_full_quartic_or_ambient_claim=True,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
