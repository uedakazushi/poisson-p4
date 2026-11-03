#!/usr/bin/env python3
"""Verify whole-open pair rank2 with saved literal QQ polynomial identities.

No GB, resultant search, factorization, rational reconstruction or finite-field
inference is performed. Every inverse is an explicit polynomial modulo C or Q23.
"""
from pathlib import Path
import sympy as s
import json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/pair-finite-subresultant-v78/whole-pair-cleanup-certificate.json'
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/pair-empty-QQ-replay-v78'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();t,H,v=s.symbols('t H v')
 with(OUT/'finite-source-replay.log').open('w')as log:
  subprocess.run([sys.executable,str(ROOT/'code/s61_nonsquare_pair_finite_bad_locus_20261103_v63.py'),str(OUT/'finite-source')],check=True,stdout=log)
 source=json.loads((OUT/'finite-source/certificate.json').read_text());data=json.loads(DATA.read_text())
 def load(terms,vars):return s.Poly(sum(s.Rational(c)*s.prod(x**k for x,k in zip(vars,ex))for ex,c in terms),*vars,domain=s.QQ)
 loadHV=lambda terms:load(terms,(H,v));loadv=lambda terms:load(terms,(v,))
 eqs=[load(q,(t,v,H))for q in source['linear_t_equations']];ab=[]
 for q in eqs:
  pp=s.Poly(q.as_expr(),t);ab.append((s.Poly(pp.coeff_monomial(t),H,v,domain=s.QQ),s.Poly(pp.coeff_monomial(1),H,v,domain=s.QQ)))
 N01,A,B=list(map(loadHV,source['elimination_polynomials']));R=loadHV(source['projection_certificate']['R']);f=(v+3)*(3*v+4)
 cubic=data['cubic_branch'];C=loadv(cubic['C']);large=data['degree23_branch'];Q=loadv(large['Q'])
 assert C.degree()==3 and Q.degree()==23
 quotient,remainder=R.div(s.Poly(f*f*C.as_expr()**2*Q.as_expr(),H,v,domain=s.QQ));assert remainder.is_zero and quotient.total_degree()==0 and not quotient.is_zero
 # C branch: t is a known unit, so beta3=-f*C forces alpha3=0.
 assert s.Poly(s.cancel(-ab[2][1].as_expr()/f),v,domain=s.QQ)==C
 RC=loadHV(cubic['alpha4_N01_resultant']);U,V=list(map(loadHV,cubic['resultant_Bezout_alpha4_N01']))
 assert U*ab[2][0]+V*N01==RC and RC.degree(H)==0
 RCv=s.Poly(RC.as_expr(),v,domain=s.QQ);RCinv=loadv(cubic['resultant_inverse_mod_C']);assert (RCv*RCinv).rem(C)==s.Poly(1,v,domain=s.QQ)
 # Q branch: a literal member of (A,B) is linear in H, with unit LC modulo Q.
 S1=loadHV(data['linear_subresultant']);UA,UB=list(map(loadHV,data['subresultant_A_B_Bezout']));assert UA*A+UB*B==S1
 linear=s.Poly(S1.as_expr(),H);assert linear.degree()==1;LC=s.Poly(linear.coeff_monomial(H),v,domain=s.QQ);cc=s.Poly(linear.coeff_monomial(1),v,domain=s.QQ)
 LCinv=loadv(large['linear_leading_inverse']);Hb=loadv(large['Hbar']);assert (LC*LCinv).rem(Q)==s.Poly(1,v,domain=s.QQ);assert (Hb+cc*LCinv).rem(Q).is_zero
 def evaluate(expr,Tb=None):
  poly=s.Poly(expr,t,H,v,domain=s.QQ);hp=[s.Poly(1,v,domain=s.QQ)];tp=[s.Poly(1,v,domain=s.QQ)];vp=[s.Poly(1,v,domain=s.QQ)]
  for _ in range(poly.degree(H)):hp.append((hp[-1]*Hb).rem(Q))
  for _ in range(poly.degree(v)):vp.append((vp[-1]*s.Poly(v,v,domain=s.QQ)).rem(Q))
  if Tb is not None:
   for _ in range(poly.degree(t)):tp.append((tp[-1]*Tb).rem(Q))
  result=s.Poly(0,v,domain=s.QQ)
  for ex,c in poly.terms():
   assert Tb is not None or ex[0]==0
   result+=(tp[ex[0]]*hp[ex[1]]*vp[ex[2]]).rem(Q)*c
  return result.rem(Q)
 alpha=evaluate(ab[2][0].as_expr());beta=evaluate(ab[2][1].as_expr());alpha_inv=loadv(large['alpha4_inverse']);Tb=loadv(large['Tbar']);assert (alpha*alpha_inv).rem(Q)==s.Poly(1,v,domain=s.QQ);assert (alpha*Tb+beta).rem(Q).is_zero
 assert evaluate(A.as_expr()).is_zero and evaluate(B.as_expr()).is_zero
 for al,be in ab:assert(evaluate(al.as_expr())*Tb+evaluate(be.as_expr())).rem(Q).is_zero
 rows=[[load(q,(t,v,H)).as_expr()for q in row]for row in source['natural_rows']]
 actual=evaluate(rows[0][5]*rows[1][0]-rows[1][5]*rows[0][0],Tb);saved=loadv(large['remaining_pair_minor_col0']);assert actual==saved
 inverse=loadv(data['degree23_remaining_minor_inverse']);assert(actual*inverse).rem(Q)==s.Poly(1,v,domain=s.QQ)
 result=dict(all_exact_QQ_checks_passed=True,whole_pair_rank2_on_a_P_J_open=True,open='a*P*J!=0',special_f_zero_excluded_by_frozen_finite_proof=True,R_factor_C3_Q23_case_cover_exact=True,cubic_branch_exact_polynomial_Bezout_exclusion=True,degree23_branch_H_t_unit_inverses_and_actual_minor_inverse=True,auxiliary_five_rank5_corollary_via_accepted_binary_D=True,raw24_unit_on_FULLU_corollary_via_accepted_gcd_and_binary_packets=True,raw24_corollary_open='a*P*F*B*C*J!=0',original17_kernel_bridge_not_reaudited_here=True,no_original_product41_full57_or_ambient_claim=True,no_finite_field_promotion_or_properness_assumption=True,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
