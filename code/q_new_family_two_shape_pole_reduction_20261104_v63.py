#!/usr/bin/env python3
"""ONE tiny tracked-Fn remainder test and conditional exact cofactor reduction."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_POLE_REDUCTION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 r,s,R,Z0,Z1,Z2=S.symbols('r s R Z0 Z1 Z2');vars=(r,s,R);P=lambda x:S.Poly(x,*vars,domain=S.QQ)
 path=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';reportpath=path.parent/'report.json';assert json.loads(reportpath.read_text())['all_passed'];data=json.loads(path.read_text())
 pairs=data['rational_coefficient_pairs'];Fn=P(S.sympify(data['Fn']));Dc=P(r*(r-1)**2*((10*r-1)*s-10*r))
 def numerator(name):
  p=pairs[name];return S.Rational(p['constant'])*S.Poly.from_dict({tuple(e):S.Rational(a,b)for a,b,*e in p['numerator']['coefficients']},vars,domain=S.QQ)
 A=numerator('A1');B=numerator('A2');C=numerator('C1');D=numerator('C2');A2=A*A
 qB,remB=B.div(Fn);qA,remA=A2.div(Fn);assert qB*Fn+remB==B and qA*Fn+remA==A2
 keys=sorted(set(remB.monoms())|set(remA.monoms()));first=next(m for m in keys if remA.coeff_monomial(m));k=remB.coeff_monomial(first)/remA.coeff_monomial(first)
 difference=remB-k*remA;found=difference.is_zero
 def record(p):return dict(terms=len(p.terms()),multidegree=list(p.degree_list()),coefficients=[[str(c.p),str(c.q),*m]for m,c in p.terms()])
 result=dict(variables=list(map(str,vars)),tracked_factor='Fn',Fn=record(Fn),A1_numerator_with_constant=record(A),A2_numerator_with_constant=record(B),
  remainder_A1_squared=record(remA),remainder_A2=record(remB),candidate_constant=str(k),constant_pole_cancellation_found=found)
 summary=None
 if found:
  quotient,rem=(B-k*A2).div(Fn);assert rem.is_zero and quotient*Fn==B-k*A2
  # On E1=0, replace the old quintic E2 by a new equation with only ONE Fn pole.
  z1,z2,z0=Z1,Z2,Z0
  e1=A.as_expr()*Dc.as_expr()+z0*C.as_expr()*Fn.as_expr()-z1*Fn.as_expr()*Dc.as_expr()
  old=B.as_expr()*Dc.as_expr()+z0*D.as_expr()*Fn.as_expr()**2-z2*Fn.as_expr()**2*Dc.as_expr()
  new=quotient.as_expr()*Dc.as_expr()+k*z1*A.as_expr()*Dc.as_expr()+z0*(D.as_expr()*Fn.as_expr()-k*C.as_expr()*A.as_expr())-z2*Fn.as_expr()*Dc.as_expr()
  allvars=(r,s,R,Z0,Z1,Z2);Pe1=S.Poly(e1,*allvars,domain=S.QQ);Pold=S.Poly(old,*allvars,domain=S.QQ);Pnew=S.Poly(new,*allvars,domain=S.QQ)
  identity=S.Poly(Fn.as_expr()*new-old+k*A.as_expr()*e1,*allvars,domain=S.QQ);assert identity.is_zero
  assert len(Pnew.terms())<=20000
  summary=dict(old_shape_bidegree=list(Pold.degree_list()[:2]),new_shape_bidegree=list(Pnew.degree_list()[:2]),new_R_degree=Pnew.degree(R),new_terms=len(Pnew.terms()),
   paired_shape_bidegrees=[list(Pe1.degree_list()[:2]),list(Pnew.degree_list()[:2])],Sylvester_in_s_size=Pe1.degree(s)+Pnew.degree(s),
   resultant_r_degree_upper=Pe1.degree(r)*Pnew.degree(s)+Pnew.degree(r)*Pe1.degree(s))
  result.update(exact_quotient=record(quotient),cofactor_identity='Fn*E2_new = E2_old - k*A1_numerator*E1',
   E1=dict(variables=list(map(str,allvars)),**record(Pe1)),E2_old=dict(variables=list(map(str,allvars)),**record(Pold)),E2_new=dict(variables=list(map(str,allvars)),**record(Pnew)),
   equivalence_scope='previous dense unitopen, including Fn!=0; no equation elimination or global extension',new_equation_summary=summary)
 else:
  witness=next(m for m in keys if difference.coeff_monomial(m));minor=remA.coeff_monomial(first)*remB.coeff_monomial(witness)-remB.coeff_monomial(first)*remA.coeff_monomial(witness);assert minor
  result.update(no_QQ_constant_can_cancel_the_double_Fn_pole=True,coefficient_two_row_witness=[list(first),list(witness)],coefficient_minor=str(minor))
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n');assert(out/'certificate.json').stat().st_size<=32*1024*1024;signal.alarm(0)
 rep=dict(id=STEM,all_passed=True,one_tiny_tracked_polynomial_remainder_job=True,constant_k=str(k)if found else None,constant_pole_cancellation_found=found,
  exact_cofactor_reduction_verified=found,new_equation_summary=summary,no_gcd_GB_resultant_source26_newpoints=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),path,reportpath]},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k]for k in['all_passed','constant_k','constant_pole_cancellation_found','new_equation_summary','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
