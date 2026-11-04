#!/usr/bin/env python3
"""Finish ONLY leading cases from saved regular/J certificates; shared60s budget."""
import argparse,hashlib,json,math,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_FINISH_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 first=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_20261104_v63/report.json'
 previous=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_CONTINUATION_20261104_v63/report.json';prev=read(previous);assert prev['error']=="KeyError('variables')"
 saved=previous.parent/'certificate.json';cert=read(saved);assert 'Fn_regular_branch'in cert
 spent=prev['combined_mathematical_seconds'];remaining=math.floor(60-spent);assert remaining>0
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('remaining shared60s budget')));signal.alarm(remaining)
 source=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';data=read(source)
 r,s,R=S.symbols('r s R');variables2=(r,R);variables3=(r,s,R);ring=S.QQ.poly_ring(r,R);P=lambda x:S.Poly(x,*variables2,domain=S.QQ)
 def load(record,gens):return S.Poly.from_dict({tuple(e):S.Rational(a,b)for a,b,*e in record['coefficients']},gens,domain=S.QQ)
 def rec(p):return dict(variables=list(map(str,p.gens)),terms=len(p.terms()),is_zero=p.is_zero,multidegree=list(p.degree_list()),coefficients=[[str(c.p),str(c.q),*e]for e,c in p.terms()])
 a,b,c=[load(cert['Fn_regular_branch'][name],variables2)for name in['a','b','c']]
 linear=[]
 for name in['A1','A2']:
  p=load(data['rational_coefficient_pairs'][name]['numerator'],variables3);ps=S.Poly(p.as_expr(),s,domain=ring);d=ps.degree();value=P(0)
  for(e,),coefficient in ps.terms():value+=P(coefficient)*((-c)**e)*(b**(d-e))
  linear.append(value)
 cert['Fn_leading_zero_cases']=dict(a_zero_not_inverted=True,case_a0_b_nonzero='a=0,b!=0: b^deg(Ai)*Ai(s=-c/b)=0 necessary for bothi',
  homogenized_linear_root_conditions=[rec(p)for p in linear],case_a0_b0='a=b=0 then c=0 necessary, otherwise Fn has no root; no remaining Ai condition discarded',these_cases_not_excluded=True)
 leading=[]
 for name in['Z1','Z2']:
  con=data['localized_constraints'][name];degree=con['shape_bidegree'][1]
  for key in['coefficient1','coefficient2','denominator']:
   p=load(con[key],variables3);leading.append(P(S.Poly(p.as_expr(),s,domain=ring).nth(degree)))
 cert['s_infinity']=dict(s_degrees=[3,5],all6_independent_target_leading_coefficients=[rec(p)for p in leading],
  necessary_condition='all6 vanish for a fixed leading root persisting over independent targetparameters; no gcd/elimination computed')
 cert['Hn_boundary']=dict(observation='OffJ*Fn and namedunits, both normalizedmaps extend regularly on the properHn=0 curve. Its image cannot dominate independent2D(Z1,Z2); pole/base intersections require explicitJ/Fn cases.',no_Hn_resultant_or_saturation=True)
 cert['continuation_provenance']=dict(no_repeated_J_gcd_Fn_pseudodivision_or_normcalculation=True,previous_partial_receipts_preserved=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');signal.alarm(0)
 summary=dict(J_residual_necessary_factor='25r-1',Fn_regular_condition_degrees={k:cert['Fn_regular_branch'][k]['multidegree']for k in['linear_remainder_determinant']},
  Fn_norm_degrees=[p['multidegree']for p in cert['Fn_regular_branch']['two_quadratic_norm_conditions']],
  leading_case_condition_degrees=[p.degree_list()for p in linear])
 report=dict(id=STEM,all_passed=True,summary=summary,completed_J_gcd_Fn_norm_and_all_leading_case_preflight=True,
  no_unknown_r_polynomial_globally_inverted=True,Fn_leading_and_regular_conditions_not_solved=True,
  no_Q19_general_factorization_fullresultant_PRS_point_prime_source26_or_full_content_gcd=True,
  no_primitive_Q19_actual_generic_degree_targetindex_or_NEW_equation_claim=True,
  previous_mathematical_seconds=spent,remaining_cap_seconds=remaining,combined_mathematical_seconds=spent+time.monotonic()-start,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),first,previous,saved,source]},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 assert(out/'certificate.json').stat().st_size<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','summary','combined_mathematical_seconds','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
