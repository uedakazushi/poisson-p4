#!/usr/bin/env python3
"""Finish from saved verified J/pseudo data under REMAINING shared60s budget."""
import argparse,hashlib,json,math,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_CONTINUATION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 prefixdir=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_20261104_v63';prefixreport=prefixdir/'report.json';prefixpath=prefixdir/'certificate.json'
 prior=read(prefixreport);prefix=read(prefixpath);assert prior['error']=="TypeError('Value after * must be an iterable, not Poly')"and prior['J_boundary_completed']
 assert 'Fn_A1_pseudo'in prefix and 'Fn_A2_pseudo'in prefix
 remaining=math.floor(60-prior['seconds']);assert remaining>0
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('remaining shared60s budget')));signal.alarm(remaining)
 source=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';data=read(source)
 r,s,R=S.symbols('r s R');variables2=(r,R);P=lambda z:S.Poly(z,*variables2,domain=S.QQ);Fn=S.sympify(data['Fn']);F=S.Poly(Fn,s);a,b,c=[P(F.nth(i))for i in[2,1,0]]
 def load(rec):return S.Poly.from_dict({tuple(e):S.Rational(n,d)for n,d,*e in rec['coefficients']},tuple(S.Symbol(z)for z in rec['variables']),domain=S.QQ)
 def rec(p):return dict(variables=list(map(str,p.gens)),terms=len(p.terms()),is_zero=p.is_zero,multidegree=list(p.degree_list()),coefficients=[[str(z.p),str(z.q),*e]for e,z in p.terms()])
 cert=dict(J_boundary=prefix['J_boundary'],all_saved_pseudo_identities_reused_without_recomputation=True);front=[];error=None
 def save(stage):
  front.append(dict(stage=stage,continuation_seconds=time.monotonic()-start));(out/'frontier.json').write_text(json.dumps(front,indent=2)+'\n');(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 try:
  (u1,v1),(u2,v2poly)=[(load(prefix['Fn_A'+str(i)+'_pseudo']['u']),load(prefix['Fn_A'+str(i)+'_pseudo']['v']))for i in[1,2]]
  determinant=u1*v2poly-u2*v1;norms=[a*v*v-b*u*v+c*u*u for u,v in[(u1,v1),(u2,v2poly)]]
  # Verify the UNIVERSAL scalar norm identity once; actual coefficient substitution is exact polynomial arithmetic.
  aa,bb,cc,uu,vv,x=S.symbols('aa bb cc uu vv x');ff=aa*x*x+bb*x+cc;ll=uu*x+vv
  assert S.Poly(uu**2*ff-(aa*ll+bb*uu-2*aa*vv)*ll-(aa*vv**2-bb*uu*vv+cc*uu**2),x,aa,bb,cc,uu,vv).is_zero
  cert['Fn_regular_branch']=dict(required_unit='a=LC_s(Fn)!=0 ONLY in regular branch',a=rec(a),b=rec(b),c=rec(c),linear_remainder_determinant=rec(determinant),two_quadratic_norm_conditions=[rec(p)for p in norms],
   common_root_necessary_conditions='determinant=norm1=norm2=0',universal_norm_identity_verified=True,
   norm_identity='u²Fn-(a*(us+v)+b*u-2a*v)*(us+v)=a*v²-b*u*v+c*u²',
   leading_a_zero_case_is_retained=True,no_global_inversion_of_unknown_r_polynomial=True)
  save('Fn_determinant_and_two_norm_conditions_completed')
  ring=S.QQ.poly_ring(r,R);linear_root_hom=[]
  for name in['A1','A2']:
   original=load(data['rational_coefficient_pairs'][name]['numerator']);ps=S.Poly(original.as_expr(),s,domain=ring);degree=ps.degree();value=P(0)
   for(e,),coef in ps.terms():value+=P(coef)*((-c)**e)*(b**(degree-e))
   linear_root_hom.append(value)
  cert['Fn_leading_zero_cases']=dict(a_zero_not_inverted=True,
   case_a0_b_nonzero='a=0,b!=0: b^deg(Ai)*Ai(s=-c/b)=0 necessary for i1,2',homogenized_linear_root_conditions=[rec(p)for p in linear_root_hom],
   case_a0_b0='a=b=0 then c=0 necessary, otherwise Fn has no root; Ai(s)=0 constraints not discarded',
   these_finite_boundary_r_cases_NOT_excluded=True)
  save('Fn_all_leading_zero_cases_explicitly_retained')
  leading=[]
  for name in['Z1','Z2']:
   con=data['localized_constraints'][name];degree=con['shape_bidegree'][1]
   for key in['coefficient1','coefficient2','denominator']:
    p=load(dict(variables=['r','s','R'],**con[key]));leading.append(P(S.Poly(p.as_expr(),s,domain=ring).nth(degree)))
  cert['s_infinity']=dict(s_degrees=[3,5],all6_independent_target_leading_coefficients=[rec(p)for p in leading],
   necessary_condition='all6 vanish for a fixed leading root persisting over independentZ0,Z1,Z2; no gcd/infinity elimination executed')
  cert['Hn_boundary']=dict(observation='Away from J*Fn and known units, both normalized maps extend regularly on properHn=0 curve. Its image cannot dominate the2D independentZ1,Z2 base. Pole/base-locus intersections require the explicitJ/Fn cases.',no_Hn_resultant_or_saturation=True)
  save('boundary_preflight_completed_without_further_gcd_or_elimination')
 except Exception as exc:error=repr(exc);save('remaining_budget_partial_frontier')
 signal.alarm(0)
 report=dict(id=STEM,all_passed=error is None,error=error,
  continuation_only_from_verified_saved_prefix=True,no_repeated_J_gcd_or_Fn_pseudodivision=True,
  original_prefix_seconds=prior['seconds'],remaining_cap_seconds=remaining,combined_mathematical_seconds=prior['seconds']+time.monotonic()-start,
  J_necessary_residual_boundary='25r-1=0 after knownr,T units, genericR-onlyBezoutdenominators retained',
  Fn_regular_norm_and_leading_cases_completed='Fn_leading_zero_cases'in cert,
  no_Q19_general_factorization_fullresultant_PRS_point_prime_or_source26=True,
  no_primitive_Q19_actual_degree_targetindex_or_NEW_equation_claim=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),prefixreport,prefixpath,source]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024-10000
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','error','combined_mathematical_seconds','Fn_regular_norm_and_leading_cases_completed','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
