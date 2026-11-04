#!/usr/bin/env python3
"""Remove the known square-root source ambiguity by an exact Sp4 gauge."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from weight1_structure import J,sp_basis
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s cap')));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False)
 src=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
 inv=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 d=json.loads(src.read_text());iv=json.loads(inv.read_text());weights=[x[0]for x in iv['Cartan_weights']]
 u,v,c,X,W,B,xi,r,t,q=s.symbols('u v c X W B xi r s q');h=t*t-25*r*q*q
 # Finite torus weights in the accepted representation, not a parameter fit.
 Hstd=s.diag(1,-1,0,0);coef=J.inv()*Hstd;Hrep=s.zeros(16);re=s.zeros(4);k=0
 for i in range(4):
  for j in range(i,4):
   re+=coef[i,j]*sp_basis()[k];Hrep+=coef[i,j]*s.Matrix(iv['all_ten_representation_matrices'][k]);k+=1
 assert re==Hstd and Hstd.T*J+J*Hstd==s.zeros(4) and Hrep==s.diag(*weights)
 K=s.QQ.frac_field(xi,r,t,q,B)
 mod=s.Poly(u*u-1600*r,u,domain=K)
 subs={v:-32*u*t/h,c:xi*B*u/40,X:t/(48*u),W:(q+1)/3072}
 result=[];checks=[]
 for i,expr in enumerate(d['primitive_beta16']):
  z=s.cancel(s.sympify(expr).subs(subs,simultaneous=True)*u**weights[i]);num,den=s.fraction(z)
  n=s.Poly(num,u,domain=K).rem(mod);de=s.Poly(den,u,domain=K).rem(mod)
  a,b=n.nth(0),n.nth(1);cc,dd=de.nth(0),de.nth(1)
  norm=cc*cc-1600*r*dd*dd;assert norm!=0
  even=s.cancel((a*cc-1600*r*b*dd)/norm);odd=s.cancel((b*cc-a*dd)/norm);assert odd==0
  result.append(str(even));checks.append(dict(index=i,Cartan_weight=weights[i],quadratic_remainder_odd_part_zero=True))
 assert sum(z!='0'for z in result)==10
 cert=dict(variables=['xi','r','s','q','B'],h='s^2-25*r*q^2',
  auxiliary_extension='u^2=1600*r',source_substitution={str(k):str(v)for k,v in subs.items()},
  exact_orbit_gauge='g_u=diag(u,u^-1,1,1) in Sp4; beta_bar_i=u^(w1_i)*beta_i',
  beta_bar16=result,all16_checks=checks,
  source_domain='xi*r*s*B*(s^2-25*r*q^2)!=0 together with previously stated source chart units',
  known_tau_ambiguity_removed_in_rational_orbit_representative=True,
  generic_one_orbit_per_source_requires_actual_target_index1_gate=True,
  full_hypersurface_rationality_or_generic_stabilizer_trivial_proved=False)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),ROOT/'code/weight1_structure.py',src,inv]
 report=dict(id=ID,all_passed=True,exact_rational16_vector_on_five_known_fixed_source_coordinates=True,
  nonzero_beta_coordinates=10,no_new_point_invariant_expansion_or_source26=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
