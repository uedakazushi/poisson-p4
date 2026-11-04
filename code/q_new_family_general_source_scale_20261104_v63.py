#!/usr/bin/env python3
"""Exact source scaling for every homogeneous Sp4 invariant, not a probe."""
import argparse, hashlib, json, resource, signal, time
from pathlib import Path
import sympy as s
from weight1_structure import J,sp_basis
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s cap')));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False)
 source=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
 inv=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 d=json.loads(source.read_text());q=json.loads(inv.read_text());c=s.Symbol('c')
 e=d['exact_c_exponents_beta16'];weights=q['Cartan_weights'];assert len(e)==16
 assertions=[]
 for i,expr in enumerate(d['primitive_beta16']):
  f=s.sympify(expr);g=s.cancel(f/c**e[i]);assert c not in g.free_symbols
  assert 2*e[i]-1==-sum(weights[i])
  assertions.append(dict(index=i,c_power=e[i],torus_weight=2*e[i]-1,coefficient_independent_of_c=True))
 # Same standard symplectic basis and original16-dimensional representation.
 h=s.diag(-1,1,-1,1);coef=J.inv()*h;rho=[s.Matrix(z)for z in q['all_ten_representation_matrices']]
 gs=sp_basis();assert len(gs)==len(rho)==10
 reconstructed=s.zeros(4);H=s.zeros(16);k=0
 for i in range(4):
  for j in range(i,4):
   reconstructed+=coef[i,j]*gs[k];H+=coef[i,j]*rho[k];k+=1
 assert reconstructed==h and h.T*J+J*h==s.zeros(4)
 assert H==s.diag(*[2*x-1 for x in e])
 cert=dict(source_coordinates=d['source_variables'],source_c_scaling_exponents=e,
  standard_cocharacter='diag(lambda^-1,lambda,lambda^-1,lambda)',
  cocharacter_standard_exponents=[-1,1,-1,1],cocharacter_W16_exponents=[2*x-1 for x in e],
  all16_coordinate_checks=assertions,
  exact_global_identity='beta(u,v,lambda^2*c,X,W,B)=lambda*rho(g_lambda)*beta(u,v,c,X,W,B)',
  every_homogeneous_Sp4_invariant_degree_d='I(beta(c))=c^(d/2)*I(beta(1)), d even',
  known_fixed_source_coordinates='xi=40c/(Bu),r=u^2/1600,s=48uX,R,Z0',
  consequence_degree14='I14=xi^7*h(r,s,R,Z0) for some rational h; if I14 nonzero on NEW then h nonzero',
  index_one_criterion_conditional='If spectral ratios and ZN recover the entire shape field, I14/F4^3 recovers xi and hence the known fixed source field',
  ZN_primitivity_or_target_index_one_proved=False)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'code/weight1_structure.py',ROOT/'research/notes'/(ID+'.md'),source,inv]
 report=dict(id=ID,all_passed=True,all16_exact_c_exponents_and_original_torus_representation_verified=True,
  general_all_degrees_covariance_not_low_word_extrapolation=True,
  no_new_source_point_invariant_expansion_resultant_or_source26=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=ap.parse_args();main(a.output.resolve())
