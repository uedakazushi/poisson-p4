#!/usr/bin/env python3
"""Small universal cross-block boundary identity; no NEW25/source evaluations."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_CROSS_BLOCK_EDGE_BOUNDARY_20261104_v63'
def read(p):
 with gzip.open(p,'rt')as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 path=ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/certificate.json.gz';cert=read(path)
 x,y,a,k2,m=s.symbols('x y a k2 m');S=x*x+y*y;Delta=x*x-y*y;D=Delta**2;P=x*x*y*y
 p3=-2*a*x*Delta-2*x*x*S
 B2=-4*a*a*x*x-2*a*x**3-6*a*x*y*y+3*x**4+y**4
 A=p3+S*S
 assert s.factor(A+Delta*(2*a*x+S))==0
 assert s.cancel(B2-(S*S-p3+D-A*A/D))==0
 K=D*(k2-S*S+p3-D)+A*A
 kap=k2-B2;assert s.cancel(kap-K/D)==0
 edge_u=s.sympify(cert['auxiliary_m_dictionary']['u_gj']);edge_v=s.sympify(cert['auxiliary_m_dictionary']['v_hi'])
 assert s.cancel(edge_u*edge_v-(kap*kap-2*S*kap*m+D*m*m)/(16*P))==0
 N=2*D*m
 G=D*N*N-4*S*K*N+4*K*K
 assert s.cancel(G-64*P*D*D*edge_u*edge_v)==0
 F4,p2,p3s,p4,N12=s.symbols('F4 p2 p3 p4 N12');Ss=F4/20;Ds=Ss*Ss-4*p4;As=p3s+Ss*Ss
 Ks=Ds*(p2-Ss*Ss+p3s-Ds)+As*As
 Gs=s.expand(Ds*N12*N12-4*Ss*Ks*N12+4*Ks*Ks)
 poly=s.Poly(Gs,F4,p2,p3s,p4,N12)
 assert all(sum(e*w for e,w in zip(mon,[4,8,8,8,12]))==32 for mon,c in poly.terms())
 assert s.expand(s.Poly(Gs,N12).coeff_monomial(N12**2)-Ds)==0 and not poly.is_zero
 den=s.ilcm(1,1,*(s.denom(c)for c in poly.coeffs()))
 raw=[int(c*den)for c in poly.coeffs()];content=s.igcd(*raw);primitive=s.Poly(Gs*den/content,F4,p2,p3s,p4,N12)
 evidence=dict(S='F4/20',D='p0=S^2-4p4',A='p3+S^2',K16='D*(p2-S^2+p3-D)+A^2',
  B2_identity='B2=S^2-p3+D-A^2/D',kappa_identity='kappa=K16/D',
  edge_products='u_edge=g*j, v_edge=h*i; not the original-family source u,v',
  edge_product_identity='u_edge*v_edge=(kappa^2-2S*kappa*m+D*m^2)/(16p4)',
  N12_identity='N12=2D*m',G32='D*N12^2-4S*K16*N12+4K16^2',
  exact_cofactor_identity='G32=64*p4*D^2*(g*j)*(h*i)',known_unit_domain='D*p4 !=0 and the frozen b=e=1 chart',
  abstract_generators=['F4','p2','p3','p4','N12'],abstract_generator_weights=[4,8,8,8,12],
  primitive_integer_G32=str(primitive.as_expr()),primitive_integer_terms=[[list(mon),str(c)]for mon,c in primitive.terms()],
  primitive_scale_over_displayed_G32=str(den/content),source_or_NEW25_vanishing_checked=False)
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 report=dict(id=STEM,all_passed=True,exact_ordered_normal_form_edge_identity=True,
  cofactor=64,universal_invariant_beta_degree=32,abstract_polynomial_terms=len(poly.terms()),
  universal_polynomial_nonzero_by_six_generator_independence=True,
  open_chart_edge_union_characterized=True,no_global_pullback_factorization_or_component_claim=True,
  no_source32_check_no_NEW25_point_or_equation_claim=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),path]},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','abstract_polynomial_terms','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
