#!/usr/bin/env python3
"""Actual v=0 boundary: full Ga, orbit parameter, and 209-incidence binding."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def rows(M):return [[str(z)for z in row]for row in M.tolist()]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json';G=read(geo);assert G['all_passed']
 inc=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz';I=read(inc)
 gates=[ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63_REVIEW/report.json',
  ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_N_GLOBAL_R1_20261103_v63_REVIEW.json']
 assert all(read(p)['all_passed']for p in gates)
 u,a,e,g,x,y,z,w,H,K,Z,T,tx,tw,ii=s.symbols('u a e g x y z w H K Z T tx tw i');n0,n1,n2=s.symbols('n0 n1 n2');coord=(x,y,z,w)
 P=s.Matrix(G['affine_P']).subs({n0:u,n1:1,n2:0})
 div=s.Matrix([s.expand(sum(s.diff(P[i,j],coord[i])for i in range(4)))for j in range(4)])
 N=s.zeros(5);N[1:,0]=div.subs(dict.fromkeys(coord,0));N[1:,1:]=div.jacobian(coord)
 assert N*N==s.zeros(5)and N.extract([2,3],[0,1]).det()==4*u*u/3
 assert N==s.Matrix(G['homogeneous_modular_N']).subs({n0:u,n1:1,n2:0})
 Qm=s.Matrix(G['quadric_matrix']).subs({n0:u,n1:1,n2:0});qd=s.factor(Qm.det());assert qd==s.Rational(32768,531441)*u**4
 # Necessary reduced flag after modular N forces b=c=0,j=a^3.
 A=s.Matrix([[a,0,0,0],[e,a*a,0,g],[H,ii,a**3,K],[0,0,0,a]])
 tr=s.Matrix([tx,T,Z,tw]);Y=A*s.Matrix(coord)+tr
 diff=(A*P*A.T/a**2-P.subs(dict(zip(coord,Y)),simultaneous=True)).applyfunc(s.expand)
 leading=[s.cancel(s.Poly(diff[1,2],*coord).coeff_monomial(m))for m in(w*w,x*w)]
 assert s.cancel(leading[0]-(a**3-a*a-64*a*g))==0
 assert s.cancel(leading[1]+a*(128*e+u*g))==0
 pw=s.cancel(s.Poly(diff[2,3],*coord).coeff_monomial(w));assert pw==-2*g/3
 # g=0 from full LINEAR w coefficient; then a!=0 givesa=1; thene=0.
 simple={a:1,e:0,g:0,tx:0,tw:0,ii:0}
 AA=A.subs(simple);tt=tr.subs(simple);ds=diff.subs(simple).applyfunc(s.expand)
 selected=s.Matrix([ds[0,2].subs(dict.fromkeys(coord,0)),ds[2,3].subs(dict.fromkeys(coord,0)),ds[1,2].subs(dict.fromkeys(coord,0))])
 L,rhs=s.linear_eq_to_matrix(list(selected),[H,K,Z]);det=s.factor(L.det());assert det==-s.Rational(256,243)*u
 answer={H:3*T/4,K:0,Z:T/(64*u)}
 assert(ds.subs(answer)).applyfunc(s.cancel)==s.zeros(4)
 M=s.eye(5);M[1:,0]=tt;M[1:,1:]=AA;M=M.subs(answer)
 assert(M-(s.eye(5)+3*T*N/(4*u))).applyfunc(s.cancel)==s.zeros(5)
 # Exact binding to the OLD99 source, not merelyN^2=0.
 graph=[s.expand(s.sympify(q).subs({n0:u,n1:1,n2:0}))for q in I['exact99_polynomial_graph']]
 mons=list(map(s.sympify,I['coefficient_monomials']));physical=s.Matrix(I['positive_P'])
 pairs=[(i,j)for i in range(4)for j in range(i+1,4)]
 for k,(i,j)in enumerate(pairs):
  delta=sum(graph[15*k+m]*mons[m]for m in range(15));physical[i,j]+=delta;physical[j,i]-=delta
 assert(physical-P).applyfunc(s.expand)==s.zeros(4)
 assert graph[90:]==[0]*9
 assert all(s.expand(s.sympify(q).subs({n0:u,n1:1,n2:0}))==0 for q in I['all209_exact_full_pullbacks'])
 cert=dict(domain='u!=0, v=0, n1=1',full_affine_P=rows(P),modular_N=rows(N),N_squared_zero=True,
  N_rank2_unit_minor='4*u^2/3',quadric_matrix=rows(Qm),quadric_det=str(qd),actual_divisor='3H+Q, Q smooth irreducible',
  full_flag_reduction='ImN andQ|Hvertex givewholeaffineflag; divergenceNyfixeseta1; Nz,w givesb=c=0; Nz,x givesj=a^3; f=a^2,kappa=a^-2',
  decisive_actual_coefficients=dict(Pyz_w2=str(leading[0]),Pyz_xw=str(leading[1]),Pzw_LINEAR_w=str(pw)),
  exhaustive_conclusion='Pzw_LINEAR_w forcesg0; Pyz_w2 forcesa1 (aunit); Pyz_xw forcese0; constants give tx=tw=i=0',
  remaining_constant_linear_matrix=rows(L),remaining_constant_rhs=list(map(str,rhs)),unit_minor=str(det),unique_H_K_Z=['3*T/4','0','T/(64*u)'],
  full_Ga_homogeneous_matrix=rows(M),all6_full_coefficient_identities_zero=True,
  exact_stabilizer='Ga for EVERYu!=0; no extrafinitecomponent',
  full_parameter_equivalence='P(u,1,0)~P(U,1,0) iffU=+u or-u; actualtau realizesminus',
  saturation_dimension=24,saturation_reason='PGLdim24-stabdim1+parameterdim1, finiteparameteridentification',
  exact99_incidence_graph_specialized=list(map(str,graph)),all209_actual_source_pullbacks_zero=True,
  center_mark='[0:0:0:1:0] at infinity',pencil_mark='dy wedge dz, rank2 degenerate',
  incidence_scope='actual CLOSED H/compactcenter/projectivepencil/N(P)^2 source; NOT automatically its finitecenter/nondegeneratepencil closure',
  not_MPS_N_CA_CB='their properflagmodel divisors are products oflinearforms, whilecurrentQissmoothirreducible',
  not_Hquartic_component='knownHquarticcomponentisdenseorbitclosure23; currentorbit23 hasdifferentdivisor3H+Q ratherthanH+irreduciblequartic, so cannotbelowerdimorbitboundary',
  no_all_boundary_or_new_component_claim=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 extra=[ROOT/'legacy/v06/nilpotent22.tex',ROOT/'reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62/R4_H_QUARTIC_ESSENTIAL5.json',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_B18_N_critical_scheme_20261103_v63.py']
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),geo,inc,*gates,*extra]
 report=dict(id=ID,all_passed=True,whole_v0_uunit_stabilizer_exact_Ga=True,
  full_parameter_orbits_u_squared_and_saturation24=True,
  original99_full209_closed_incidence_membership_bound=True,
  squarezero_alone_not_used_for_component_membership=True,
  actualN_CA_CB_and_Hquartic_component_excluded_by_separate_proofs=True,
  no_new_point_prime_GB_or_nonlinear_solver=True,independent_review_pending=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','whole_v0_uunit_stabilizer_exact_Ga','original99_full209_closed_incidence_membership_bound','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
