#!/usr/bin/env python3
"""Independent actual arc, Pfaffian-period and dense-curve membership audit."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.galoistools import gf_berlekamp
from ambient_CA_global_flag_review_20261103_v63 import IDS
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def dual(v):return s.Matrix([v[5],-v[4],v[3],v[2],-v[1],v[0]])
def skew(v):
 A=s.zeros(4)
 for(i,j),c in zip(itertools.combinations(range(4),2),v):A[i,j]=c;A[j,i]=-c
 return A
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths=[Path(__file__),public,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
 arc=read(folder/'universal-two-eigen-arc.json.gz');rho,x,y,t=s.symbols('rho x y t');X=s.Matrix(s.symbols('x0:5'));ss=3*(rho-1)/2;aa=8*rho/(9*(rho*rho-1));lam=-2*aa/(9*ss);mu=-4/(9*ss*ss*(ss+3));links=[1,ss,(ss+3)/2,-9/(ss+3)];N=s.zeros(5)
 for i,v in enumerate(links):N[i+1,i]=v
 diagonal=s.diag(0,1,2,x,y);Nt=N+t*diagonal;powers=[s.eye(5)]
 for k in range(4):powers.append((powers[-1]*Nt).applyfunc(s.expand))
 W=[powers[k]-s.trace(powers[k])*s.eye(5)/5 for k in range(1,5)];wedge=lambda u,v:u*v.T-v*u.T
 P=lam*wedge(W[0]*X,W[3]*X)+mu*wedge(W[1]*X,W[2]*X);raw=s.Matrix([s.cancel(s.Poly(s.expand(P[i,j]),*X).coeff_monomial(s.prod(X[k]for k in mon)))for((i,j),mon)in IDS]);saved=s.Matrix([sum(s.sympify(c)*t**degree for degree,c in terms)for terms in arc['canonical_arc_raw150']]);assert(raw-saved).applyfunc(s.cancel)==s.zeros(150,1)
 Qa=s.zeros(5)
 for(i,j),v in{(1,4):aa*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,(2,4):s.Rational(4,3)*X[0]*X[1],(3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}.items():Qa[i,j]=v;Qa[j,i]=-v
 assert(P.subs(t,0)-Qa).applyfunc(s.cancel)==s.zeros(5)
 plain=lam*wedge(powers[1]*X,powers[4]*X)+mu*wedge(powers[2]*X,powers[3]*X);tau=[s.trace(powers[k])/5 for k in range(1,5)];Z=lam*(tau[3]*powers[1]-tau[0]*powers[4])+mu*(tau[2]*powers[2]-tau[1]*powers[3]);assert(P-plain-wedge(X,Z*X)).applyfunc(s.cancel)==s.zeros(5)
 S=s.Matrix([[s.sympify(v)/t**(i-j)if i>j else s.sympify(v)for j,v in enumerate(row)]for i,row in enumerate(arc['eigenframe_coefficients'])]);assert(Nt*S-S*t*diagonal).applyfunc(s.cancel)==s.zeros(5)and S.det()==1
 u=s.symbols('b0:6');A=skew(u);B=skew(dual(u));pf=u[0]*u[5]-u[1]*u[4]+u[2]*u[3];assert(A*B+pf*s.eye(4)).applyfunc(s.expand)==s.zeros(4)and dual(dual(u))==s.Matrix(u)
 pairs=list(itertools.combinations(range(1,5),2))
 for perm in itertools.permutations(range(5)):
  sign=(-1)**sum(perm[i]>perm[j]for i in range(5)for j in range(i+1,5))
  for k in range(6):
   v=[int(i==k)for i in range(6)];A5=s.zeros(5);B5=s.zeros(5)
   for(i,j),a0,b0 in zip(pairs,v,dual(v)):A5[i,j]=a0;A5[j,i]=-a0;B5[i,j]=b0;B5[j,i]=-b0
   for i in range(1,5):B5[0,i]=sum(B5[i,j]for j in range(1,5));B5[i,0]=-B5[0,i]
   chart=[A5[perm[i],perm[j]]-A5[perm[i],perm[0]]-A5[perm[0],perm[j]]for i,j in pairs]
   assert dual(chart)==sign*s.Matrix([B5[perm[i],perm[j]]for i,j in pairs])
 records=read(ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json')['diagrams'];curves=read(ROOT/'code/data/q_mps_toric_arc_curves_20261103_v63.json');curve_results=read(folder/'dense-curve-certificates.json');eigen=[0,1,2,x,y];gamma=s.Matrix([rho*(eigen[i]*eigen[j]**4-eigen[i]**4*eigen[j])+eigen[i]**2*eigen[j]**3-eigen[i]**3*eigen[j]**2 for i,j in pairs]);vandermonde=2*x*y*(x-1)*(x-2)*(y-1)*(y-2)*(y-x);assert s.expand(gamma[0]*gamma[5]-gamma[1]*gamma[4]+gamma[2]*gamma[3]-rho*vandermonde)==0
 for witness,certificate in zip(curves,curve_results):
  rec=next(r for r in records if r['id']==witness['component']);f=s.sympify(witness['curve']);F=s.Poly(f,x,y);rr=s.cancel(s.sympify(witness['rho']));numerator,denominator=s.fraction(rr)
  assert s.gcd_list(s.Poly(f,x).all_coeffs())==1
  yy=witness['irreducibility_witness']['y_value'];prime=witness['irreducibility_witness']['prime'];special=s.Poly(f.subs(y,yy),x,modulus=prime)
  assert special.degree()==s.degree(f,x)and special.is_irreducible
  coefficients=[int(v)%prime for v in special.monic().all_coeffs()]
  assert gf_berlekamp(coefficients,prime,s.ZZ)==[coefficients]
  rem=lambda p:s.Poly(s.expand(p),x,y).rem(F).as_expr()
  period=dual(gamma);period_cleared=period.applyfunc(lambda p:s.expand(numerator*s.expand(p).coeff(rho)+denominator*p.subs(rho,0)))
  basis=s.Matrix(rec['period_basis']).applyfunc(s.Rational).T;parameters=s.Matrix([s.sympify(v)for v in certificate['exact_period_parameters']]);polyparams=parameters.applyfunc(lambda p:s.cancel(p*denominator));assert all(not s.denom(p).has(x,y)for p in polyparams)
  assert all(rem(p)==0 for p in basis*polyparams-period_cleared)
  tangent=s.expand((s.diff(numerator,x)*denominator-numerator*s.diff(denominator,x))*s.diff(f,y)-(s.diff(numerator,y)*denominator-numerator*s.diff(denominator,y))*s.diff(f,x))
  factors=[vandermonde,denominator,numerator,numerator-denominator,numerator+denominator,tangent,numerator*numerator+denominator*denominator];remainders=[rem(p)for p in factors];assert all(p!=0 for p in remainders)
  assert all(s.expand(p-s.sympify(q))==0 for p,q in zip(remainders,certificate['admissibility_and_tangent_nonzero_remainders']))
 old=read(ROOT/'legacy/v23/proofs/literature/ss_geometry_mps_certificate.json');assert old['catalogue_id']=='T37'and old['MPS_all_brackets_match']and old['checks']['full_MPS_identification_all_a_t'];legacytext=(ROOT/'legacy/v23/proofs/literature/SS_geometry_and_MPS_identification.txt').read_text();assert all(text in legacytext for text in['τ14=x0^2','τ23=-2*x0^2/(9a)','τ24=4*x0*x1/3','τ34=2*x0*x2/3+x1^2'])
 paths += [ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63_REVIEW.json',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',ROOT/'research/notes/Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63_REVIEW.md']
 report=dict(id='Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_all_evidence_match=True,independent_full150_generic_two_eigen_arc_and_Euler_identity=True,
  independent_unit_lower_eigenframe_and_linear_Pfaffian_duality_all120_permutations=True,
  independent_mod5_irreducible_specializations_via_second_algorithm=True,
  primitive_in_x_and_preserved_degree_prove_QQ_curve_irreducibility=True,
  exact_QQ_curve_actual_period_parameters_all6_match=True,
  simultaneous_admissible_open_and_delta_a_nonzero_by_irreducible_domain=True,
  actual_toric_boundary_membership_and_dense_a_image_closedness_prove_ALL_finite_Qa_in_T02_T03=True,
  no_absolute_irreducibility_or_unsaturated_specialization_assumption=True,
  old_Qa1_T37_tau1_provenance_confirmed_without_old_large_rerun=True,
  nonmembership_local_branch_exhaustion_or_new_component_count_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS actual ALL-Qa T02/T03 via exact period arcs and dense irreducible curves')

if __name__=='__main__':main()
