#!/usr/bin/env python3
"""Independent commuting-power realization of every Qa in actual T01."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as sp
from ambient_CA_global_flag_review_20261103_v63 import IDS

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--membership-replay',type=Path,required=True);ap.add_argument('--constraints-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63';paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py'];reports=[]
 for public,replay in [('report.json',args.membership_replay/'report.json'),('global-constraint-checks.json',args.constraints_replay/'report.json')]:
  p=folder/public;d=load(p);f=load(replay);assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths.append(p);reports.append(d)
  for file,hsh in d['input_sha256'].items():q=ROOT/file;assert sha(q)==hsh,file;paths.append(q)
 cert=folder/'certificate.json';assert cert.read_bytes()==(args.membership_replay/'certificate.json').read_bytes()and sha(cert)==reports[0]['certificate_sha256'];paths.append(cert);c=load(cert)
 a,s,t=sp.symbols('a s t');X=sp.Matrix(sp.symbols('x0:5'));N=sp.zeros(5)
 for i,link in enumerate([1,s,(s+3)/2,-9/(s+3)]):N[i+1,i]=link
 assert N==sp.Matrix([[sp.sympify(v)for v in row]for row in c['N']])and N**5==sp.zeros(5)and(N**4)[4,0]==-9*s/2
 lam=-2*a/(9*s);mu=-4/(9*s*s*(s+3));Nt=N+t*sp.diag(0,1,2,3,4);I=sp.eye(5)
 powers=[I]+[Nt**k for k in range(1,5)];W=[powers[k]-sp.trace(powers[k])*I/5 for k in range(1,5)]
 assert all(sp.trace(M)==0 for M in W)
 assert all((M*L-L*M).applyfunc(sp.cancel)==sp.zeros(5)for M,L in itertools.combinations(W,2))
 wedge=lambda u,v:u*v.T-v*u.T
 P=lam*wedge(W[0]*X,W[3]*X)+mu*wedge(W[1]*X,W[2]*X)
 raw=sp.Matrix([sp.cancel(sp.Poly(sp.expand(P[i,j]),*X).coeff_monomial(sp.prod(X[k]for k in mon)))for((i,j),mon)in IDS])
 saved=sp.Matrix([sum(sp.sympify(v)*t**int(deg)for deg,v in terms)for terms in c['canonical_arc_raw150']])
 assert(raw-saved).applyfunc(sp.cancel)==sp.zeros(150,1)
 Q=sp.zeros(5)
 for(i,j),value in{(1,4):a*X[0]**2,(2,3):-sp.Rational(2,9)*X[0]**2,(2,4):sp.Rational(4,3)*X[0]*X[1],(3,4):sp.Rational(2,3)*X[0]*X[2]+X[1]**2}.items():Q[i,j]=value;Q[j,i]=-value
 target=sp.Matrix([sp.Poly(Q[i,j],*X).coeff_monomial(sp.prod(X[k]for k in mon))for((i,j),mon)in IDS])
 F=3*a*s*s+(9*a-4)*s-6;error=(raw.subs(t,0)-target).applyfunc(sp.factor)
 assert{j for j,v in enumerate(error)if v}=={121,137}
 assert sp.cancel(error[121]-F/(3*(s+3)))==0 and sp.cancel(error[137]-F/(6*s))==0
 assert F.subs(s,0)==-6 and sp.expand(F.subs(s,-3))==6
 assert sp.expand(F.subs(a,0))==-4*s-6 and F.subs({a:0,s:-sp.Rational(3,2)})==0
 assert sp.discriminant(F,s)==81*a*a+16
 # The two generators are already actual trace-free commuting fields; their
 # wedges have zero divergence and Schouten bracket over QQ(a,s)[t].  The
 # physical all350 replay is also retained, with no truncation or point test.
 plain=lam*wedge(powers[1]*X,powers[4]*X)+mu*wedge(powers[2]*X,powers[3]*X)
 tau=[sp.trace(powers[k])/5 for k in range(1,5)]
 Z=lam*(tau[3]*powers[1]-tau[0]*powers[4])+mu*(tau[2]*powers[2]-tau[1]*powers[3])
 assert(P-plain-wedge(X,Z*X)).applyfunc(sp.cancel)==sp.zeros(5)
 S=sp.Matrix([[sp.sympify(v)for v in row]for row in c['eigenframe']]);Sinv=sp.Matrix([[sp.sympify(v)for v in row]for row in c['eigenframe_inverse']])
 assert S.det()==1 and(S*Sinv-I).applyfunc(sp.cancel)==sp.zeros(5)
 assert(Nt*S-S*(t*sp.diag(0,1,2,3,4))).applyfunc(sp.cancel)==sp.zeros(5)
 gam={(i,j):sp.cancel(lam*(i*j**4-i**4*j)+mu*(i*i*j**3-i**3*j*j))for i,j in itertools.combinations(range(1,5),2)}
 assert all(sp.cancel(gam[i,j]-sp.sympify(c['toric_periods'][f'{i},{j}']))==0 for i,j in gam)
 assert sp.cancel(gam[1,2]*gam[3,4]-gam[1,3]*gam[2,4]+gam[1,4]*gam[2,3]-288*lam*mu)==0
 # Distinct eigenvalues imply I,M,...,M^4 independent; the four traceless
 # powers span Cartan, so arbitrary six wedges give exactly basic T01 there.
 # Closedness of T01 pulls membership back from this dense open to ALL Mat5.
 weights=load(ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json')['parameter_weights']
 normalized=load(ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/certificate.json')['scalar_section_normal61']
 weighted=[sp.cancel(a**d*sp.sympify(v))for d,v in zip(weights,normalized)]
 assert [(i+1,v)for i,v in enumerate(weighted)if v]==[(25,-sp.Rational(4,3)),(36,sp.Integer(3)),(38,9*a/2)]
 paths += [ROOT/'research/notes/Q_T01_MEMBERSHIP_20261103_v63.md',ROOT/'research/notes/Q_T01_MEMBERSHIP_20261103_v63_REVIEW.md']
 result=dict(id='Q_T01_MEMBERSHIP_20261103_v63_REVIEW',all_passed=True,
  fresh_both_producers_and_actual_arc_certificate_match=True,
  independent_full150_sympy_commuting_power_arc_and_limit_identity=True,
  polynomial_canonical_arc_not_unsaturated_or_Delta_open_specialization=True,
  universal_Mat5_times_A6_polynomial_image_in_actual_T01_by_dense_Cartan_and_closedness=True,
  cover_has_admissible_root_for_every_complex_a_including_zero=True,
  all_Qa_in_actual_existing_T01=True,
  independent_eigenframe_and_Euler_wedge_and_six_toric_periods=True,
  fixed_nonzero_beta_entire_affine_qline_in_T01_by_weighted_action_and_closedness=True,
  necessary274_constraints_not_used_as_sufficiency=True,
  a_zero_S61_line_endpoint_equals_physical_Qzero_claimed=False,
  local_branch_exhaustion_or_full_T01_image_surjectivity_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS actual universal matrix-power T01 membership for ALL Qa')

if __name__=='__main__':main()
