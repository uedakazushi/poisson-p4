#!/usr/bin/env python3
"""Independent closed powers/full150 arcs and four MPS period memberships."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from ambient_CA_global_flag_review_20261103_v63 import IDS
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha
from q_mps_arc_membership_review_20261103_v63 import dual,skew

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths=[Path(__file__),public,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/q_mps_arc_membership_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
 witnesses=read(ROOT/'code/data/q_mps_extended_toric_arc_witnesses_20261103_v63.json');arcs=read(folder/'full150-arcs.json.gz');periodcert=read(folder/'period-certificates.json');catalog=read(ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json')['diagrams']
 def wc(A,B):return[A[i][p]*B[j][q]-B[i][p]*A[j][q]+(A[i][q]*B[j][p]-B[i][q]*A[j][p]if p!=q else 0)for((i,j),(p,q))in IDS]
 for witness,arc,certificate in zip(witnesses,arcs,periodcert):
  parameter=s.Symbol(witness['free']);K=s.QQ.frac_field(parameter);R,t=ring('t',K);cv=lambda e:K.from_sympy(s.cancel(e));rho=s.sympify(witness['rho']);alpha=s.sympify(witness['alpha']);beta=s.sympify(witness['beta']);di=list(map(s.sympify,witness['eigenvalues']));ss=3*(rho-1)/2;aa=8*rho/(9*(rho*rho-1));mu=-4/(9*ss*ss*(ss+3));links=[s.S.One,ss,(ss+3)/2,-9/(ss+3)]
  assert cv(s.diff(rho,parameter))and cv(s.diff(aa,parameter))
  assert all(cv(e)for e in[rho,rho-1,rho+1,*[di[i]-di[j]for i in range(5)for j in range(i)]])
  # Closed divided-difference formula for powers of a lower bidiagonal matrix.
  # Recurrence verifies it independently of the producer's matrix products.
  powers=[]
  for k in range(5):
   M=[[R.zero for j in range(5)]for i in range(5)]
   for i in range(5):
    for j in range(i+1):
     dist=i-j
     if dist>k:continue
     homogeneous=sum((s.prod(di[l]for l in mon)for mon in itertools.combinations_with_replacement(range(j,i+1),k-dist)),s.S.Zero)
     M[i][j]=R(cv(s.prod(links[l]for l in range(j,i))*homogeneous))*t**(k-dist)
   powers.append(M)
  for k in range(4):
   assert all(powers[k+1][i][j]==R(cv(di[i]))*t*powers[k][i][j]+(R(cv(links[i-1]))*powers[k][i-1][j]if i else R.zero)for i in range(5)for j in range(5))
  tau=[sum((powers[k][i][i]for i in range(5)),R.zero)/5 for k in range(1,5)];W=[[[powers[k][i][j]-(tau[k-1]if i==j else R.zero)for j in range(5)]for i in range(5)]for k in range(1,5)]
  terms=[(0,1,R(cv(mu*alpha))*t*t),(0,2,R(cv(mu*beta))*t),(0,3,R(cv(mu*rho))),(1,2,R(cv(mu)))];raw=[R.zero]*150;plain=[R.zero]*150;Z=[[R.zero for j in range(5)]for i in range(5)]
  for i,j,coefficient in terms:
   raw=[u+coefficient*v for u,v in zip(raw,wc(W[i],W[j]))];plain=[u+coefficient*v for u,v in zip(plain,wc(powers[i+1],powers[j+1]))]
   for l in range(5):
    for m in range(5):Z[l][m]+=coefficient*(tau[j]*powers[i+1][l][m]-tau[i]*powers[j+1][l][m])
  identity=[[R.one if i==j else R.zero for j in range(5)]for i in range(5)];assert all(u-v-w==R.zero for u,v,w in zip(raw,plain,wc(identity,Z)))
  saved=[R.from_dict({(degree,):cv(s.sympify(c))for degree,c in row})for row in arc['canonical_raw150']];assert raw==saved
  target=[K.zero]*150
  for key,value in{((1,4),(0,0)):aa,((2,3),(0,0)):-s.Rational(2,9),((2,4),(0,1)):s.Rational(4,3),((3,4),(0,2)):s.Rational(2,3),((3,4),(1,1)):1}.items():target[IDS[key]]=cv(value)
  assert all(p.get((0,),K.zero)==v for p,v in zip(raw,target))
  S=s.Matrix([[s.sympify(v)for v in row]for row in arc['eigenframe_coefficients']]);assert all(cv((di[i]-di[j])*S[i,j]+links[i-1]*S[i-1,j])==K.zero for i in range(5)for j in range(i))
  pairs=list(itertools.combinations(range(1,5),2));gamma=s.Matrix([alpha*(di[i]*di[j]**2-di[i]**2*di[j])+beta*(di[i]*di[j]**3-di[i]**3*di[j])+rho*(di[i]*di[j]**4-di[i]**4*di[j])+di[i]**2*di[j]**3-di[i]**3*di[j]**2 for i,j in pairs]);period=dual(gamma)
  rec=next(v for v in catalog if v['id']==witness['component']);basis=s.Matrix(rec['period_basis']).applyfunc(s.Rational).T;coordinates=s.Matrix(list(map(s.sympify,certificate['period_coordinates'])));assert all(cv(v)==K.zero for v in basis*coordinates-period)
  pf=gamma[0]*gamma[5]-gamma[1]*gamma[4]+gamma[2]*gamma[3];V=s.prod(di[i]for i in range(1,5))*s.prod(di[j]-di[i]for i in range(1,5)for j in range(i+1,5));assert cv(pf-rho*V)==K.zero and cv(pf)
  assert cv(s.diff(aa,parameter)-s.sympify(arc['a_derivative']))==K.zero
 a,b,r=s.symbols('alpha beta rho');C=skew([a,b,r,1,0,0]);assert C[0,1]*C[2,3]-C[0,2]*C[1,3]+C[0,3]*C[1,2]==r
 paths += [ROOT/'reports/Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63_REVIEW.json',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',ROOT/'research/notes/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63_REVIEW.md']
 report=dict(id='Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_all_evidence_match=True,independent_closed_bidiagonal_powers_and_full150_all4_actual_arcs=True,
  all_Qa_limits_Euler_identities_and_eigenframe_recurrences_independently_checked=True,
  fresh_all4_physical_divergence25_and_Jacobi350_receipts_retained=True,
  extra_t2_W1_W2_and_t_W1_W3_vanish_and_match_common_t5_toric_coefficients=True,
  all6_actual_period_basis_coordinates_and_nonzero_Pfaffians_verified=True,
  rational_rho_and_a_density_closedness_give_ALL_finite_Qa_memberships=['T02','T03','T16','T17'],
  internal_MPS_T16_T17_distinguished_from_R_boundary_strata=True,
  nonmembership_new_count_or_local_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS all4 full150 extended arcs and ALL-Qa MPS T02,T03,T16,T17 memberships')

if __name__=='__main__':main()
