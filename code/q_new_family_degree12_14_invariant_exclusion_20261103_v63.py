#!/usr/bin/env python3
"""Complete degree12/14 invariant exclusion for the new first-jet hypersurface."""
import argparse
from collections import defaultdict
import gzip
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from q_new_family_degree8_invariant_exclusion_20261103_v63 import matrix,sha,qrows,ROOT
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 previous=ROOT/'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/report.json';old=json.loads(previous.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(previous.parent/rel)==want,rel
 prevcert=json.loads((previous.parent/'certificate.json').read_text());basis8=prevcert['independent_five_columns'];assert len(basis8)==5
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';data=json.loads(invpath.read_text())
 B=matrix(data['primitive24_by16_basis']);rs=list(B.T.rref()[1]);left=B[rs,:].inv();weights=list(map(tuple,data['Cartan_weights']));O=matrix(data['symplectic_form16']);oi=O.inv();rep=[matrix(M)for M in data['all_ten_representation_matrices']];gs=sp_basis();gram=s.Matrix([[s.trace(A*C)for C in gs]for A in gs]);gi=gram.inv()
 # A separate complete character pass through14, all integer operations.
 characters=[defaultdict(int)for _ in range(15)];characters[0][0,0]=1
 for a,b in weights:
  for d in range(1,15):
   for(i,j),c in list(characters[d-1].items()):characters[d][i+a,j+b]+=c
 rho=(2,1);weyl=[]
 for perm in[(0,1),(1,0)]:
  for signs in itertools.product([-1,1],repeat=2):
   wr=tuple(signs[i]*rho[perm[i]]for i in range(2));weyl.append((tuple(rho[i]-wr[i]for i in range(2)),signs[0]*signs[1]*(-1 if perm==(1,0)else 1)))
 charcert=[]
 for d,char in enumerate(characters):
  assert sum(char.values())==s.binomial(15+d,d)
  terms=[dict(weight=list(delta),sign=sign,multiplicity=char[delta])for delta,sign in weyl];dimension=sum(v['sign']*v['multiplicity']for v in terms)
  assert dimension=={0:1,4:1,8:5,12:9,14:3}.get(d,0)
  charcert.append(dict(degree=d,total_monomials=int(s.binomial(15+d,d)),Weyl_terms=terms,invariant_dimension=dimension))
 def words(d):return [w for w in itertools.product(range(2),repeat=d)if w==min(w[i:]+w[:i]for i in range(d))]
 w4=words(4);w6=words(6);bracketpairs=list(itertools.combinations(range(5),2));DV=[s.Matrix.hstack(*[O*G[:,k]for G in rep])for k in range(16)]
 def matrices(b):
  mu=s.Matrix([(b.T*O*G*b)[0]/2 for G in rep]);co=gi*mu;theta=sum((c*G for c,G in zip(co,rep)),s.zeros(16));V=s.Matrix.hstack(*[O*G*b for G in rep]);U=V*gi;H=theta+oi*U*V.T
  return mu,co,theta,V,U,H
 def traces(b,length):
  mu,co,theta,V,U,H=matrices(b);mats=[theta,H];result=[]
  for word in words(length):
   M=s.eye(16)
   for j in word:M=M*mats[j]
   result.append(s.factor(s.trace(M)))
  if length==4:result.append(s.factor((mu.T*co)[0]**2))
  return result
 def gradients(b):
  mu,co,theta,V,U,H=matrices(b);mats=[theta,H];dc=gi*V.T;derivatives=[]
  for k in range(16):
   dT=sum((dc[j,k]*G for j,G in enumerate(rep)),s.zeros(16));dH=dT+oi*(DV[k]*U.T+U*DV[k].T);derivatives.append([dT,dH])
  result=[]
  for index in basis8:
   if index==6:result.append(4*(mu.T*co)[0]*V*co);continue
   word=w4[index];K=[s.zeros(16),s.zeros(16)]
   for r,letter in enumerate(word):
    M=s.eye(16)
    for j in word[r+1:]+word[:r]:M=M*mats[j]
    K[letter]+=M
   result.append(s.Matrix([sum(K[l][i,j]*derivatives[k][l][j,i]for l in range(2)for i in range(16)for j in range(16))for k in range(16)]))
  return result
 model=S61();z,u,beta=primitive_beta(model);one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
 beta24=[beta[i,j].coeff(zz)for i,j in itertools.combinations(range(4),2)for zz in z];betaU=s.Matrix([[s.expand(c).coeff(uu)for uu in u]for c in beta24])
 mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz';payload=json.loads(gzip.decompress(mappath.read_bytes()));dag=DAG(payload['input_names']);dag.nodes=payload['nodes']
 points=[((1,1,1),(0,0,0,0)),((1,1,2),(0,0,0,0)),((2,1,1),(0,0,0,0)),((1,2,3),(0,0,0,0)),((1,1,1),(1,0,0,0)),((1,1,1),(0,1,0,0)),((2,3,5),(1,2,0,1)),((3,2,1),(0,0,1,0)),((1,3,2),(0,0,0,1))]
 values12=[];values14=[];samples=[];gradient_saved=[];first=None
 for n,base in points:
  point=dict(zip(['n0','n1','n2'],n));point.update(scale=1);point.update(dict(zip(['base1','base2','base3','base4'],base)));point.update({f'g{i}':v for i,v in enumerate([1,0,0,1,0,0,0,0,0,0])})
  v=dag.evaluate(point);U=s.Matrix([s.Rational(v[j])for j in payload['output_u1_through_u61']]);b24=betaU*U[one,:];b=left*b24[rs,:];assert B*b==b24
  pivot=v[payload['normalization']['value_pivot_node']];pf=v[payload['normalization']['value_pfaffian_node']];assert pivot and pf
  row12=traces(b,6);values12.append(row12)
  grads=gradients(b);row14=[s.factor((grads[i].T*oi*grads[j])[0])for i,j in bracketpairs];values14.append(row14)
  if first is None:first=b;firstgrads=grads
  gradient_saved.append([list(map(str,g))for g in grads])
  samples.append(dict(point=point,value_pivot=str(pivot),value_pfaffian=str(pf),primitive_beta24=list(map(str,b24))))
 # Independent exact univariate interpolation checks all five degree8 gradients.
 direction=s.Matrix([1,-1,2,0,-2,1,0,1,1,0,-1,2,0,1,-1,0]);nodes=list(range(-4,5));diffweights=s.finite_diff_weights(1,nodes,0)[1][-1];evaluated=[traces(first+t*direction,4)for t in nodes]
 interpolation=[]
 for j,index in enumerate(basis8):
  numerical=s.factor(sum(c*row[index]for c,row in zip(diffweights,evaluated)));analytic=s.factor((firstgrads[j].T*direction)[0]);assert numerical==analytic
  interpolation.append(str(analytic))
 E12=s.Matrix(values12);E14=s.Matrix(values14);columns12=list(DomainMatrix.from_Matrix(E12).rref()[1]);columns14=list(DomainMatrix.from_Matrix(E14).rref()[1]);assert len(columns12)==9 and len(columns14)==3
 det12=DomainMatrix.from_Matrix(E12[:,columns12]).convert_to(s.QQ).det();det14=DomainMatrix.from_Matrix(E14[:3,columns14]).convert_to(s.QQ).det();assert det12 and det14
 certificate=dict(complete_character_dimensions_through14=charcert,degree12_trace_words=w6,degree14_bracket_degree8_basis_columns=basis8,degree14_bracket_pairs=bracketpairs,
  constant_Poisson_bracket='gradient(f)^T*omega^-1*gradient(g)',all_actual_new_family_samples=samples,degree12_evaluation=qrows(E12),degree14_evaluation=qrows(E14),degree8_gradients=gradient_saved,
  degree12_basis_columns=columns12,degree12_exact9minor=str(det12),degree14_basis_columns=columns14,degree14_first3_exact3minor=str(det14),
  gradient_formula='dtheta=sum_j(Gram^-1*V^T)[j,k]*rho(A_j);dH=dtheta+omega^-1*(dV_k*U^T+U*dV_k^T), U=V*Gram^-1; cyclic trace product rule',
  independent_polynomial_directional_derivative_check=dict(direction=list(direction),nodes=nodes,weights=list(map(str,diffweights)),values=interpolation,degree_bound=8),
  full_degree12_dimension=9,full_degree14_dimension=3,restrictions_injective_at_both_degrees=True,defining_equation_degree_atleast16=True,actual_defining_equation_not_constructed=True)
 # SymPy Integers from list(direction) must be normalized for JSON.
 certificate['independent_polynomial_directional_derivative_check']['direction']=list(map(int,direction))
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 inputs=[Path(__file__),previous,previous.parent/'certificate.json',invpath,mappath,ROOT/'code/weight1_structure.py',ROOT/'research/notes/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63.md']
 report=dict(id='Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63',all_passed=True,full_invariant_dimensions={'12':9,'14':3},actual_new_image_restriction_ranks={'12':9,'14':3},analytic_gradients_independently_interpolated=True,all_degree12_and14_hypersurface_equations_excluded=True,defining_homogeneous_invariant_degree_lower_bound=16,actual_equation_or_full_exhaustion_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FULL DEGREE12/14 INVARIANT RESTRICTION INJECTIVE; DEGREE>=16 PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
