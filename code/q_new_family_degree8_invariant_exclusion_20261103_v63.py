#!/usr/bin/env python3
"""Exact degree-eight invariant basis and injective restriction to the new beta image."""
import argparse
from collections import defaultdict
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def matrix(rows):return s.Matrix([[s.Rational(c)for c in row]for row in rows])
def qrows(M):return [[str(c)for c in row]for row in M.tolist()]

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
 out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 data=json.loads(invpath.read_text());B=matrix(data['primitive24_by16_basis']);weights=list(map(tuple,data['Cartan_weights']));omega=matrix(data['symplectic_form16']);rep=[matrix(M)for M in data['all_ten_representation_matrices']]
 rows=list(B.T.rref()[1]);left=B[rows,:].inv();generators=sp_basis();flat=s.Matrix.hstack(*[s.Matrix(list(G))for G in generators]);flatrows=list(flat.T.rref()[1]);flatleft=flat[flatrows,:].inv()
 assert omega.T==-omega and omega.det()!=0
 assert all(G.T*omega+omega*G==s.zeros(16)for G in rep)
 # These are pullback derivatives, hence an anti-representation.
 for i,A in enumerate(generators):
  for j,C in enumerate(generators):
   comm=A*C-C*A;coordinates=flatleft*s.Matrix(list(comm))[flatrows,:]
   assert flat*coordinates==s.Matrix(list(comm))
   assert rep[i]*rep[j]-rep[j]*rep[i]==-sum((c*G for c,G in zip(coordinates,rep)),s.zeros(16))
 gram=s.Matrix([[s.trace(A*C)for C in generators]for A in generators]);graminverse=gram.inv();omegainverse=omega.inv()
 # Complete characters Sym^d(W16*), counting all monomials with multiplicity.
 degree_max=10;characters=[defaultdict(int)for _ in range(degree_max+1)];characters[0][0,0]=1
 for a,b in weights:
  for d in range(1,degree_max+1):
   for(i,j),count in list(characters[d-1].items()):characters[d][i+a,j+b]+=count
 rho=(2,1);weyl=[]
 for perm in[(0,1),(1,0)]:
  for signs in itertools.product([-1,1],repeat=2):
   wr=tuple(signs[i]*rho[perm[i]]for i in range(2));det=signs[0]*signs[1]*(-1 if perm==(1,0)else 1)
   weyl.append((tuple(rho[i]-wr[i]for i in range(2)),det))
 multiplicities=[]
 for d,char in enumerate(characters):
  assert sum(char.values())==s.binomial(15+d,d)
  # Independent invariance check on every weight multiplicity.
  assert all(char[(a,b)]==char[(-a,b)]==char[(a,-b)]==char[(b,a)]for a,b in list(char))
  terms=[dict(weight=list(delta),sign=sign,weight_multiplicity=char[delta])for delta,sign in weyl]
  value=sum(row['sign']*row['weight_multiplicity']for row in terms)
  assert value=={0:1,4:1,8:5}.get(d,0)
  multiplicities.append(dict(degree=d,total_monomials=int(s.binomial(15+d,d)),weight_zero=char[0,0],Weyl_terms=terms,invariant_dimension=value))
 words=[word for word in itertools.product(range(2),repeat=4)if word==min(word[i:]+word[:i]for i in range(4))]
 def evaluate_invariants(b):
  mu=s.Matrix([(b.T*omega*G*b)[0]/2 for G in rep]);coeff=graminverse*mu
  moment=sum((c*A for c,A in zip(coeff,generators)),s.zeros(4))
  theta=sum((c*G for c,G in zip(coeff,rep)),s.zeros(16));V=s.Matrix.hstack(*[omega*G*b for G in rep])
  H=theta+omegainverse*V*graminverse*V.T
  matrices=[theta,H];values=[]
  for word in words:
   M=s.eye(16)
   for i in word:M=M*matrices[i]
   values.append(s.factor(s.trace(M)))
  t2=s.factor((mu.T*graminverse*mu)[0]);assert t2==s.factor(s.trace(moment*moment))
  values.append(s.factor(t2*t2));return values,t2,s.factor(moment.det())
 model=S61();z,u,beta=primitive_beta(model);one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
 beta24=[beta[i,j].coeff(zz)for i,j in itertools.combinations(range(4),2)for zz in z]
 betaU=s.Matrix([[s.expand(c).coeff(uu)for uu in u]for c in beta24])
 mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz';payload=json.loads(gzip.decompress(mappath.read_bytes()));dag=DAG(payload['input_names']);dag.nodes=payload['nodes']
 points=[((1,1,1),(0,0,0,0)),((1,1,2),(0,0,0,0)),((2,1,1),(0,0,0,0)),((1,2,3),(0,0,0,0)),((1,1,1),(1,0,0,0))]
 values=[];samples=[]
 for n,base in points:
  point=dict(zip(['n0','n1','n2'],n));point.update(scale=1);point.update(dict(zip(['base1','base2','base3','base4'],base)));point.update({f'g{i}':v for i,v in enumerate([1,0,0,1,0,0,0,0,0,0])})
  allvalues=dag.evaluate(point);U=s.Matrix([s.Rational(allvalues[j])for j in payload['output_u1_through_u61']]);b24=betaU*U[one,:];b=left*b24[rows,:];assert B*b==b24
  row,t2,det=evaluate_invariants(b);values.append(row)
  pivot=allvalues[payload['normalization']['value_pivot_node']];pf=allvalues[payload['normalization']['value_pfaffian_node']];assert pivot and pf and t2
  samples.append(dict(point=point,value_pivot=str(pivot),value_pfaffian=str(pf),primitive_beta24=list(map(str,b24)),invariant_values=list(map(str,row)),moment_trace2=str(t2),moment_determinant=str(det),moment_det_over_trace2square=str(s.factor(det/t2**2))))
 evaluation=s.Matrix(values);_,columns=DomainMatrix.from_Matrix(evaluation).rref();assert len(columns)==5
 determinant=DomainMatrix.from_Matrix(evaluation[:,list(columns)]).convert_to(s.QQ).det();assert determinant
 assert samples[0]['moment_det_over_trace2square']!=samples[1]['moment_det_over_trace2square']
 certificate=dict(
  invariant_dimensions_through10=multiplicities,
  complete_character_recipe='Product over all16 weights of (1-t*x^weight)^(-1), retaining degrees0..10; trivial multiplicity=sum_{w in W(C2)}det(w)*char_d(rho-wrho), rho=(2,1).',
  positive_roots=[[2,0],[0,2],[1,1],[1,-1]],
  degree8_invariant_definitions=dict(omega='pinned nondegenerate invariant alternatingform',mu_A='b^T*omega*rho(A)*b/2',trace_gram='Tr(A_i*A_j)',F='Tr(mu_matrix^2)/2',theta='rho(mu_matrix)',H='omega^-1*Hessian(F)=theta+omega^-1*V*trace_gram^-1*V^T',V='columns omega*rho(A_i)*b',trace_words=words,last_invariant='Tr(mu_matrix^2)^2'),
  universal_invariance_checks=dict(all10_symplectic=True,all100_pullback_anti_Lie_brackets=True,trace_gram_nondegenerate=True),
  actual_new_image_samples=samples,full_evaluation_matrix=qrows(evaluation),independent_five_columns=list(columns),exact5_by5_determinant=str(determinant),
  full_degree8_invariant_dimension=5,restriction_to_new_beta_image_injective=True,
  defining_hypersurface_degree_atleast12=True,
  no_equation_det_moment_minus_constant_times_quartic_squared=True,
  actual_hypersurface_equation_not_constructed=True)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 inputs=[Path(__file__),invpath,mappath,ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/report.json',ROOT/'code/weight1_structure.py',ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'research/notes/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63.md']
 report=dict(id='Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63',all_passed=True,exact_character_dimension8=5,exact_character_dimension10=0,explicit_degree8_invariants_independent_on_actual_new_beta_image=5,all_degree8_hypersurface_equations_excluded=True,defining_homogeneous_invariant_degree_lower_bound=12,actual_equation_or_exhaustion_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FULL DEGREE8 INVARIANT RESTRICTION INJECTIVE; DEFINING DEGREE>=12 PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
