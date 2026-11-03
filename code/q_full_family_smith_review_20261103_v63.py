#!/usr/bin/env python3
"""Independent physical126 tangent/orbit rebuild and all-a Smith audit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,action

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def matrix(rows):return s.Matrix([[s.sympify(v)for v in row]for row in rows])

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63';report=folder/'report.json';archive=folder/'certificate.json.gz'
 d=load(report);fresh=load(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 assert archive.read_bytes()==(args.replay/'certificate.json.gz').read_bytes()and sha(archive)==d['certificate_sha256']
 paths=[Path(__file__),report,archive,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,hsh in d['input_sha256'].items():p=ROOT/file;assert sha(p)==hsh,file;paths.append(p)
 with gzip.open(archive,'rt')as f:c=json.load(f)
 a=s.Symbol('a');X=s.symbols('x0:5');P=s.zeros(5)
 for(i,j),v in{(1,4):a*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,(2,4):s.Rational(4,3)*X[0]*X[1],(3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}.items():P[i,j]=v;P[j,i]=-v
 q=s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in mon))for((i,j),mon)in IDS])
 assert q==s.Matrix(list(map(s.sympify,c['raw_Qa150'])))
 C=matrix(c['canonical_basis']);free=c['canonical_free_raw_indices'];assert C[free,:]==s.eye(126)
 div=s.zeros(25,150)
 for((i,j),mon),col in IDS.items():
  for k in set(mon):
   other=list(mon);other.remove(k);v=mon.count(k)
   if k==i:div[5*j+other[0],col]+=v
   if k==j:div[5*i+other[0],col]-=v
 assert div.rank()==24 and div*C==s.zeros(25,126)and div*q==s.zeros(25,1)
 triples=list(itertools.combinations(range(5),3));mons=list(itertools.combinations_with_replacement(range(5),3))
 monindex={tuple(sorted(m)):i for i,m in enumerate(mons)}
 def rows(R):
  out=s.zeros(350,1)
  for row,(i,j,k)in enumerate(triples):
   f=s.expand(sum(P[i,t]*s.diff(R[j,k],X[t])+R[i,t]*s.diff(P[j,k],X[t])+P[j,t]*s.diff(R[k,i],X[t])+R[j,t]*s.diff(P[k,i],X[t])+P[k,t]*s.diff(R[i,j],X[t])+R[k,t]*s.diff(P[i,j],X[t])for t in range(5)))
   for exp,v in s.Poly(f,*X).terms():
    if v:
     mon=tuple(n for n in range(5)for _ in range(exp[n]));assert len(mon)==3
     out[35*row+monindex[mon]]=v
  return out
 J=s.zeros(350,126)
 for col in range(126):
  R=s.zeros(5)
  for((i,j),mon),index in IDS.items():
   if C[index,col]:R[i,j]+=C[index,col]*s.prod(X[k]for k in mon)
  R=R-R.T;J[:,col]=rows(R)
 assert J==matrix(c['full350_by126_Jacobi'])and rows(P)==s.zeros(350,1)
 positions=list(itertools.product(range(5),repeat=2));O=s.Matrix.hstack(*[action(i,j)*q for i,j in positions])
 assert O==matrix(c['full150_by25_orbit'])and J*O[free,:]==s.zeros(350,25)
 assert sum(((i-2)*O[:,positions.index((i,i))]for i in range(5)),s.zeros(150,1))==5*q
 weights=[i+j-sum(mon)for((i,j),mon)in IDS];degree=[5-weights[i]for i in free]
 target=[10-(sum(tri)-sum(mon))for tri in triples for mon in mons]
 assert degree==c['canonical_relative_weights']and target==c['Jacobi_row_relative_weights']
 assert all(J[i,j]==0 for i in range(350)for j in range(126)if target[i]!=degree[j])
 ranks={};nonunits=[];totals={}
 for key,source in [('Jacobi_blocks',J),('orbit_blocks',O),('initial_nonpositive_group_blocks',O)]:
  total=0
  for block in c[key]:
   B=source[block['rows'],block['columns']];U=matrix(block['left']);V=matrix(block['right']);saved=matrix(block['matrix'])
   assert B==saved
   result=(U*B*V).applyfunc(s.expand)
   diagonal=list(map(s.sympify,block['diagonal']))
   expected=s.zeros(B.rows,B.cols)
   for i,v in enumerate(diagonal):expected[i,i]=v
   assert result==expected
   du=s.factor(U.det());dv=s.factor(V.det())
   assert du==s.sympify(block['left_determinant'])and dv==s.sympify(block['right_determinant'])
   assert du!=0 and dv!=0 and not du.has(a)and not dv.has(a)
   rank=sum(v!=0 for v in diagonal);assert rank==block['rank'];total+=rank
   if key=='Jacobi_blocks':
    nonunits.extend((degree[block['columns'][0]],s.monic(v,a))for v in diagonal if v.has(a))
   else:assert all(not v.has(a)for v in diagonal)
  totals[key]=total
 assert totals=={'Jacobi_blocks':75,'orbit_blocks':17,'initial_nonpositive_group_blocks':7}
 assert nonunits==[(6,a-s.Rational(2,9)),(7,a),(8,a)]
 # These evaluations check stored special receipts; the polynomial Smith
 # equalities above, with constant determinants, prove the ALL-a assertion.
 for value,record in d['exact_ranks_at_special_parameters'].items():
  v=s.Rational(value);rank=sum(s.sympify(e).subs(a,v)!=0 for b in c['Jacobi_blocks']for e in b['diagonal'])
  expected=73 if v==0 else 74 if v==s.Rational(2,9) else 75
  assert rank==expected==record['Jacobi_rank']and record['affine_orbit_dimension']==17
  assert record['affine_tangent_dimension']==126-rank and record['projective_tangent_dimension']==125-rank
  assert record['transverse_tangent_dimension']==109-rank and record['projective_orbit_dimension']==16
  ranks[value]=rank
 paths += [ROOT/'research/notes/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63.md',ROOT/'research/notes/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW.md']
 result=dict(id='Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_compressed_certificate_byte_match=True,
  independent_canonical126_divergence_basis_and_actual350_by126_physical_Jacobi=True,
  independent_full25_GL_orbit_actions_and_weighted_scalar_relation=True,
  all_polynomial_Smith_left_matrix_right_identities_and_constant_unit_determinants_checked=True,
  nonunit_invariant_factors_by_weight={'6':'a-2/9','7':'a','8':'a'},
  all_a_affine_orbit_rank=17,all_a_projective_orbit_rank=16,
  exact_Jacobi_rank='73 at a=0; 74 at a=2/9; 75 otherwise',special_ranks=ranks,
  old_selected_group_minor_extra_factors_are_frame_artifacts=True,
  nonpositive7_group_frame_all_a=True,flatness_or_higher_ideal_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS physical full126 and all-a Jacobi/orbit polynomial Smith certificates')

if __name__=='__main__':main()
