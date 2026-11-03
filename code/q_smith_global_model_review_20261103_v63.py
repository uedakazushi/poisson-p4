#!/usr/bin/env python3
"""Universal350, polynomial frames and exact finite37/35 DAG audit."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from ambient_CA_global_flag_review_20261103_v63 import IDS,action
from q_weighted_finite_model_review_20261103_v63 import derivative,add_product,clean
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1];a=s.Symbol('a');BASE=s.QQ.poly_ring(a);RA,aa=ring('a',s.QQ)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()))if p.suffix=='.gz'else json.loads(p.read_text())
def matrix(data):return s.Matrix([[s.sympify(v)for v in row]for row in data])
def dm(M):return DomainMatrix.from_Matrix(M).convert_to(BASE).to_dense()
def scalar(v):return RA.from_expr(s.sympify(v))
def plus(*values):
 out={}
 for value in values:
  for mon,c in value.items():out[mon]=out.get(mon,RA.zero)+c
 return{m:c for m,c in out.items()if c}
def times(left,right):
 out={}
 for m,c in left.items():
  for n,d in right.items():
   mon=tuple(sorted(m+n));assert len(mon)<=2;out[mon]=out.get(mon,RA.zero)+c*d
 return{m:c for m,c in out.items()if c}
def scaled(value,c):return{m:c*v for m,v in value.items()if c*v}
def normalize(node):return(node[0],tuple(node[1]))if node[0]in('+','*')else tuple(node)
def source_polynomials(nodes):
 values=[]
 for i,(op,data)in enumerate(nodes):
  if op=='q':v=scalar(data);value={():v}if v else{}
  elif op=='v':value={():aa}if data==0 else{(data-1,):RA.one}
  elif op=='+':assert all(j<i for j in data);value=plus(*(values[j]for j in data))
  elif op=='*':
   assert all(j<i for j in data);value={():RA.one}
   for j in data:value=times(value,values[j])
  else:raise ValueError(op)
  values.append(value)
 return values
def weights(nodes,inputs):
 out=[];deps=[]
 for i,(op,data)in enumerate(nodes):
  if op=='q':w=None if not Fraction(data)else 0;dep=0
  elif op=='v':w=inputs[data];dep=1<<data
  elif op in('+','*'):
   assert all(j<i for j in data);dep=0
   for j in data:dep|=deps[j]
   if op=='+':choices={out[j]for j in data if out[j]is not None};assert len(choices)<=1;w=next(iter(choices),None)
   else:w=None if any(out[j]is None for j in data)else sum(out[j]for j in data)
  elif op=='inv':assert data<i and out[data]==0 and deps[data]==1;w=0;dep=1
  else:raise ValueError(op)
  out.append(w);deps.append(dep)
 return out
def closure(nodes,outputs):
 used=set(outputs);pending=list(used)
 while pending:
  op,data=nodes[pending.pop()]
  for j in(data if op in('+','*')else[data]if op=='inv'else[]):
   if j not in used:used.add(j);pending.append(j)
 return sorted(used)
def dag_composition(source,model,pivots):
 dag=DAG(model['input_names']);dag.nodes=list(map(normalize,model['nodes']));dag.lookup={v:i for i,v in enumerate(dag.nodes)}
 dag.zero=dag.lookup[('q','0')];dag.one=dag.lookup[('q','1')]
 dag.inputs={name:dag.lookup[('v',i)]for i,name in enumerate(dag.names)}
 assigned=model['positive108_reconstruction_nodes']
 def clone(outputs,replacement):
  values={}
  for i in closure(source['nodes'],outputs):
   op,data=source['nodes'][i]
   if op=='q':v=dag.q(data)
   elif op=='v':v=dag.inputs['a']if data==0 else replacement[data-1]
   elif op=='+':v=dag.add(*(values[j]for j in data))
   elif op=='*':v=dag.mul(*(values[j]for j in data))
   else:raise ValueError(op)
   values[i]=v
  return values
 outputs=source['raw150_reconstruction_nodes']+source['all350_residual_nodes']+[p['source_node']for p in pivots]
 composed=clone(outputs,assigned)
 assert[composed[i]for i in outputs]==model['raw150_reconstruction_nodes']+model['all350_residual_nodes']+model['Smith_transformed75_residual_nodes']
 chosen=model['eliminated_Smith_pivots'];columns={p['column']for p in chosen}
 assert len(columns)==len(chosen)
 for i,c in enumerate(model['free_weighted_columns']):assert assigned[c]==dag.inputs[f'z{i}']
 for d in range(1,13):
  current=[p for p in chosen if p['weight']==d]
  if not current:continue
  replace=[dag.zero if p in columns and source['input_weights'][p+1]>=d else node for p,node in enumerate(assigned)]
  records=[next(q for q in pivots if q['column']==p['column'])for p in current]
  values=clone([p['source_node']for p in records],replace)
  for p in records:
   diagonal=s.sympify(p['diagonal'])
   inverse=dag.q(1/diagonal)if p['kind']=='unit'else dag.mul(dag.q(1/s.cancel(diagonal/a)),dag.inv(dag.inputs['a']))
   answer=dag.neg(dag.mul(inverse,values[p['source_node']]))
   assert answer==assigned[p['column']],('exact substitution node',p['column'])

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63';report=folder/'report.json';d=read(report);f=read(args.replay/'report.json')
 assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 paths=[Path(__file__),report,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',ROOT/'code/q_weighted_finite_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
 frames=read(folder/'frames.json.gz');source=read(folder/'source109.json.gz');C=matrix(frames['canonical_basis']);q=s.Matrix(list(map(s.sympify,frames['raw_Qa150'])))
 old=read(ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz');free=old['canonical_free_raw_indices'];J=matrix(old['full350_by126_Jacobi']);O=s.Matrix.hstack(*[action(i,j)*q for i,j in itertools.product(range(5),repeat=2)])
 assert C==matrix(old['canonical_basis'])and q==s.Matrix(list(map(s.sympify,old['raw_Qa150'])))
 B=matrix(frames['raw150_positive108_basis']);back=matrix(frames['positive108_from_raw150']);w=frames['positive108_weights'];target=frames['Jacobi_row_weights']
 assert B.shape==(150,108)and dm(back)*dm(B)==dm(s.eye(108))and dm(back)*dm(q)==dm(s.zeros(108,1))
 aindex=IDS[((1,4),(0,0))];assert q[aindex]==a and B[aindex,:]==s.zeros(1,108)
 degree=old['canonical_relative_weights'];rawdegree=[5-i-j+sum(mon)for((i,j),mon)in IDS]
 assert all(not B[i,j]or rawdegree[i]==w[j]for i in range(150)for j in range(108))
 ranks={'unit':0,'a':0,'u':0};offset=0
 for orbit,block in zip(frames['orbit_frames'],frames['Jacobi_complement_Smith_blocks']):
  cols=orbit['canonical_columns'];S=matrix(orbit['left']);complement=matrix(orbit['complement']);rank=orbit['removed_group_columns'];n=len(cols);assert block['weight']==orbit['weight']
  expected=s.zeros(n,n-rank)
  for i in range(n-rank):expected[rank+i,i]=1
  assert dm(S)*dm(complement)==dm(expected)
  U=matrix(block['left']);V=matrix(block['right']);saved=matrix(block['matrix']);rr=block['rows'];cc=block['columns']
  assert cc==list(range(offset,offset+len(cc)));offset+=len(cc)
  assert dm(C[:,cols])*dm(complement)*dm(V)==dm(B[:,cc])
  assert dm(J[rr,cols])*dm(complement)==dm(saved)
  diagonal=s.zeros(saved.rows,saved.cols)
  for i,v in enumerate(block['diagonal']):diagonal[i,i]=s.sympify(v)
  assert dm(U)*dm(saved)*dm(V)==dm(diagonal)
  for M in[S,U,V]:
   det=BASE.to_sympy(dm(M).det());assert det and not det.has(a)
  for v in block['diagonal']:
   v=s.sympify(v)
   if v:ranks['unit'if not v.has(a)else'a'if v.subs(a,0)==0 else'u']+=1
 assert offset==108 and ranks=={'unit':72,'a':2,'u':1}
 H=matrix(frames['group_cut_matrix']);G=matrix(frames['group_Lie_coefficient_matrix']);assert dm(H)*dm(B)==dm(s.zeros(17,108))
 positions=list(itertools.product(range(5),repeat=2))
 assert all(not G[k,j]or positions[k][0]-positions[k][1]==frames['group_Lie_weights'][j]for k in range(25)for j in range(17))
 assert dm(H)*dm(O)*dm(G)==dm(matrix(frames['group_response']))
 zero=s.Matrix.hstack(q.diff(a),B,O*G)[free,:];assert zero==matrix(frames['full126_action_Qa_matrix'])
 rowweights=frames['full126_action_row_weights'];colweights=frames['full126_action_column_weights']
 assert rowweights==degree and colweights==[0]+w+[-v for v in frames['group_Lie_weights']]and sum(rowweights)==sum(colweights)
 detprod=s.S.One
 for block in frames['full126_action_Qa_blocks']:
  det=BASE.to_sympy(dm(zero[block['rows'],block['columns']]).det());assert det==s.sympify(block['determinant'])and det and not det.has(a);detprod*=det
 assert s.cancel(detprod-s.sympify(frames['full126_action_reordered_determinant']))==0
 assert all(not zero[i,j]or rowweights[i]==colweights[j]for i in range(126)for j in range(126))
 # Independent universal physical coefficient expansion before substituting
 # any a-dependent frame, using five-variable polynomial differentiation.
 pairs=list(itertools.combinations(range(5),2));triples=list(itertools.combinations(range(5),3));mons3=list(itertools.combinations_with_replacement(range(5),3));P={ij:{}for ij in pairs}
 for(ij,mon),index in IDS.items():P[ij][tuple(mon.count(k)for k in range(5)),(index,)]=Fraction(1)
 def entry(i,j):return P[i,j]if i<j else{m:-v for m,v in P[j,i].items()}if i>j else{}
 deriv={(i,j,k):derivative(entry(i,j),k)for i in range(5)for j in range(5)for k in range(5)}
 physical=[]
 for i,j,k in triples:
  result=defaultdict(Fraction)
  for t in range(5):
   add_product(result,entry(i,t),deriv[j,k,t]);add_product(result,entry(j,t),deriv[k,i,t]);add_product(result,entry(k,t),deriv[i,j,t])
  for mon in mons3:
   exp=tuple(mon.count(v)for v in range(5));physical.append({p:c for(e,p),c in clean(result).items()if e==exp})
 assert physical==[{(i,j):Fraction(c)for i,j,c in row}for row in frames['physical_raw350_quadratics']]
 values=source_polynomials(source['nodes']);expectedraw=[]
 for i in range(150):
  p={():scalar(q[i])}if q[i]else{}
  p.update({(j,):scalar(B[i,j])for j in range(108)if B[i,j]});expectedraw.append(p)
 assert[values[i]for i in source['raw150_reconstruction_nodes']]==expectedraw
 equations=[]
 for row in physical:
  expected=plus(*(scaled(times(expectedraw[i],expectedraw[j]),scalar(c))for(i,j),c in row.items()));equations.append(expected)
 assert[values[i]for i in source['all350_residual_nodes']]==equations
 pivots=frames['source_Smith_pivots'];assert pivots==source['transformed_Smith_pivots']
 for p in pivots:
  polynomial=values[p['source_node']];col=p['column'];dweight=p['weight']
  assert polynomial.get((col,),RA.zero)==scalar(p['diagonal'])and ()not in polynomial
  assert all(sum(w[i]for i in mon)==dweight for mon in polynomial)
  assert all(all(w[i]<dweight for i in mon)for mon in polynomial if mon!=(col,))
 for n in[37,35]:
  model=read(folder/f'model{n}.json.gz');ww=weights(model['nodes'],model['input_weights'])
  assert len(model['input_names'])==n and len(model['free_weighted_columns'])==n-1
  assert sum(op=='inv'for op,_ in model['nodes'])==(n==35)
  assert all(ww[j]in(None,dw)for j,dw in zip(model['positive108_reconstruction_nodes'],w))
  assert all(ww[j]in(None,dw)for j,dw in zip(model['all350_residual_nodes'],target))
  assert all(ww[j]in(None,dw)for j,dw in zip(model['raw150_reconstruction_nodes'],rawdegree))
  dag_composition(source,model,pivots)
 paths += [ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.md']
 result=dict(id='Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW',all_passed=True,
  fresh_all_reports_and_four_compressed_archives_match=True,independent_universal_raw350_physical_quadratics=True,
  independent_full_source109_DAG_polynomial_semantics_in_QQ_a_w=True,
  all_polynomial_complement_and_new_Smith_unit_identities_checked=True,
  independent_all75_exact_diagonal_plus_strictly_lower_positive_weight_remainders=True,
  both_final_DAGs_all150_350_75_are_exact_source_substitutions=True,
  all72_or74_reconstruction_assignments_match_exact_triangular_source_recipe=True,
  model37_polynomial_all_a_and_model35_only_a_inverse=True,
  full126_action_Jacobian_including_qprime_Bprime_w_has_weight_zero_constant_unit=True,
  reordered_Qa_action_determinant=str(detprod),canonical_affine_formal_product_with17_at_all_Qa=True,
  global_action_etale_along_entire_identity_times_slice=True,
  global_group_cross_section_projective16_literal_at_arbitrary_point_or_full_branches_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS universal350/source109 polynomials, all-a37/a-unit35 exact DAGs and full action unit')

if __name__=='__main__':main()
