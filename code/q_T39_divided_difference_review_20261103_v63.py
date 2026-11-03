#!/usr/bin/env python3
"""Independent nodewise universal cofactor induction for the actual T39 open."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from q_local_q51_review_20261103_v63 import evaluate,firstjets

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def gz(p):
 with gzip.open(p,'rt')as f:return json.load(f)
def ancestry(nodes,outputs):
 used=set(outputs);todo=list(outputs)
 while todo:
  op,data=nodes[todo.pop()];children=data if op in['+','*']else[data]if op=='inv'else[]
  for i in children:
   if i not in used:used.add(i);todo.append(i)
 return used
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63';path=directory/'report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/q_local_q51_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 m=gz(directory/'cofactor-model.json.gz');parent=gz(ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz');param=gz(ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/T39-nine-parameter-model.json.gz')
 assert m['nodes'][:9647]==parent['nodes']and m['unchanged_parent_node_prefix']==9647
 assert all(n[0]!='inv'for n in m['nodes'][9647:])
 dag=DAG(m['input_names']);dag.nodes=[(n[0],tuple(n[1]))if n[0]in['+','*']else tuple(n)for n in m['nodes']];dag.lookup={v:i for i,v in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1'];dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
 columns=m['dependent_parent_input_indices'];free=param['parent_input_indices'];assert columns==[i for i in range(35)if i not in free]
 images=m['graph_parent35_input_nodes'];N=m['normalized_graph_generator_nodes']
 assert N==[dag.sub(dag.inputs[m['input_names'][i]],images[i])for i in columns]
 R,*xs=ring(','.join(param['input_names']),s.QQ);pv=evaluate(param['nodes'],xs,R,True);H=[pv[i]for i in param['parent35_input_nodes']]
 needed=ancestry(m['nodes'],images);values={};free_map={j:i for i,j in enumerate(free)}
 for i in sorted(needed):
  op,data=m['nodes'][i]
  if op=='q':v=R(s.Rational(data))
  elif op=='v':assert data in free_map;v=xs[free_map[data]]
  elif op=='+':v=sum((values[j]for j in data),R.zero)
  elif op=='*':
   v=R.one
   for j in data:v*=values[j]
  elif op=='inv':assert values[data].is_ground;v=R.one/values[data]
  values[i]=v
 assert [values[i]for i in images]==H
 # Dense26-vectors, reconstructed directly at EVERY parent ancestor, verify
 # the universal product telescope and inverse identity structurally.
 eq=m['selected_parent_equation_nodes'];rows=m['selected_parent350_equation_indices'];assert eq==[parent['all350_residual_nodes'][i]for i in rows]
 used=ancestry(parent['nodes'],eq);bar={};coef={};positions={i:j for j,i in enumerate(columns)};initial_nodes=len(dag.nodes)
 for i in sorted(used):
  op,data=parent['nodes'][i]
  if op=='q':b=dag.q(data);v=[dag.zero]*26
  elif op=='v':
   b=images[data];v=[dag.one if positions.get(data)==j else dag.zero for j in range(26)]
   if data not in positions:assert b==i
  elif op=='+':b=dag.add(*(bar[k]for k in data));v=[dag.add(*(coef[k][j]for k in data))for j in range(26)]
  elif op=='*':
   b=dag.mul(*(bar[k]for k in data));v=[]
   for j in range(26):v.append(dag.add(*(dag.mul(coef[k][j],*data[:t],*[bar[q]for q in data[t+1:]])for t,k in enumerate(data)if coef[k][j]!=dag.zero)))
  elif op=='inv':
   b=dag.inv(bar[data]);assert dag.nodes[b][0]=='q'and Fraction(dag.nodes[b][1]);v=[dag.neg(dag.mul(x,i,b))for x in coef[data]]
  bar[i]=b;coef[i]=v
 assert len(used)==2846 and [bar[i]for i in eq]==m['graph_selected_equation_nodes']
 assert [coef[i]for i in eq]==m['cofactor26_matrix_nodes']and len(dag.nodes)==initial_nodes
 # Independent full350x35 AD on the ORIGINAL parent at the actual ray point.
 actual=[Fraction(2,9)]+[Fraction(0)]*34
 for name,c in zip(['z11','z12','z23','z33'],[4,6,-4,1]):actual[parent['input_names'].index(name)]=Fraction(c)
 inputs=[(a,tuple(Fraction(int(i==j))for i in range(35)))for j,a in enumerate(actual)];jets=firstjets(parent['nodes'],inputs)
 jac=s.Matrix([jets[i][1]for i in parent['all350_residual_nodes']]);jc=load(ROOT/'reports/Q_T39_SCHEME_OPEN_20261103_v63/certificate.json')
 assert jac==s.Matrix(jc['full350x35_Jacobian_at_Rstar_eps1'])and all(jets[i][0]==0 for i in parent['all350_residual_nodes'])
 val=dag.evaluate(dict(zip(parent['input_names'],actual)));M=s.Matrix([[val[i]for i in row]for row in m['cofactor26_matrix_nodes']]);assert M==jac[rows,columns]
 det=str(DomainMatrix.from_Matrix(M).det());assert det==d['exact_cofactor_determinant_at_Rstar1']
 weights=[]
 for op,data in m['nodes']:
  if op=='q':w=None if not Fraction(data)else 0
  elif op=='v':w=parent['input_weights'][data]
  elif op=='+':
   ww={weights[j]for j in data if weights[j]is not None};assert len(ww)<=1;w=next(iter(ww),None)
  elif op=='*':w=None if any(weights[j]is None for j in data)else sum(weights[j]for j in data)
  elif op=='inv':assert weights[data]==0;w=0
  weights.append(w)
 for i,row in enumerate(m['cofactor26_matrix_nodes']):
  for j,node in enumerate(row):assert weights[node]in[None,parent['residual_weight_bounds'][rows[i]]-parent['input_weights'][columns[j]]]
 degree=sum(parent['residual_weight_bounds'][i]for i in rows)-sum(parent['input_weights'][i]for i in columns);assert degree==264
 paths += [ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.json',
  ROOT/'research/notes/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW.md']
 result=dict(id='Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW',all_passed=True,
  fresh_full_QQ9_untruncated_eps_and_dense_QQ_producer_replay_matches=True,
  graph35_polynomial_reprojection_independently_matches_verified_A9_model=True,
  every2846_parent_ancestor_universal_difference_identity_structurally_rebuilt=True,
  all676_provided_cofactor_nodes_match_independent_dense26_vector_induction=True,
  no_new_inverse_and_parent9647_prefix_verified=True,
  independent_original_full350x35_Jacobian_and_all676_Rstar_matrix_values_match=True,
  exact_cofactor_determinant_at_Rstar1=det,universal_weight=264,exact_curve_power=44,
  actual_localized_IDEAL_equality_I_C_equals_graph_N_C_accepted=True,
  normality_or_component_classification_needed=False,full35_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS universal exact F=M*N, actual T39 scheme open and eps^44')

if __name__=='__main__':main()
