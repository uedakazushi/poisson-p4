#!/usr/bin/env python3
"""Independent T40 QQ9 kernel and universal divided-difference open audit."""
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
from q_T39_divided_difference_review_20261103_v63 import ancestry

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def gz(p):
 with gzip.open(p,'rt')as f:return json.load(f)
def weights(nodes,inputs):
 out=[]
 for op,data in nodes:
  if op=='q':w=None if not Fraction(data)else 0
  elif op=='v':w=inputs[data]
  elif op=='+':
   ws={out[i]for i in data if out[i]is not None};assert len(ws)<=1;w=next(iter(ws),None)
  elif op=='*':w=None if any(out[i]is None for i in data)else sum(out[i]for i in data)
  elif op=='inv':assert out[data]==0;w=0
  out.append(w)
 return out
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--finite-replay',type=Path,required=True)
 ap.add_argument('--cofactor-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();paths=[Path(__file__),ROOT/'code/q_local_q51_review_20261103_v63.py',ROOT/'code/q_T39_divided_difference_review_20261103_v63.py'];reports=[]
 for stem,replay in [('Q_T40_FINITE_QUOTIENT',args.finite_replay),('Q_T40_DIVIDED_DIFFERENCE_OPEN',args.cofactor_replay)]:
  directory=ROOT/f'reports/{stem}_20261103_v63';p=directory/'report.json';d=load(p);f=load(replay);assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
  paths.append(p);reports.append(d)
  for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file;paths.append(q)
  for file,h in d['evidence_sha256'].items():q=directory/file;assert sha(q)==h,file;paths.append(q)
 finite,cofactor=reports;parent=gz(ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz');source=gz(ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/graph-model.json.gz');par=gz(ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/T40-nine-parameter-model.json.gz');cert=load(ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/certificate.json');m=gz(ROOT/'reports/Q_T40_DIVIDED_DIFFERENCE_OPEN_20261103_v63/cofactor-model.json.gz')
 center=Fraction(-8,9);inputjets=[(center if j==0 else Fraction(0),tuple(Fraction(int(i==j))for i in range(35)))for j in range(35)];jets=firstjets(source['nodes'],inputjets)
 jac=s.Matrix([jets[i][1]for i in source['Q52_graph_relation_nodes']]);assert jac==s.Matrix(cert['Q52_graph_Jacobian45x35'])
 assert jac[cert['Q52_graph_selected_rows'],cert['Q52_graph_selected_columns']].det()==s.Rational(cert['Q52_graph_minor'])==s.Rational(83385908498332845,3378528649216)
 assert s.Matrix([jets[i][1]for i in parent['all350_residual_nodes']]).rank()==1
 R,*xs=ring(','.join(par['input_names']),s.QQ);pv=evaluate(par['nodes'],xs,R,True);H=[pv[i]for i in par['parent35_input_nodes']]
 assert H[0]==R(s.Rational(-8,9))and par['input_names']==['z22','z26','z27','z28','z29','z30','z31','z32','z33']
 ov=evaluate(source['nodes'],H,R,True);assert all(ov[i]==R.zero for i in[*source['parent_all350_residual_nodes'],*source['Q52_graph_relation_nodes']])
 assert [ov[i]for i in source['raw150_parent_reconstruction_nodes']]==[pv[i]for i in par['raw150_nodes']]
 for j,index in enumerate(par['parent_input_indices']):
  assert H[index]==xs[j];raw=parent['positive_raw_indices'][parent['free_positive_columns'][index-1]];assert ov[source['raw150_parent_reconstruction_nodes'][raw]]==xs[j]
 sw=weights(source['nodes'],parent['input_weights'])
 for block in par['blocks']:
  M=s.Matrix(block['matrix']);inv=s.Matrix(block['inverse']);assert M*inv==s.eye(M.rows)and M.det()==s.Rational(block['determinant'])!=0
  assert M==jac[block['rows'],block['columns']]
  assert all(sw[source['Q52_graph_relation_nodes'][i]]==block['weight']for i in block['rows'])
  assert all(parent['input_weights'][i]==block['weight']for i in block['columns'])
 assert [len(b['rows'])for b in par['blocks']]==[2,4,5,6,5,3]
 dag=DAG(parent['input_names']);dag.nodes=[(op,tuple(data))if op in['+','*']else(op,data)for op,data in m['nodes']];dag.lookup={v:i for i,v in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1'];dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
 assert m['nodes'][:9647]==parent['nodes']and all(n[0]!='inv'for n in m['nodes'][9647:]);columns=m['dependent_parent_input_indices'];free=par['parent_input_indices'];assert columns==[i for i in range(35)if i not in free]
 images=m['graph_parent35_input_nodes'];assert m['normalized_graph_generator_nodes']==[dag.sub(dag.inputs[parent['input_names'][i]],images[i])for i in columns]
 graph={};fm={j:i for i,j in enumerate(free)}
 for i in sorted(ancestry(m['nodes'],images)):
  op,data=m['nodes'][i]
  if op=='q':v=R(s.Rational(data))
  elif op=='v':assert data in fm;v=xs[fm[data]]
  elif op=='+':v=sum((graph[j]for j in data),R.zero)
  elif op=='*':
   v=R.one
   for j in data:v*=graph[j]
  elif op=='inv':assert graph[data].is_ground;v=R.one/graph[data]
  graph[i]=v
 assert [graph[i]for i in images]==H
 equations=m['selected_parent_equation_nodes'];rows=m['selected_parent350_equation_indices'];assert equations==[parent['all350_residual_nodes'][i]for i in rows]
 used=ancestry(parent['nodes'],equations);bar={};vectors={};positions={i:j for j,i in enumerate(columns)};original_count=len(dag.nodes)
 for i in sorted(used):
  op,data=parent['nodes'][i]
  if op=='q':b=dag.q(data);v=[dag.zero]*26
  elif op=='v':b=images[data];v=[dag.one if positions.get(data)==j else dag.zero for j in range(26)]
  elif op=='+':b=dag.add(*(bar[k]for k in data));v=[dag.add(*(vectors[k][j]for k in data))for j in range(26)]
  elif op=='*':
   b=dag.mul(*(bar[k]for k in data));v=[dag.add(*(dag.mul(vectors[k][j],*data[:t],*[bar[q]for q in data[t+1:]])for t,k in enumerate(data)if vectors[k][j]!=dag.zero))for j in range(26)]
  elif op=='inv':b=dag.inv(bar[data]);assert dag.nodes[b][0]=='q';v=[dag.neg(dag.mul(x,i,b))for x in vectors[data]]
  bar[i]=b;vectors[i]=v
 assert len(used)==2801 and [vectors[i]for i in equations]==m['cofactor26_matrix_nodes']and[bar[i]for i in equations]==m['graph_selected_equation_nodes']and len(dag.nodes)==original_count
 actual=[Fraction(str(sum(h.values(),s.QQ.zero)))for h in H];assert actual==list(map(Fraction,cert['dense_parent35_point']))
 densejets=firstjets(parent['nodes'],[(a,tuple(Fraction(int(i==j))for i in range(35)))for j,a in enumerate(actual)]);fulljac=s.Matrix([densejets[i][1]for i in parent['all350_residual_nodes']]);assert fulljac[:,columns]==s.Matrix(cert['dense_full350x26_dependent_Jacobian'])
 values=dag.evaluate(dict(zip(parent['input_names'],actual)));M=s.Matrix([[values[i]for i in row]for row in m['cofactor26_matrix_nodes']]);assert M==fulljac[rows,columns]
 det=str(DomainMatrix.from_Matrix(M).det());assert det==cofactor['exact_cofactor_determinant_at_dense_T40_point']
 ww=weights(m['nodes'],parent['input_weights'])
 for r,row in enumerate(m['cofactor26_matrix_nodes']):
  for c,node in enumerate(row):assert ww[node]in[None,parent['residual_weight_bounds'][rows[r]]-parent['input_weights'][columns[c]]]
 degree=sum(parent['residual_weight_bounds'][i]for i in rows)-sum(parent['input_weights'][i]for i in columns);assert degree==275
 paths += [ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.json',ROOT/'reports/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW.json',ROOT/'reports/Q_GROUP_SLICE_GLOBAL_UNIT_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_T40_GRAPH_OPEN_20261103_v63_REVIEW.md']
 result=dict(id='Q_T40_GRAPH_OPEN_20261103_v63_REVIEW',all_passed=True,fresh_finite_and_cofactor_producers_match=True,
  independent_QQ9_substitution_into_original_full350_plus45_graph_equations_zero=True,
  all150_raw_and_physical9_left_inverse_and6_exact_constant_inverse_blocks_checked=True,
  actual_R_by_Q52_quotient_equals_QQ9=True,center='-8/9',center_full_Jacobian_rank=1,center_regular_claimed=False,
  all2801_parent_ancestor_dense26_cofactor_vectors_structurally_rebuilt=True,all676_actual_cofactor_node_IDs_match=True,
  independent_dense_parent_fullJacobian_dependent26_columns_and_exact_minor_checked=True,
  no_new_inverse_and_weight275_checked=True,actual_scheme_open_I_C40_equals_N40_C40_accepted=True,
  nonzero_weighted_lambda_ray_all_regular9_canonical26_projective25=True,
  new_branch_through_Qstar_or_full_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS T40 exactA9 quotient and universal actual scheme open, separate center')

if __name__=='__main__':main()
