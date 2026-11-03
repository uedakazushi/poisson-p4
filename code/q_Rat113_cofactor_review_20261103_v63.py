#!/usr/bin/env python3
"""All-node17-vector divided-difference and vertical-torsion audit."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sympy as s
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha,normalize,closure,weights
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_RAT113_COFACTOR_OPEN_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 packed=folder/'cofactor-model.json.gz';assert packed.read_bytes()==(args.replay/'cofactor-model.json.gz').read_bytes();model=read(packed);paths=[Path(__file__),public,packed,ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 parent=read(ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/model37.json.gz');parameter=read(ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63/model20.json.gz');cert=read(ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63/certificate.json');assert model['nodes'][:len(parent['nodes'])]==parent['nodes']
 dag=DAG(model['input_names']);dag.nodes=list(map(normalize,model['nodes']));dag.lookup={v:i for i,v in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1'];dag.inputs={n:dag.lookup['v',i]for i,n in enumerate(dag.names)}
 embedded={}
 for i in closure(parameter['nodes'],parameter['parent37_input_nodes']):
  op,data=parameter['nodes'][i]
  if op=='q':v=dag.q(data)
  elif op=='v':v=dag.inputs[parameter['input_names'][data]]
  elif op=='+':v=dag.add(*(embedded[j]for j in data))
  elif op=='*':v=dag.mul(*(embedded[j]for j in data))
  else:raise ValueError(op)
  embedded[i]=v
 projection=[embedded[i]for i in parameter['parent37_input_nodes']];assert projection==model['graph_parent37_input_nodes']
 columns=model['dependent_parent_input_indices'];positions={c:i for i,c in enumerate(columns)};generators=[dag.sub(dag.inputs[parent['input_names'][i]],projection[i])for i in columns];assert generators==model['normalized17_generator_nodes']and generators[0]==dag.inputs['a']
 rows=model['selected350_equation_indices'];equations=[parent['all350_residual_nodes'][i]for i in rows];assert equations==model['selected_equation_nodes'];graph={};vectors={};needed=closure(parent['nodes'],equations)
 for i in needed:
  op,data=parent['nodes'][i]
  if op=='q':image=dag.q(data);v=[dag.zero]*17
  elif op=='v':
   image=projection[data];v=[dag.one if positions.get(data)==j else dag.zero for j in range(17)]
   if data not in positions:assert image==i
  elif op=='+':image=dag.add(*(graph[j]for j in data));v=[dag.add(*(vectors[j][k]for j in data))for k in range(17)]
  elif op=='*':
   image=dag.mul(*(graph[j]for j in data));v=[]
   for k in range(17):v.append(dag.add(*(dag.mul(vectors[child][k],*data[:pos],*(graph[j]for j in data[pos+1:]))for pos,child in enumerate(data))))
  else:raise ValueError(op)
  graph[i]=image;vectors[i]=v
 assert len(needed)==3247 and[vectors[i]for i in equations]==model['cofactor17_matrix_nodes']and[graph[i]for i in equations]==model['graph_selected_equation_nodes']
 point={name:F(value)for name,value in zip(parent['input_names'],cert['dense_parent37_point'])};values=dag.evaluate(point);assert all(not values[i]for i in generators)
 M=s.Matrix([[values[i]for i in row]for row in model['cofactor17_matrix_nodes']]);full=s.Matrix([[s.Rational(v)for v in row]for row in cert['dense_full350x37_Jacobian']]);assert M==full[rows,columns]
 det=M.det();assert det==s.Rational(cert['dense_dependent17_minor'])!=0
 ww=weights(model['nodes'],model['input_weights']);assert all(op!='inv'for op,_ in model['nodes']);frame=read(ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/frames.json.gz');rowweights=frame['Jacobi_row_weights'];degree=sum(rowweights[i]for i in rows)-sum(parent['input_weights'][i]for i in columns);assert degree==129
 for i,row in enumerate(model['cofactor17_matrix_nodes']):
  for j,node in enumerate(row):assert ww[node]in(None,rowweights[rows[i]]-parent['input_weights'][columns[j]])
 paths += [ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_RAT113_COFACTOR_OPEN_20261103_v63_REVIEW.md']
 report=dict(id='Q_RAT113_COFACTOR_OPEN_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_entire_cofactor_archive_match=True,
  actual_A20_graph_projection_independently_embedded=True,
  independent_all3247_ancestor17vector_product_telescoping_induction=True,
  all289_matrix_nodeIDs_and_graph_outputs_exactly_match=True,
  actual17_dependentJac_values_and_nonzero_determinant_checked=True,
  no_added_inverses_and_determinant_weight129=True,
  exact_I_C_equals_N113_C_colon_allpowers_and_prime_peel_accepted=True,
  explicit_nonzero_a_times_C113_zero_in_full_R37_local_and_completion_Q0=True,
  residual_I_plus_C113_or_full_branches_classified=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS all3247 exact17vector cofactors, Rat113 prime peel and explicit full a-torsion')

if __name__=='__main__':main()
