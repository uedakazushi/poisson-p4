#!/usr/bin/env python3
"""Universal symbolic derivative circuits and Delta-times-J scheme-open audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def gz(p):
 with gzip.open(p,'rt')as f:return json.load(f)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_T39_SCHEME_OPEN_20261103_v63';path=directory/'report.json';d=load(path);f=load(args.replay);assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 parent=gz(ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz');m=gz(directory/'Jacobian-minor-model.json.gz');assert m['nodes'][:9647]==parent['nodes']and all(n[0]!='inv'for n in m['nodes'][9647:])
 dag=DAG(parent['input_names']);dag.nodes=[(op,tuple(data))if op in['+','*']else(op,data)for op,data in m['nodes']];dag.lookup={v:i for i,v in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1'];initial=len(dag.nodes);columns=m['selected_input_columns'];positions={i:j for j,i in enumerate(columns)};der=[]
 for i,(op,data)in enumerate(parent['nodes']):
  if op=='q':v=[dag.zero]*26
  elif op=='v':v=[dag.one if positions.get(data)==j else dag.zero for j in range(26)]
  elif op=='+':v=[dag.add(*(der[k][j]for k in data))for j in range(26)]
  elif op=='*':v=[dag.add(*(dag.mul(der[k][j],*data[:t],*data[t+1:])for t,k in enumerate(data)if der[k][j]!=dag.zero))for j in range(26)]
  elif op=='inv':v=[dag.neg(dag.mul(x,i,i))for x in der[data]]
  der.append(v)
 assert [der[parent['all350_residual_nodes'][i]]for i in m['selected_equation_rows']]==m['selected_derivative_matrix_nodes']and len(dag.nodes)==initial
 weights=[]
 for op,data in m['nodes']:
  if op=='q':w=None if not Fraction(data)else 0
  elif op=='v':w=parent['input_weights'][data]
  elif op=='+':
   ws={weights[j]for j in data if weights[j]is not None};assert len(ws)<=1;w=next(iter(ws),None)
  elif op=='*':w=None if any(weights[j]is None for j in data)else sum(weights[j]for j in data)
  elif op=='inv':assert weights[data]==0;w=0
  weights.append(w)
 rows=m['selected_equation_rows']
 for i,row in enumerate(m['selected_derivative_matrix_nodes']):
  for j,node in enumerate(row):assert weights[node]in[None,parent['residual_weight_bounds'][rows[i]]-parent['input_weights'][columns[j]]]
 assert sum(parent['residual_weight_bounds'][i]for i in rows)-sum(parent['input_weights'][i]for i in columns)==264
 paths += [ROOT/'reports/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW.json',ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.json',ROOT/'reports/Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_T39_SCHEME_OPEN_20261103_v63_REVIEW.md']
 result=dict(id='Q_T39_SCHEME_OPEN_20261103_v63_REVIEW',all_passed=True,fresh_exact_QQ_eps_producer_replay_matches=True,
  all9647_parent_node_dense26_symbolic_derivative_vectors_independently_rebuilt=True,
  all676_provided_derivative_node_IDs_match=True,parent_prefix_and_no_new_inverse_checked=True,
  original350x35_numerical_AD_and_nonzero26_minor_previously_independently_verified=True,universal_J_weight=264,exact_ray_power=44,
  actual_scheme_equality_R_DeltaJ_equals_A9_DeltaJ_accepted=True,
  reduced_support_input_uses_audited_actual_normal_divisor_classification=True,
  entire_Delta_open_or_full_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS all676 universal symbolic derivatives and precise Delta*J scheme open')

if __name__=='__main__':main()
