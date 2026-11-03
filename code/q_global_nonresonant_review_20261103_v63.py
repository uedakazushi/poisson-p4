#!/usr/bin/env python3
"""Complete relation and exact DAG quotient audit of global34 Q model."""
import argparse
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha,normalize,closure,weights
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_SMITH_GLOBAL_NONRESONANT_MODEL_20261103_v63';p=folder/'report.json';d=read(p);f=read(args.replay/'report.json')
 assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths=[Path(__file__),p,ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file;paths.append(q)
 for name,h in d['evidence_sha256'].items():q=folder/name;assert sha(q)==h and q.read_bytes()==(args.replay/name).read_bytes();paths.append(q)
 model=read(folder/'model34.json.gz');elimination=read(folder/'elimination.json.gz');parent=read(ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/model35.json.gz')
 a=s.Symbol('a');K=s.QQ.frac_field(a);R,*z=ring(','.join(parent['input_names'][1:]),K);G=elimination['complete_relation_parent_node'];removed=elimination['parent_input_index'];values={}
 relevant=[G]+[parent['all350_residual_nodes'][r]for r,c in elimination['relation_parent350_linear_combination']]
 for i in closure(parent['nodes'],relevant):
  op,data=parent['nodes'][i]
  if op=='q':v=R(s.Rational(data))
  elif op=='v':v=R(K.gens[0])if data==0 else z[data-1]
  elif op=='+':v=sum((values[j]for j in data),R.zero)
  elif op=='*':
   v=R.one
   for j in data:v*=values[j]
  elif op=='inv':assert values[data].is_ground;v=R.one/values[data]
  else:raise ValueError(op)
  values[i]=v
 full=values[G];saved=R.from_dict({tuple(m):K.from_sympy(s.sympify(c))for m,c in elimination['complete_polynomial_terms']});assert full==saved
 combination=sum((R(K.from_sympy(s.sympify(c)))*values[parent['all350_residual_nodes'][r]]for r,c in elimination['relation_parent350_linear_combination']),R.zero)
 assert full==combination
 mon=tuple(int(j==removed-1)for j in range(34));diagonal=s.sympify(elimination['complete_diagonal']);assert full[mon]==K.from_sympy(diagonal)
 assert all(not m[removed-1]or m==mon for m in full)
 assert len(full)==d['complete_high_order_relation_terms']==259 and max(map(sum,full))==5
 for c in full.values():s.Poly(K.to_sympy(c),a,domain=s.QQ)
 rhs=full-R(K.from_sympy(diagonal))*z[removed-1]
 assert rhs==R.from_dict({tuple(m):K.from_sympy(s.sympify(c))for m,c in elimination['rhs_terms']})
 dag=DAG(model['input_names']);dag.nodes=list(map(normalize,model['nodes']));dag.lookup={v:i for i,v in enumerate(dag.nodes)};dag.zero=dag.lookup[('q','0')];dag.one=dag.lookup[('q','1')];dag.inputs={name:dag.lookup[('v',i)]for i,name in enumerate(dag.names)}
 mapping=model['parent_input_reconstruction_nodes'];composition=model['complete_parent_node_composition_map'];assert len(composition)==len(parent['nodes'])==8879
 def clone(required,replacement):
  values={}
  for i in closure(parent['nodes'],required):
   op,data=parent['nodes'][i]
   if op=='q':v=dag.q(data)
   elif op=='v':v=replacement[data]
   elif op=='+':v=dag.add(*(values[j]for j in data))
   elif op=='*':v=dag.mul(*(values[j]for j in data))
   elif op=='inv':v=dag.inv(values[data])
   values[i]=v
  return values
 actual=clone(list(range(len(parent['nodes']))),mapping);assert[actual[i]for i in range(len(parent['nodes']))]==composition
 for key in['raw150_reconstruction_nodes','all350_residual_nodes','positive108_reconstruction_nodes','Smith_transformed75_residual_nodes']:assert model[key]==[composition[i]for i in parent[key]]
 zero=list(mapping);zero[removed]=dag.zero;lower=clone([G],zero)[G]
 u=dag.add(dag.mul(dag.q(9),dag.inputs['a']),dag.q(-2));scale=s.sympify(elimination['diagonal_over_u'])
 solution=dag.neg(dag.mul(dag.q(1/scale),dag.inv(u),lower));assert solution==mapping[removed]==model['eliminated_input_expression_node']
 retained=model['parent_retained_input_indices'];assert retained==[i for i in range(35)if i!=removed]
 assert all(mapping[i]==dag.inputs[dag.names[j]]for j,i in enumerate(retained))
 ww=weights(model['nodes'],model['input_weights']);assert ww[solution]==6
 scalars=[];invargs=[]
 for op,data in model['nodes']:
  if op=='q':v=K.from_sympy(s.Rational(data))
  elif op=='v':v=K.gens[0]if data==0 else K.zero
  elif op=='+':v=sum((scalars[j]for j in data),K.zero)
  elif op=='*':
   v=K.one
   for j in data:v*=scalars[j]
  elif op=='inv':invargs.append(s.factor(K.to_sympy(scalars[data])));v=K.one/scalars[data]
  scalars.append(v)
 assert set(invargs)=={a,9*a-2}and len(invargs)==2
 paths += [ROOT/'reports/Q_SMITH_PIVOT_LIFT_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_SMITH_GLOBAL_NONRESONANT_MODEL_20261103_v63_REVIEW.md']
 report=dict(id='Q_SMITH_GLOBAL_NONRESONANT_MODEL_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_both_compressed_archives_match=True,
  independent_complete259_term_G_and_its_original350_polynomial_linear_combination=True,
  coefficients_in_QQ_a_and_only_unknown_term_a_minus_two_ninths_times_z30=True,
  all8879_parent_node_compositions_exactly_verified=True,
  eliminated_input_matches_exact_unit_linear_solution=True,
  all150_350_and108_outputs_retained=True,only_nonconstant_inverses=['a','9*a-2'],
  actual_entire_ideal_quotient_isomorphism_on_a_times_9a_minus2_open_accepted=True,
  higher_jet_truncation_flatness_or_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS complete G/original350 ideal identity and all8879 exact global34 composition nodes')

if __name__=='__main__':main()
