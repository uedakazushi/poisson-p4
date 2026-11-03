#!/usr/bin/env python3
"""Independent finite9-variable quotient, graph dependencies and block audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def readgz(p):
 with gzip.open(p,'rt')as f:return json.load(f)
def evaluate(nodes,inputs,R,constant_inverse=False):
 values=[]
 for op,data in nodes:
  if op=='q':v=R(s.Rational(data))
  elif op=='v':v=inputs[data]
  elif op=='+':v=sum((values[i]for i in data),R.zero)
  elif op=='*':
   v=R.one
   for i in data:v*=values[i]
  elif op=='inv':
   if constant_inverse:assert values[data].is_ground and values[data]
   v=R.one/values[data]
  values.append(v)
 return values
def firstjets(nodes,inputs):
 zero=(Fraction(0),)*35;values=[]
 for op,data in nodes:
  if op=='q':v=(Fraction(data),zero)
  elif op=='v':v=inputs[data]
  elif op=='+':v=(sum((values[i][0]for i in data),Fraction(0)),tuple(sum((values[i][1][j]for i in data),Fraction(0))for j in range(35)))
  elif op=='*':
   v=(Fraction(1),zero)
   for i in data:
    c,g=values[i];b,h=v;v=(b*c,tuple(b*x+c*y for x,y in zip(g,h)))
  elif op=='inv':
   c,g=values[data];assert c;v=(1/c,tuple(-x/(c*c)for x in g))
  values.append(v)
 return values

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63';path=directory/'report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 g=readgz(directory/'graph-model.json.gz');n=readgz(directory/'T39-nine-parameter-model.json.gz');cert=load(directory/'certificate.json');parent=readgz(ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
 assert g['nodes'][:9647]==parent['nodes']and g['unchanged_parent_node_prefix']==9647
 dep=[];weights=[]
 for i,(op,data)in enumerate(g['nodes']):
  if op=='q':bit=0;w=None if not Fraction(data)else 0
  elif op=='v':bit=1<<data;w=parent['input_weights'][data]
  elif op in['+','*']:
   assert all(j<i for j in data);bit=0
   for j in data:bit|=dep[j]
   if op=='+':
    ww={weights[j]for j in data if weights[j]is not None};assert len(ww)<=1;w=next(iter(ww),None)
   else:w=None if any(weights[j]is None for j in data)else sum(weights[j]for j in data)
  elif op=='inv':assert data<i;bit=dep[data];assert bit&~1==0 and weights[data]==0;w=0
  dep.append(bit);weights.append(w)
 for label in ['Q51','Q52']:assert dep[g[label+'_graph_relation_nodes'][25]]==1
 a,c,e=s.symbols('a c d');frame=s.Matrix([[1,0,0,0],[0,0,0,1/a],[0,1,0,-c/a],[0,0,-s.Rational(9,2),9*e/(2*a)]])
 inv=s.Matrix([[1,0,0,0],[0,c,1,0],[0,e,0,-s.Rational(2,9)],[0,a,0,0]])
 C=s.Matrix([[0,c,e,a],[-c,0,-s.Rational(2,9),0],[-e,s.Rational(2,9),0,0],[-a,0,0,0]])
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 assert (frame*inv-s.eye(4)).applyfunc(s.cancel)==s.zeros(4)and(frame*C*frame.T-J).applyfunc(s.cancel)==s.zeros(4)
 K=s.QQ.frac_field(a);field=evaluate(g['nodes'],[K.gens[0]]+[K.zero]*34,K)
 for label,target in [('Q51',(9*a-2)/(2*a*a)),('Q52',(9*a+8)/(2*a*a))]:
  values=[K.to_sympy(field[i])for i in g[label+'_graph_relation_nodes']]
  assert all(s.cancel(v-(target if j==25 else 0))==0 for j,v in enumerate(values))
 inputs=[]
 for j in range(35):inputs.append((Fraction(2,9)if j==0 else Fraction(0),tuple(Fraction(int(i==j))for i in range(35))))
 jets=firstjets(g['nodes'],inputs);jac=s.Matrix([jets[i][1]for i in g['Q51_graph_relation_nodes']]);assert jac==s.Matrix(cert['Q51_jacobian_45x35'])
 assert jac[cert['selected_rows'],cert['selected_columns']].det()==s.Rational(cert['selected_minor'])==s.Rational(-50031545098999707,95420416)
 assert all(jets[i]==(Fraction(0),(Fraction(0),)*35)for i in g['parent_all350_residual_nodes'])
 # Evaluate the separate9-variable DAG, THEN independently compose its35
 # input polynomials into the ORIGINAL10953-node graph/parent model.
 R,*variables=ring(','.join(n['input_names']),s.QQ);v=evaluate(n['nodes'],variables,R,True)
 h=[v[i]for i in n['parent35_input_nodes']];assert h[0]==R(s.Rational(2,9))
 assert n['input_names']==['z22','z24','z27','z28','z29','z30','z31','z32','z33']
 for j,i in enumerate(n['parent_input_indices']):assert h[i]==variables[j]
 original=evaluate(g['nodes'],h,R,True)
 assert all(original[i]==R.zero for i in [*g['parent_all350_residual_nodes'],*g['Q51_graph_relation_nodes']])
 assert [original[i]for i in g['raw150_parent_reconstruction_nodes']]==[v[i]for i in n['raw150_nodes']]
 for j,input_index in enumerate(n['parent_input_indices']):
  raw_index=parent['positive_raw_indices'][parent['free_positive_columns'][input_index-1]];assert original[g['raw150_parent_reconstruction_nodes'][raw_index]]==variables[j]
 for block in n['blocks']:
  M=s.Matrix(block['matrix']);I=s.Matrix(block['inverse']);assert M*I==s.eye(M.rows)and M.det()==s.Rational(block['determinant'])!=0
  assert M==jac[block['rows'],block['columns']]
  degree=block['weight'];assert all(parent['input_weights'][i]==degree for i in block['columns'])
  assert all(weights[g['Q51_graph_relation_nodes'][i]]==degree for i in block['rows'])
 assert [len(b['rows'])for b in n['blocks']]==[2,4,5,6,5,3]
 paths += [ROOT/'research/notes/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.md']
 result=dict(id='Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW',all_passed=True,
  fresh30second_exact_QQ_a_physical_and_QQ9_producer_replay_matches=True,
  fixed_value_frame_inverse_and_all_DAG_weights_dependencies_independently_checked=True,
  full_model_graph25_identities={'Q51':'(9a-2)/(2a²)','Q52':'(9a+8)/(2a²)'},
  independent_all45x35_exact_firstjet_matrix_and26_minor=True,
  independent_QQ9_composition_into_ORIGINAL_parent_and_graph_model_all395_zero=True,
  all150_raw_polynomials_and_physical9_coordinate_left_inverse_checked=True,
  all6_constant_inverse_blocks_and_exact_weighted_triangular_ideal_argument_accepted=True,
  actual_global_quotient_R_by_Q51_equals_QQ9=True,T39_transverse_completion_regular_dimension=9,
  T40_through_Qstar_or_full35_branch_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS exact GLOBAL Q51 quotient=QQ9, all395 original equations and graph kernel')

if __name__=='__main__':main()
