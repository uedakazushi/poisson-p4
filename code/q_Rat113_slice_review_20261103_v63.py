#!/usr/bin/env python3
"""Independent rank3 Schur scheme, exact A20 and full-parent component audit."""
import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from ambient_CA_global_flag_review_20261103_v63 import IDS
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha,normalize,weights
from q_local_q51_review_20261103_v63 import evaluate
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

def jets(nodes,point):
 values=[]
 for op,data in nodes:
  if op=='q':value={-1:F(data)}if F(data)else{}
  elif op=='v':value={data:F(1)};value.update({-1:F(point[data])}if point[data]else{})
  elif op=='+':
   value={}
   for j in data:
    for k,c in values[j].items():value[k]=value.get(k,F(0))+c
   value={k:c for k,c in value.items()if c}
  elif op=='*':
   value={-1:F(1)}
   for j in data:
    right=values[j];a=value.get(-1,F(0));b=right.get(-1,F(0));new={-1:a*b}if a*b else{}
    for k in(set(value)|set(right))-{ -1 }:
     c=a*right.get(k,F(0))+b*value.get(k,F(0))
     if c:new[k]=c
    value=new
  else:raise ValueError(op)
  values.append(value)
 return values
def rational_matrix(data):return s.Matrix([[s.Rational(v)for v in row]for row in data])

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths=[Path(__file__),public,ROOT/'code/q_local_q51_review_20261103_v63.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
 parent=read(ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/model37.json.gz');source=read(folder/'source-Schur.json.gz');model=read(folder/'model20.json.gz');fixed=read(folder/'fixed-Casimir-linear.json.gz');cert=read(folder/'certificate.json')
 assert source['nodes'][:len(parent['nodes'])]==parent['nodes']
 dag=DAG(parent['input_names']);dag.nodes=list(map(normalize,source['nodes']));dag.lookup={v:i for i,v in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1'];dag.inputs={n:dag.lookup['v',i]for i,n in enumerate(dag.names)}
 raw=parent['raw150_reconstruction_nodes'];mons=list(itertools.combinations_with_replacement(range(5),2));M=[]
 for i in range(5):
  for mon in mons:M.append([dag.zero if i==j else raw[IDS[((i,j),mon)]]if i<j else dag.neg(raw[IDS[((j,i),mon)]])for j in range(5)])
 assert M==source['contraction75x5_nodes'];rr=source['pivot_rows'];cc=source['pivot_columns'];assert rr==[45,30,50]and cc==[2,3,4]
 A=[[M[i][j]for j in cc]for i in rr];c=A[2][0];assert A==[[dag.q(F(2,9)),dag.zero,dag.zero],[dag.zero,dag.q(F(-2,9)),dag.zero],[c,dag.zero,dag.one]]
 E=[[M[i][j]for j in range(2)]for i in rr];K=[[dag.one,dag.zero],[dag.zero,dag.one],[dag.neg(dag.mul(dag.q(F(9,2)),E[0][j]))for j in range(2)],[dag.mul(dag.q(F(9,2)),E[1][j])for j in range(2)],[dag.sub(dag.mul(dag.q(F(9,2)),c,E[0][j]),E[2][j])for j in range(2)]]
 assert K==source['kernel5x2_nodes']
 roots=[dag.add(*(dag.mul(row[k],K[k][j])for k in range(5)))for row in M for j in range(2)];assert roots==source['Schur150_nodes']and roots[121]==dag.neg(dag.inputs['a'])
 nodeweights=weights(source['nodes'],source['input_weights']);assert all(nodeweights[i]in(None,v)for i,v in zip(roots,source['Schur150_weights']))
 cc0=s.Symbol('c');pivot=s.Matrix([[s.Rational(2,9),0,0],[0,-s.Rational(2,9),0],[cc0,0,1]]);assert pivot.det()==-s.Rational(4,81)
 zerojets=jets(source['nodes'],[F(0)]*37);Jac=s.Matrix([[s.Rational(zerojets[i].get(j,F(0)))for j in range(37)]for i in roots]);assert Jac==rational_matrix(cert['Schur150x37_Jacobian_at_Q0'])
 selected=cert['Schur_selected_rows'];dependent=cert['dependent_parent_input_columns'];assert 121 in selected and 0 in dependent and Jac[selected,dependent].det()==s.Rational(cert['Schur17_minor'])!=0
 assert DomainMatrix.from_Matrix(Jac).rank()==17
 assert [len(b['rows'])for b in model['triangular_blocks']]==[1,2,3,4,4,2]
 for block in model['triangular_blocks']:
  B=Jac[block['rows'],block['columns']];inverse=rational_matrix(block['inverse']);assert B==rational_matrix(block['matrix'])and B*inverse==s.eye(B.rows)
  assert all(parent['input_weights'][j]==block['weight']for j in block['columns'])and all(source['Schur150_weights'][i]==block['weight']for i in block['rows'])
 R,*v=ring(','.join(model['input_names']),s.QQ);values=evaluate(model['nodes'],v,R);H=[values[i]for i in model['parent37_input_nodes']];assert H[0]==0
 assert H==[R.from_expr(s.sympify(p))for p in cert['parent37_polynomials']]
 for j,i in enumerate(model['parent37_input_indices']):assert H[i]==v[j]
 pvalues=evaluate(parent['nodes'],H,R);actualraw=[pvalues[i]for i in parent['raw150_reconstruction_nodes']]
 assert actualraw==[values[i]for i in model['raw150_reconstruction_nodes']]and actualraw==[R.from_expr(s.sympify(p))for p in cert['raw150_polynomials']]
 assert all(not pvalues[i]for i in parent['all350_residual_nodes'])
 MM=[]
 for i in range(5):
  for mon in mons:MM.append([R.zero if i==j else actualraw[IDS[((i,j),mon)]]if i<j else-actualraw[IDS[((j,i),mon)]]for j in range(5)])
 E=[[MM[i][j]for j in range(2)]for i in rr];c=MM[50][2];KK=[[R.one,R.zero],[R.zero,R.one],[-s.Rational(9,2)*v for v in E[0]],[s.Rational(9,2)*v for v in E[1]],[s.Rational(9,2)*c*E[0][j]-E[2][j]for j in range(2)]]
 assert all(sum((row[k]*KK[k][j]for k in range(5)),R.zero)==0 for row in MM for j in range(2))
 assert [p for row in KK for p in row]==[values[i]for i in model['kernel5x2_nodes']]
 assert all(not values[i]for i in model['all150_Schur_nodes']+model['all350_residual_nodes'])and sum(map(len,actualraw))==69
 # Fixed-W exact relative Poincare splitting over arbitrary QQ algebras.
 cubicmons=[m for m in itertools.combinations_with_replacement(range(5),3)if any(i>=2 for i in m)];T=s.zeros(150,31);L=s.zeros(31,150);lookup={m:i for i,m in enumerate(cubicmons)}
 for pair,active,sign in[((2,3),4,1),((2,4),3,-1),((3,4),2,1)]:
  for mon in mons:
   i=IDS[pair,mon];m=tuple(sorted(mon+(active,)));T[i,lookup[m]]=sign*m.count(active);L[lookup[m],i]+=s.Rational(sign,1+sum(k>=2 for k in mon))
 assert T==rational_matrix(fixed['H31_to_raw150'])and L==rational_matrix(fixed['raw150_to_H31_relative_integration'])and L*T==s.eye(31)
 div=s.zeros(25,150)
 for((i,j),mon),column in IDS.items():
  for k,target,sign in[(i,j,1),(j,i,-1)]:
   if k in mon:remaining=list(mon);remaining.remove(k);div[5*target+remaining[0],column]+=sign*mon.count(k)
 active=fixed['active45_raw_indices'];assert div*T==s.zeros(25,31)and div[:,active]==rational_matrix(fixed['divergence25_by45'])and div[:,active].rank()==14
 frames=read(ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/frames.json.gz')
 for polynomial in frames['physical_raw350_quadratics']:
  terms={}
  for i,j,c in polynomial:
   for k in range(31):
    if T[i,k]:
     for l in range(31):
      if T[j,l]:mon=tuple(sorted([k,l]));terms[mon]=terms.get(mon,0)+c*T[i,k]*T[j,l]
  assert all(v==0 for v in terms.values())
 point=[sum((F(str(c))for c in p.values()),F(0))for p in H]
 dense=jets(parent['nodes'],point);full=s.Matrix([[s.Rational(dense[i].get(j,F(0)))for j in range(37)]for i in parent['all350_residual_nodes']]);assert full==rational_matrix(cert['dense_full350x37_Jacobian'])and DomainMatrix.from_Matrix(full).rank()==17
 assert full[cert['dense_selected_equation_rows'],dependent].det()==s.Rational(cert['dense_dependent17_minor'])!=0
 paths += [ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_RAT113_GLOBAL_SLICE_20261103_v63_REVIEW.md']
 report=dict(id='Q_RAT113_GLOBAL_SLICE_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_all_evidence_match=True,independent_symbolic_contraction_K_and150_Schur_DAG_identities=True,
  universal_constant_rank3_minor='-4/81',literal_minus_a_Schur_equation=True,
  independent150_by37_centerJac_and_all6_constant_block_inverses=True,
  exact_QQ20_actual_parent350_and150_Schur_and_free_coordinate_leftinverse=True,
  fixed_Casimir_T31_L_inverse_divergence14_and_all350_physical_quadratic_identities=True,
  scheme_recognition_over_arbitrary_QQ_algebras_including_nilpotents_accepted=True,
  independent_dense_actual350_by37_Jac_rank17_and_minor=True,
  A20_graph_is_actual_minimal_component_and_smooth_Rat113_branch=True,
  full_R37_local_and_completed_a_nonflat_at_Qzero_by_vertical_minimal_Ass_accepted=True,
  full_Qzero_regular_or_all_branches_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS actual scheme Rat113=A20,relative31 splitting and full-parent vertical minimal component')

if __name__=='__main__':main()
