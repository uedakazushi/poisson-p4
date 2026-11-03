#!/usr/bin/env python3
"""Independent Levi-Civita derivative circuits, integer minor and geometric scope."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import math
import sympy as s
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def node_tuple(n):return(n[0],tuple(n[1]))if n[0]in['+','*']else tuple(n)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63';path=directory/'report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 with gzip.open(directory/'matrix-model.json.gz','rt')as f:m=json.load(f)
 parentpath=ROOT/m['parent_model']
 with gzip.open(parentpath,'rt')as f:parent=json.load(f)
 assert m['nodes'][:len(parent['nodes'])]==parent['nodes']and m['unchanged_parent_node_prefix']==len(parent['nodes'])==9647
 assert all(n[0]!='inv'for n in m['nodes'][9647:])
 pairs=list(itertools.combinations(range(5),2));mon2=list(itertools.combinations_with_replacement(range(5),2));ids={(ij,mon):15*i+j for i,ij in enumerate(pairs)for j,mon in enumerate(mon2)}
 mon4=list(map(tuple,m['quartic_monomial_exponents']));form={(h,a):defaultdict(Fraction)for h in range(5)for a in mon4}
 # F=(1/8)sum epsilon_ijklm*x_i*P_jk*P_lm. Differentiate
 # every ambient coordinate BEFORE restriction to x3=x4=0.
 for order in itertools.permutations(range(5)):
  i,j,k,l,t=order;sign=(-1)**sum(order[a]>order[b]for a in range(5)for b in range(a+1,5));sign*=1 if j<k else-1;sign*=1 if l<t else-1
  for a in mon2:
   for b in mon2:
    powers=[i,*a,*b]
    for derivative in set(powers):
     remaining=powers.copy();remaining.remove(derivative)
     if max(remaining)>2:continue
     exp=tuple(remaining.count(q)for q in range(3));raw=tuple(sorted([ids[(tuple(sorted([j,k])),a)],ids[(tuple(sorted([l,t])),b)]]))
     form[derivative,exp][raw]+=Fraction(sign*powers.count(derivative),8)
 form={key:{mon:int(c)for mon,c in row.items()if c}for key,row in form.items()}
 saved=[{tuple(mon):c for mon,c in row}for row in m['raw150_quadratic_formulas']]
 assert [form[h,a]for h in range(5)for a in mon4]==saved
 dag=DAG(m['input_names']);dag.nodes=list(map(node_tuple,m['nodes']));dag.lookup={n:i for i,n in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
 raw=parent['raw150_reconstruction_nodes'];start=len(dag.nodes)
 for h in range(5):
  for j,a in enumerate(mon4):
   node=dag.add(*(dag.mul(dag.q(c),raw[i],raw[k])for(i,k),c in sorted(form[h,a].items())))
   assert node==m['partial_coefficient_nodes'][h][j]
   assert dag.mul(dag.q(162),node)==m['scaled_partial_coefficient_nodes'][h][j]
 assert len(dag.nodes)==start
 values=[];weights=[];inputs=[Fraction(2,9)]+[Fraction(0)]*34
 for name,c in zip(['z11','z12','z23','z33'],[4,6,-4,1]):inputs[m['input_names'].index(name)]=Fraction(c)
 for op,data in dag.nodes:
  if op=='q':v=Fraction(data);w=None if not v else 0
  elif op=='v':v=inputs[data];w=m['input_weights'][data]
  elif op=='+':
   v=sum((values[i]for i in data),Fraction(0));ww={weights[i]for i in data if weights[i]is not None};assert len(ww)<=1;w=next(iter(ww),None)
  elif op=='*':v=math.prod(values[i]for i in data);w=None if any(weights[i]is None for i in data)else sum(weights[i]for i in data)
  elif op=='inv':v=1/values[data];assert weights[data]==0;w=0
  values.append(v);weights.append(w)
 assert all(values[i]==0 for i in parent['all350_residual_nodes'])
 mon6=[tuple(m.count(i)for i in range(3))for m in itertools.combinations_with_replacement(range(3),6)]
 mon10=[tuple(m.count(i)for i in range(3))for m in itertools.combinations_with_replacement(range(3),10)];positions={a:i for i,a in enumerate(mon4)}
 matrix=[];rowweights=[]
 for h in range(5):
  for shift in mon6:
   row=[];rowweights.append(h-shift[1]-2*shift[2])
   for mon in mon10:
    diff=tuple(x-y for x,y in zip(mon,shift));node=m['scaled_partial_coefficient_nodes'][h][positions[diff]]if min(diff)>=0 else dag.zero
    row.append(node)
   matrix.append(row)
 assert matrix==m['matrix_nodes']and [list(a)for a in mon10]==m['matrix_columns']
 for h in range(5):
  for j,a in enumerate(mon4):assert weights[m['scaled_partial_coefficient_nodes'][h][j]]in[None,h+a[1]+2*a[2]]
 selected=m['selected_minor_rows'];M=s.Matrix([[s.Rational(values[node])for node in matrix[i]]for i in selected]);assert all(v.q==1 for v in M)
 det=int(DomainMatrix.from_Matrix(M).convert_to(s.ZZ).det());assert str(det)==d['exact_integer_minor_at_Rstar_eps_one']and det%32003==9347
 weight=sum(rowweights[i]for i in selected)+sum(a[1]+2*a[2]for a in mon10);assert weight==354 and weight//6==59
 paths += [ROOT/'research/notes/Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63_REVIEW.md']
 result=dict(id='Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63_REVIEW',all_passed=True,
  fresh_untruncated_QQ_eps_producer_replay_matches=True,independent_all75_full_Levi_Civita_ambient_derivative_formulas_match=True,
  all75_actual_DAG_coefficient_and_scaled_output_identities_structurally_checked=True,
  unchanged_parent9647_nodes_and_no_new_inverse_verified=True,all140x66_matrix_nodes_independently_rebuilt=True,
  selected_integer66_minor_independently_computed=str(det),minor_residue_mod32003=9347,universal_weight=354,
  actual_resonance_axis_identity='Delta=integer_minor*eps^59',normality_on_minor_open_geometrically_proved=True,
  reduced_support_T39_containment_uses_audited_normal_ACTUAL_divisor_classification=True,
  full35_branch_exhaustion_Delta_zero_non_normality_or_scheme_ideal_equality_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS independent actual normality-open minor, ZZ determinant and eps^59')

if __name__=='__main__':main()
