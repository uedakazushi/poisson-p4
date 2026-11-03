#!/usr/bin/env python3
"""Independent polynomial-ring cofactor-chain and T39 tangent audit."""
import argparse
from fractions import Fraction as F
import gzip
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_TRUE_INITIAL_ETA_OPEN_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json')
 assert{k:v for k,v in d.items()if k not in('seconds','certificate_sha256')}=={k:v for k,v in f.items()if k not in('seconds','certificate_sha256')}
 packed=folder/'certificate.json.gz';decoded=gzip.decompress(packed.read_bytes());assert decoded==gzip.decompress((args.replay/'certificate.json.gz').read_bytes());c=json.loads(decoded)
 paths=[Path(__file__),public,packed,ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 q=read(ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz');g=read(ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz');names=q['variables'];assert c['variables']==names
 R,*x=ring(','.join(names),s.QQ);ix={n:i for i,n in enumerate(names)}
 def poly(terms):
  value=R.zero
  for mon,v in terms:
   p=R(s.Rational(v))
   for i in mon:p*=x[i]
   value+=p
  return value
 quadrics=[poly([[[i,j],v]for i,j,v in row])for row in q['all350_quadratic_rows']]
 cubics=[poly(row['terms'])for row in g['cubic_generators']];rows=[quadrics[i]for i in g['quadratic_generators']]+cubics
 for record in c['constant_row_certificates']:
  expected=sum((s.Rational(v)*quadrics[rawindex]for i,rawindex,v in record['constant_cofactors']),R.zero)
  assert expected==poly(record['polynomial'])
 eta=x[ix['z33']]
 for name,v in [('z11',4),('z12',6),('z23',-4)]:assert eta*(x[ix[name]]-v*eta)in[poly(r['polynomial'])for r in c['constant_row_certificates']]
 mons=list(itertools.combinations_with_replacement(range(35),2));mi={m:i for i,m in enumerate(mons)};M=s.zeros(57,630)
 for i,p in enumerate([q['all350_quadratic_rows'][j]for j in g['quadratic_generators']]):
  for a,b,v in p:M[i,mi[a,b]]=s.Rational(v)
 pure={ix[n]for n in['z11','z12','z23','z33']};outside=[j for j,m in enumerate(mons)if not set(m)<=pure]
 assert len(DomainMatrix.from_Matrix(M[:,outside]).convert_to(s.QQ).rref()[1])==48
 for name,v in [('z33',1),('z11',4),('z12',6),('z23',-4)]:rows=[p.compose(x[ix[name]],R(v))for p in rows]
 assert rows==[poly(p)for p in c['prepared_74_rows']]
 assert len(c['triangular_chain'])==23
 for step in c['triangular_chain']:
  i=step['index'];pivot=s.Rational(step['pivot_coefficient']);rhs=poly(step['rhs']);selected=poly(step['selected_polynomial'])
  assert pivot and selected==rows[step['source_generator']]==pivot*(x[i]-rhs)
  assert all(not mon[i]for mon in rhs)
  cof={k:poly(p)for k,p in step['division_cofactors']};after=[p.compose(x[i],rhs)for p in rows]
  for k,(p,new)in enumerate(zip(rows,after)):assert p-new==(x[i]-rhs)*cof.get(k,R.zero)
  rows=after;assert sum(bool(p)for p in rows)==step['nonzero_rows_after']
 assert not any(rows)
 graph=[poly(c['homogeneous_graph'][n])for n in names];hfree=[ix[n]for n in c['free_variables']]+[ix['z33']]
 assert len(hfree)==9 and all(graph[i]==x[i]for i in hfree)
 for p in quadrics+cubics:
  value=R.zero
  for mon,v in p.items():
   term=R(v)
   for i,power in enumerate(mon):
    if power:term*=graph[i]**power
   value+=term
  assert not value
 equations=[poly(p)for p in c['homogeneous_graph_26_equations']]
 assert len(equations)==26 and all(x[i]-graph[i]in equations for i in range(35)if i not in hfree)
 # Fresh first-order evaluation of the previously certified ACTUAL A9 map.
 oldfolder=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63';old=read(oldfolder/'certificate.json');tm=read(oldfolder/'T39-nine-parameter-model.json.gz');values=[]
 for op,arg in tm['nodes']:
  if op=='q':v=[F(arg)]+[F(0)]*9
  elif op=='v':v=[F(0)]+[F(i==arg)for i in range(9)]
  elif op=='+':v=[sum((values[j][i]for j in arg),F(0))for i in range(10)]
  elif op=='*':
   v=[F(1)]+[F(0)]*9
   for j in arg:
    right=values[j];v=[v[0]*right[0]]+[v[0]*right[i]+right[0]*v[i]for i in range(1,10)]
  else:raise ValueError(op)
  values.append(v)
 parent=[list(values[i])for i in tm['parent35_input_nodes']];parent[0][0]-=F(2,9);assert all(not row[0]for row in parent)
 A=s.Matrix([[s.Rational(v.numerator,v.denominator)for v in row[1:]]for row in parent]);graphmatrix=s.Matrix([[p.get(tuple(int(k==j)for k in range(35)),0)for j in hfree]for p in graph])
 assert A==graphmatrix and A[hfree,:]==s.eye(9)
 H=s.Matrix([[s.Rational(v)for v in row]for row in old['Q51_jacobian_45x35']]);assert H.rank()==26 and H*A==s.zeros(45,9)
 paths += [ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_TRUE_INITIAL_ETA_OPEN_20261103_v63.md',ROOT/'research/notes/Q_TRUE_INITIAL_ETA_OPEN_20261103_v63_REVIEW.md']
 result=dict(id='Q_TRUE_INITIAL_ETA_OPEN_20261103_v63_REVIEW',all_passed=True,
  fresh_decoded_certificate_and_report_match_compressed_time_allowed=True,
  independent_all12_QQ_constant_row_certificates_and_pure4_rowspace9=True,
  independent23_steps_all74_polynomial_division_cofactor_identities=True,
  independent_reverse_all350Q_and17C_linear_graph_pullbacks_zero=True,
  fresh_actual_T39_first_order_all35_by9_matches_linear_graph=True,
  original_J_eta_equals_linear_prime_P_eta_scheme_theoretically=True,
  full_true_ordinary_initial_eta_equals_P_eta_by_J_in_initialI_in_T39_tangentP=True,
  actual_tangent_cone_eta_open='reduced smooth Gm x A8',all_other_cone_components_eta_zero=True,
  actual_completed_I_eta_linear_equality_or_global_cone_dimension_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS exact J_eta=P_eta and full true ordinary initial eta open by actual T39 sandwich')

if __name__=='__main__':main()
