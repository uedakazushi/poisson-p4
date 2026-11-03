#!/usr/bin/env python3
"""Independent polynomial-ring pullbacks of both closed A14 charts."""
import argparse
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 packed=folder/'certificate.json.gz';assert packed.read_bytes()==(args.replay/'certificate.json.gz').read_bytes();certificate=read(packed);paths=[Path(__file__),public,packed,ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 quadrics=read(ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz');cubics=read(ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz');names=quadrics['variables'];assert certificate['variables']==names
 indices={n:i for i,n in enumerate(names)};graphs=[]
 for branch in certificate['branches']:
  free=branch['free_variables'];assert len(free)==14;R,*v=ring(','.join(free),s.QQ);variables=dict(zip(free,v));graph=[]
  for terms in branch['all35_graph']:
   polynomial=R.zero
   for mon,c in terms:
    term=R(s.Rational(c))
    for i in mon:term*=variables[names[i]]
    polynomial+=term
   graph.append(polynomial)
  assert all(graph[indices[n]]==variables[n]for n in free)and graph[indices['z22']]==R.one
  assert graph[indices['z20']]==-2*graph[indices['z31']]
  assert graph[indices['z6']]==(R.zero if branch['branch']==0 else 5*variables['z27']/7)
  for row in quadrics['all350_quadratic_rows']:
   assert sum((R(s.Rational(c))*graph[i]*graph[j]for i,j,c in row),R.zero)==0
  for row in cubics['cubic_generators']:
   value=R.zero
   for mon,c in row['terms']:
    term=R(s.Rational(c))
    for i in mon:term*=graph[i]
    value+=term
   assert not value
  assert [p['variable']for p in branch['linear_eliminations']]==['z1','z4','z9','z3','z19']
  assert all(s.Rational(p['coefficient'])!=0 for p in branch['linear_eliminations'])
  graphs.append(graph)
 assert certificate['branches'][0]['free_variables']==certificate['branches'][1]['free_variables']
 paths += [ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63_REVIEW.md']
 report=dict(id='Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_entire_graph_archive_match=True,
  independent_QQ_14variable_polynomial_ring_all350Q_and17trueC_pullbacks_zero_both_charts=True,
  free14_coordinate_identity_and_z22_one_prove_closed_A14_charts=True,
  two_distinct_irreducible_homogeneous_cone_closures_dimension15=True,
  every12_homogeneous_linear_section_of_lower_ideal_J_dimension_at_least3=True,
  previous_J_only_12section_Artinian_target_impossible=True,
  full_initial_or_actual_local_dimension_lower15_or_upper12_refutation_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS both exact A14 charts and15D cones in lower ideal J only')

if __name__=='__main__':main()
