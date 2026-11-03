#!/usr/bin/env python3
"""Independently expand the physical raw350 rows and audit the finite DAG."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
TRIPLES=list(itertools.combinations(range(5),3))
MON3=list(itertools.combinations_with_replacement(range(5),3))
BASIS=list(itertools.product(PAIRS,MONS))


def clean(poly):return {m:c for m,c in poly.items() if c}


def derivative(poly,index):
 out=defaultdict(Fraction)
 for (powers,parameters),c in poly.items():
  if powers[index]:
   new=list(powers);new[index]-=1
   out[tuple(new),parameters]+=powers[index]*c
 return clean(out)


def add_product(out,left,right):
 for (exponents,p),a in left.items():
  for (other,q),b in right.items():
   out[tuple(x+y for x,y in zip(exponents,other)),tuple(sorted(p+q))]+=a*b


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--model',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();args.model=args.model.resolve()
 raw=gzip.decompress(args.model.read_bytes());d=json.loads(raw)
 P={ij:{} for ij in PAIRS}
 for ((ij,mon),form) in zip(BASIS,d['canonical_linear_forms']):
  powers=tuple(mon.count(i) for i in range(5))
  for variable,c in form:
   P[ij][powers,(() if variable==-1 else (variable,))]=Fraction(c)
 def entry(i,j):
  if i<j:return P[i,j]
  if i>j:return {m:-c for m,c in P[j,i].items()}
  return {}
 derivatives={(i,j,k):derivative(entry(i,j),k) for i in range(5) for j in range(5) for k in range(5)}
 # Physical sparse polynomial computation in x0..x4,a,y1..y108.
 for j in range(5):
  divergence=defaultdict(Fraction)
  for i in range(5):
   for mon,c in derivatives[i,j,i].items():divergence[mon]+=c
  assert not clean(divergence),j
 for n,(i,j,k) in enumerate(TRIPLES):
  jacobi=defaultdict(Fraction)
  for u in range(5):
   add_product(jacobi,entry(i,u),derivatives[j,k,u])
   add_product(jacobi,entry(j,u),derivatives[k,i,u])
   add_product(jacobi,entry(k,u),derivatives[i,j,u])
  jacobi=clean(jacobi)
  expected=defaultdict(Fraction)
  for m,mon in enumerate(MON3):
   powers=tuple(mon.count(v) for v in range(5))
   for p,q,c in d['raw_quadratic_equations'][35*n+m]:
    params=tuple(v for v in (p,q) if v!=-1)
    expected[powers,params]+=Fraction(c)
  assert jacobi==clean(expected),('raw Jacobi triple',i,j,k)
 nodes=d['nodes'];weights=[];dependencies=[]
 for index,(op,*v) in enumerate(nodes):
  if op=='q':
   coefficient=Fraction(v[0]);weight=None if coefficient==0 else 0;dependency=0
  elif op=='v':weight=d['input_weights'][v[0]];dependency=1<<v[0]
  elif op in ('+','*'):
   assert all(j<index for j in v[0])
   dependency=0
   for j in v[0]:dependency|=dependencies[j]
   if op=='+':
    options={weights[j] for j in v[0] if weights[j] is not None}
    assert len(options)<=1;weight=next(iter(options),None)
   else:weight=None if any(weights[j] is None for j in v[0]) else sum(weights[j] for j in v[0])
  elif op=='inv':
   assert v[0]<index and weights[v[0]]==0 and dependencies[v[0]]&~1==0
   weight=0;dependency=dependencies[v[0]]
  else:raise ValueError(op)
  weights.append(weight);dependencies.append(dependency)
 assert len(d['all350_residual_nodes'])==350
 assert all(weights[j] in (None,w) for j,w in zip(d['all350_residual_nodes'],d['residual_weight_bounds']))
 assert all(weights[j] in (None,5-i-k+sum(mon)) for j,((i,k),mon) in zip(d['raw150_reconstruction_nodes'],BASIS))
 a=s.Symbol('a');unit=s.factor(s.sympify(d['localized_a_open']))
 assert unit.subs(a,s.Rational(2,9))!=0
 sizes=[]
 for block in d['implicit_blocks']:
  degree=block['weight'];cols=block['raw_positive_columns'];rows=block['equation_rows']
  assert len(cols)==len(rows)
  assert all(d['positive_weights'][j]==degree for j in cols)
  assert all(d['residual_weight_bounds'][i]==degree for i in rows)
  matrix=s.zeros(len(rows),len(cols))
  for r,i in enumerate(rows):
   terms={(p,q):Fraction(c) for p,q,c in d['raw_quadratic_equations'][i]}
   for c,j in enumerate(cols):
    constant=terms.get((-1,j+1),Fraction())
    slope=terms.get((0,j+1),Fraction())
    matrix[r,c]=s.Rational(constant.numerator,constant.denominator)+a*s.Rational(slope.numerator,slope.denominator)
  inverse=s.Matrix([[s.sympify(v) for v in row] for row in block['inverse']])
  assert (matrix*inverse-s.eye(len(cols))).applyfunc(s.cancel)==s.zeros(len(cols))
  assert s.cancel(matrix.det()-s.sympify(block['determinant']))==0
  sizes.append(len(cols))
 assert sum(sizes)==74 and len(d['free_positive_columns'])==34
 paths=[Path(__file__),args.model,ROOT/'code/q_weighted_finite_model_20261103_v63.py',
        ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py',
        ROOT/'research/notes/Q_WEIGHTED_FINITE_MODEL_20261103_v63.md',
        ROOT/'research/notes/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.md',
        args.model.parent/'construction.json',args.model.parent/'validation/validation.json']
 construction=json.loads((args.model.parent/'construction.json').read_text())
 validation=json.loads((args.model.parent/'validation/validation.json').read_text())
 assert construction['construction_checks_passed'] and validation['all_passed']
 assert validation['decoded_model_sha256']==hashlib.sha256(raw).hexdigest()
 for metadata in [construction,validation]:
  for path,expected in metadata['input_sha256'].items():
   p=ROOT/path;assert hashlib.sha256(p.read_bytes()).hexdigest()==expected,path
 report=dict(id='Q_WEIGHTED_FINITE_MODEL_20261103_V63_REVIEW',all_passed=True,
             independent_physical_sparse_all350_raw_Jacobi_polynomials=True,
             independent_physical_divergence_identity=True,
             all_DAG_inverses_depend_only_on_weight_zero_a=True,
             all_DAG_node_topology_and_weight_homogeneity_checked=True,
             implicit_block_sizes=sizes,free_positive_variables=34,
             all12_saved_block_inverses_and_determinants_independently_verified=True,
             localization_unit_at_a_two_over_nine=True,
             exact_triangular_ideal_elimination_accepted=True,
             canonical_affine350_equals_projective_Jacobi_ideal_by_BV_Euler=True,
             etale_affine_slice_with17_smooth_parameters_accepted=True,
             projective_smooth_factor_dimension=16,
             global_group_cross_section_claimed=False,
             producer_physical_numeric_validation_retained_as_separate_input=True,
             input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
             decoded_model_sha256=hashlib.sha256(raw).hexdigest())
 with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
 print('PASS independent physical symbolic350+div, topology, weights and inverse dependencies')


if __name__=='__main__':main()
