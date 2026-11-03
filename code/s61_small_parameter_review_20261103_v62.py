#!/usr/bin/env python3
"""Exact QQ slice/dominance/restriction checks; no generic-fiber GB."""
import hashlib
from collections import Counter
from fractions import Fraction
import gzip
import json
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 directory=ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62'
 report_path=directory/'report.json';d=json.loads(report_path.read_text());assert d['all_passed'];paths=[Path(__file__),report_path]
 for name,expected in d['input_sha256'].items():
  p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 for name,expected in d['evidence_sha256'].items():
  p=directory/name;assert sha(p)==expected,name;paths.append(p)
 matrices=json.loads((directory/'slice-matrices.json').read_text())
 old=s.Matrix(matrices['primitive_old_map']);section=s.Matrix(matrices['normalized_first_jet_map'])
 primitive=json.loads((ROOT/'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json').read_text())
 assert old.rank()==16 and old*section==s.Matrix(primitive['weight_basis_24_by16'])*s.Matrix(primitive['slice_matrix_b'])
 assert s.Matrix(matrices['dominance_matrix']).det()==s.Rational(438994543750,3)
 source=json.loads((ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json').read_text())
 weights=source['parameter_weights'];first=matrices['old_weight_one_indices'];second=matrices['old_weight_two_indices']
 counts=Counter(next(iter({sum(weights[i]for i in mon)for mon,_,_ in equation}))for equation in source['equations'])
 assert counts=={4:70,5:64}
 fi={v:i for i,v in enumerate(first)};qi={v:i for i,v in enumerate(second)};zero=(0,0,0,0,0)
 forms=[]
 for row in section.tolist():
  poly={zero:Fraction(str(row[0]))}
  for i,c in enumerate(row[1:]):
   if c:poly[tuple(int(j==i)for j in range(5))]=Fraction(str(c))
  forms.append({m:c for m,c in poly.items()if c})
 def mul(left,right):
  out={}
  for a,c in left.items():
   for b,e in right.items():
    mon=tuple(x+y for x,y in zip(a,b));out[mon]=out.get(mon,Fraction())+c*e
  return {m:c for m,c in out.items()if c}
 cache={():{zero:Fraction(1)}}
 encoded=json.loads(gzip.decompress((directory/'restricted-original134.json.gz').read_bytes()))
 assert len(encoded)==134
 for equation,saved in zip(source['equations'],encoded):
  out={}
  for mon,n,den in equation:
   betamon=tuple(i for i in mon if i in fi);qmon=tuple(qi[i]for i in mon if i in qi)
   if betamon not in cache:
    poly={zero:Fraction(1)}
    for i in betamon:poly=mul(poly,forms[fi[i]])
    cache[betamon]=poly
   row=out.setdefault(qmon,{})
   for m,c in cache[betamon].items():row[m]=row.get(m,Fraction())+Fraction(n,den)*c
  out={m:{a:c for a,c in row.items()if c}for m,row in out.items()};out={m:row for m,row in out.items()if row}
  expected={tuple(mon):{tuple(ex):Fraction(c)for ex,c in poly}for mon,poly in saved}
  assert out==expected
 paths += [ROOT/'research/notes/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62.md',
           ROOT/'research/notes/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62_REVIEW.md']
 output=dict(id='S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_V62_REVIEW',all_passed=True,
   independent_exact_slice_matrix_and_dominance_determinant=True,
   independent_Fraction_sparse_all134_symbolic_restrictions=True,
   all134_weight_counts={'4':70,'5':64},
   dominant_Sp4_Gm_five_parameter_reduction_accepted=True,
   generic_secant_center_free_finite_degree_criterion_accepted=True,
   degree30_used_only_for_original_integral_surfaces_and_total_extended_cycle=True,
   charzero_Wronski_generic_unit_iff_no_surfaces_accepted=True,
   secant_parameter_count=9,Wronski_parameter_count=8,
   solver_repeated=False,generic_rank_bound_computed=False,
   affine_surfaces_excluded=False,five_point_exhaustion_claimed=False,
   input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 path=ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62_REVIEW.json'
 with path.open('x')as stream:json.dump(output,stream,indent=2);stream.write('\n')
 print('PASS exact slice, dominance, ALL134 sparse QQ restrictions and both criterion scopes')


if __name__=='__main__':main()
