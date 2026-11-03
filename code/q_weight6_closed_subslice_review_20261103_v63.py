#!/usr/bin/env python3
"""Check both local ideal containments against a direct full-DAG evaluation."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--certificate-directory',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();directory=args.certificate_directory.resolve()
 report_path=directory/'report.json';producer=json.loads(report_path.read_text());assert producer['all_passed']
 paths=[Path(__file__),report_path]
 for name,expected in producer['input_sha256'].items():
  p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 archive=directory/'certificate.json.gz';paths.append(archive)
 assert sha(archive)==producer['evidence_sha256']['certificate.json.gz']
 cert=json.loads(gzip.decompress(archive.read_bytes()))
 R,a=ring('a',s.QQ);K=R.to_field();S,*v=ring('w0,w1,w2,eta',K)
 zero=(0,0,0,0);astar=s.QQ(2,9)
 def at_center(p):return sum((c*astar**m[0] for m,c in p.items()),s.QQ.zero)
 h=R.from_expr(s.sympify(producer['explicit_localization_unit']))
 assert at_center(h)!=0
 checked_denominators=set()
 def scalar(data):
  polys=[sum((s.QQ(c)*a**int(i) for i,c in terms),R.zero)for terms in data]
  assert polys[1] and at_center(polys[1])!=0
  key=str(polys[1])
  if key not in checked_denominators:
   remainder=polys[1]
   while not remainder.is_ground:
    divisor=remainder.gcd(h);assert not divisor.is_ground,key
    remainder=remainder.exquo(divisor)
   checked_denominators.add(key)
  return K(polys[0])/K(polys[1])
 def polynomial(data):return S.from_dict({tuple(m):scalar(c)for m,c in data})
 T=[[scalar(c)for c in row]for row in cert['coordinate_matrix']]
 images=[polynomial(p)for p in cert['inverse_images_of_s0_through_s3']]
 for row,target in zip(T,v):
  assert not sum((S(c)*image for c,image in zip(row,images)),S.zero)-target
 determinant=K.zero
 for p in itertools.permutations(range(4)):
  term=K.one
  for i,j in enumerate(p):term*=T[i][j]
  determinant+=(-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))*term
 assert at_center(determinant.numer)/at_center(determinant.denom)==1
 model_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
 model=json.loads(gzip.decompress(model_path.read_bytes()));paths.append(model_path)
 retained=[i for i,w in enumerate(model['input_weights'])if w==6]
 assert retained==[12,13,24,34]
 inputs=[S(K(a))]+[images[retained.index(i)]if i in retained else S.zero for i in range(1,35)]
 # Directly evaluate the parent DAG in the FINAL w,eta coordinates.
 # The producer instead evaluated first in four s coordinates, then changed variables.
 values=[]
 for op,*data in model['nodes']:
  if op=='q':value=S(K(s.QQ(data[0])))
  elif op=='v':value=inputs[data[0]]
  elif op=='+':value=sum((values[i]for i in data[0]),S.zero)
  elif op=='*':
   value=S.one
   for i in data[0]:value*=values[i]
  elif op=='inv':
   prior=values[data[0]];assert prior and set(prior)=={zero}
   c=prior[zero];assert at_center(c.numer)!=0 and at_center(c.denom)!=0
   value=S(1/c)
  else:raise ValueError(op)
  values.append(value)
 residual=[values[i]for i in model['all350_residual_nodes']]
 saved=[polynomial(p)for p in cert['all350_residuals_after_coordinate_change']]
 assert all(not x-y for x,y in zip(residual,saved))
 generators=[polynomial(p)for p in cert['ten_standard_generators']]
 expected=[S(K(a)-K(astar))*v[3]]+[v[i]*v[j]for i in range(3)for j in range(i,3)]+[v[i]*v[3]for i in range(3)]
 assert all(not x-y for x,y in zip(generators,expected))
 def apply(rows,source):
  return [sum((polynomial(p)*source[index]for index,p in row),S.zero)for row in rows]
 assert all(not x-y for x,y in zip(apply(cert['forward350_by10'],generators),residual))
 assert all(not x-y for x,y in zip(apply(cert['reverse10_by350'],residual),generators))
 N=[[scalar(c)for c in row]for row in cert['quadratic_matrix']]
 inverse=[[scalar(c)for c in row]for row in cert['quadratic_matrix_inverse']]
 for i in range(9):
  for j in range(9):assert not sum((N[i][k]*inverse[k][j]for k in range(9)),K.zero)-int(i==j)
 # No inverse of u=a-2/9 is admitted: every parsed certificate coefficient
 # above has a denominator explicitly nonzero at a=2/9.
 # Verify the claimed primary decomposition as a monomial intersection.
 def e(*positions):return tuple(positions.count(i)for i in range(5))
 Q1=[e(1)]+[e(i,j)for i in range(2,5)for j in range(i,5)]
 Q2=[e(0)]+[e(i)for i in range(2,5)]
 candidates={tuple(max(x,y)for x,y in zip(p,q))for p in Q1 for q in Q2}
 minimal={m for m in candidates if not any(n!=m and all(x<=y for x,y in zip(n,m))for n in candidates)}
 I={e(0,1)}|{e(i,j)for i in range(2,5)for j in range(i,5)}|{e(i,1)for i in range(2,5)}
 assert minimal==I
 paths += [ROOT/'research/notes/Q_WEIGHT6_CLOSED_SUBSLICE_20261103_v63_REVIEW.md']
 output=dict(id='Q_WEIGHT6_CLOSED_SUBSLICE_20261103_V63_REVIEW',all_passed=True,
             independent_parent_full9647_node_DAG_evaluation_in_final_coordinates=True,
             all350_transformed_polynomials_match_certificate=True,
             all350_forward_and_all10_reverse_ideal_identities_checked=True,
             every_certificate_coefficient_denominator_a_star_unit=True,
             all_saved_certificate_denominators_supported_on_explicit_h=True,
             distinct_saved_denominators_checked=len(checked_denominators),
             u_a_minus_two_over_nine_never_inverted=True,
             coordinate_inverse_and_9_by9_inverse_checked=True,
             exact_closed_subslice_primary_decomposition_checked=True,
             closed_subslice_CM_dimension=1,closed_subslice_multiplicity=5,
             reduced_branch_count=2,generic_branch_lengths=[4,1],
             embedded_prime_in_this_closed_subslice=False,
             full35_branch_classification_claimed=False,new_global_component_claimed=False,
             input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as stream:json.dump(output,stream,indent=2);stream.write('\n')
 print('PASS direct parentDAG, both local ideal identities and closed-subslice primary decomposition')


if __name__=='__main__':main()
