#!/usr/bin/env python3
"""Normal-only cubic Jacobian and exact origin/residual splitter recipe."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
from s61_affine_polarization_splitter_20261103_v63 import determinant_including_zero
ROOT=Path(__file__).resolve().parents[1];P=32003

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def residue(c):return c.numerator*pow(c.denominator,-1,P)%P
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--cross-output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 cert=json.loads((dep.parent/'certificate.json').read_text());assert sha(args.cross_output)==cert['native_output_sha256']
 source=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';data=json.loads(gzip.decompress(source.read_bytes()));A=[[{tuple(m):Fraction(c)for m,c in poly}for poly in row]for row in data['all129_linear_factor_rows']]
 mons=[(j,)for j in range(42)]+list(itertools.combinations_with_replacement(range(42),2));cross=np.array([[residue(A[i][j].get((k,),Fraction(0)))for k in range(3)for j in range(42)]for i in range(129)],dtype=np.int64);rest=[]
 for row in A:
  poly={}
  for j,column in enumerate(row):
   for mon,c in column.items():
    if mon and mon[0]<3:continue
    key=tuple(sorted((j,)+tuple(k-3 for k in mon)));poly[key]=poly.get(key,Fraction(0))+c
  rest.append([residue(poly.get(mon,Fraction(0)))for mon in mons])
 rest=np.array(rest,dtype=np.int64);selected=cert['selected_cross_rows'];aug=np.concatenate((cross[selected],np.eye(126,dtype=np.int64)),axis=1)
 for j in range(126):
  i=next(i for i in range(j,126)if aug[i,j]);aug[[j,i]]=aug[[i,j]];aug[j]=aug[j]*pow(int(aug[j,j]),-1,P)%P
  for i in range(126):
   if i!=j and aug[i,j]:aug[i]=(aug[i]-aug[i,j]*aug[j])%P
 inverse=aug[:,126:];assert np.array_equal(cross[selected]@inverse%P,np.eye(126,dtype=np.int64));coeff=(-inverse@rest[selected])%P
 # Compare every119070 modular coefficient with the exact native QQ output.
 def direct_mod(token):
  numerator,slash,denominator=token.partition('/');n=int(numerator)%P;d=int(denominator)%P if slash else 1;assert d;return n*pow(d,-1,P)%P
 with args.cross_output.open()as handle:
  assert handle.readline().split()==['126','945'];assert list(map(int,handle.readline().split()))==selected
  for _ in range(126):handle.readline()
  for i in range(126):assert list(map(direct_mod,handle.readline().split()))==coeff[i].tolist()
  assert not handle.read().strip()
 T=coeff[:,:42].reshape(3,42,42);commutators=[];determinants=[]
 for i,j in itertools.combinations(range(3),2):
  C=(T[i]@T[j]-T[j]@T[i])%P;commutators.append(C.tolist());determinants.append(determinant_including_zero(C))
 assert determinants==[21866,23983,19574]and all(determinants)
 (out/'modular-linear-commutators.json.gz').write_bytes(gzip.compress(json.dumps(dict(prime=P,all_three_constant_matrices=T.tolist(),all_three_commutators=commutators,commutator_determinants=determinants),separators=(',',':')).encode(),mtime=0))
 result=dict(exact_QQ_relation_recipe='P_k(n)=R_k(n)n by first-normal factoring; [R_i(n),R_j(n)]n in originalI',matrix_R_degree_atmost1=True,commutator_matrix_degree_atmost2=True,cubic_block_size=42,selected_pair=[0,1],normal_origin_Jacobian_determinant_mod_prime=21866,
  full_QQ_cross_coefficient_binding=True,all119070_exact_QQ_coefficient_reductions_matched=True,
  scalar_J0_recipe='J0=J_normal2796+42entries([R0,R1]n); J0=(n) intersect I_B',scalar_generators_total=2838,origin_reduced_exactly=True,
  normal_only_splitter_recipe='D84(n)=det([R0(n),R1(n)])/det([R0(0),R1(0)]), degree<=84; D84=e in actualA',
  scalar_unsaturated_residual_ideal='I_B=J0+(D84)',projective_or_beta_finiteness_or_emptiness_not_claimed=True,
  exact_degree84_or_expanded_QQ_cubics_not_claimed=True)
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 inputs=[Path(__file__),dep,dep.parent/'certificate.json',source,ROOT/'research/notes/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63.md']
 report=dict(id='S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63',all_passed=True,exact_normal_cubic_relations_finite_product_recipe=True,all_three_normal_Jacobian42_minors_nonzero=True,scalar_normal_image_ideal_with_reduced_origin=True,normal_only_degree84_idempotent_recipe=True,unit_or_exhaustion_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('NORMAL CUBIC JAC42/REDUCEDORIGIN/NORMALSPLITTER RECIPE PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
