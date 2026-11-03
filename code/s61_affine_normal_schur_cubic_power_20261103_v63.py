#!/usr/bin/env python3
"""Exact universal Laplace identity and cubic-power exhaustion target."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(poly,mon,c):
 value=poly.get(mon,0)+c
 if value:poly[mon]=value
 else:poly.pop(mon,None)
def multiply(a,b):
 out={}
 for u,c in a.items():
  for v,d in b.items():add(out,tuple(sorted(u+v)),c*d)
 return out
def determinant(M):
 n=len(M);out={}
 for permutation in itertools.permutations(range(n)):
  sign=(-1)**sum(permutation[i]>permutation[j]for i in range(n)for j in range(i+1,n));product={():sign}
  for i,j in enumerate(permutation):product=multiply(product,M[i][j])
  for mon,c in product.items():add(out,mon,c)
 return out
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
 deps=[ROOT/'reports/S61_AFFINE_NORMAL_SCALAR_FITTING_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_SCALAR_SCHUR_20261103_v63_REVIEW/report.json',ROOT/'reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63_REVIEW/report.json']
 for dep in deps:
  data=json.loads(dep.read_text());assert data['all_passed']
  for rel,want in data['input_sha256'].items():assert sha(ROOT/rel)==want,rel
  for rel,want in data.get('evidence_sha256',{}).items():assert sha(dep.parent/rel)==want,rel
 var=lambda i:{(i,):1};B=[[var(4*i+j)for j in range(4)]for i in range(5)];X=[[var(20+5*i+j)for j in range(5)]for i in range(5)];delta=determinant(B[1:]);left=multiply(delta,determinant(X));right={}
 for j in range(5):
  g=determinant([B[i]+[X[i][j]]for i in range(5)]);minor=determinant([[X[i][k]for k in range(5)if k!=j]for i in range(1,5)])
  for mon,c in multiply(g,minor).items():add(right,mon,(-1)**j*c)
 assert left==right and len(delta)==24 and len(left)==2880
 data=dict(universal_identity='delta*det(X)=sum_j(-1)^j det[B4,Xcolumnj]*det(Xbottom4excludingcolumnj)',generic_integer_variables=45,delta_terms=len(delta),full_identity_nonzero_terms_each_side=len(left),all_QQ_integer_coefficients_checked=True,
  ideal_inclusion='delta*I5(Phi) subset J_g',cubic_power_criterion='B=0 iff delta^3 inJ_g overRprime_K(a)',
  forward_steps=['B=0 impliesFitt0A=m^5','delta in m^4 implies delta^2 in m^8 subsetm^5=I5Phi','universalidentityimpliesdelta^3 inJ_g'],
  reverse_steps=['J_g subsetI5Phi subsetAnnA','delta^3 inJ_g impliesdelta^3annihilatesB(a)','deltaregularonB(a) impliesB(a)=0;faithfulfieldextensionimpliesB=0'],
  target_formal_parameter_degree_atmost=24,cofactor_degree_bound_not_inferred=True,actual_delta_not_expanded=True,actual_cubic_membership_NOT_computed=True)
 (out/'certificate.json').write_text(json.dumps(data,indent=2)+'\n')
 paths=[Path(__file__),*deps,ROOT/'research/notes/S61_AFFINE_NORMAL_SCHUR_CUBIC_POWER_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_SCHUR_CUBIC_POWER_20261103_v63',all_passed=True,universal45variable_integer_Laplace_identity_verified=True,exact_cubic_power_exhaustion_criterion_proved=True,membership_or_unit_itself_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')})
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('UNIVERSAL LAPLACE/DELTA_CUBED_MEMBERSHIP_CRITERION PASS;MEMBERSHIPUNCOMPUTED')
if __name__=='__main__':main()
