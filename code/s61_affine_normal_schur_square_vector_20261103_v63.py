#!/usr/bin/env python3
"""Sharp controlled Schur targets: delta^2 scalar or delta*unit vector."""
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
 out={};n=len(M)
 for p in itertools.permutations(range(n)):
  value={():(-1)**sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))}
  for i,j in enumerate(p):value=multiply(value,M[i][j])
  for mon,c in value.items():add(out,mon,c)
 return out
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
 deps=[ROOT/'reports/S61_AFFINE_NORMAL_SCHUR_CUBIC_POWER_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63_REVIEW/report.json']
 for dep in deps:
  data=json.loads(dep.read_text());assert data['all_passed']
  for rel,want in data['input_sha256'].items():assert sha(ROOT/rel)==want,rel
  for rel,want in data.get('evidence_sha256',{}).items():assert sha(dep.parent/rel)==want,rel
 # Universal five-row toy: Phi_std contains x_j*e_i for20variables.
 # The20independententriesV formfourconstantcombinationsofthose100columns.
 V=[[{(4*i+j,):1}for j in range(4)]for i in range(5)];lambdas=[]
 for i in range(5):lambdas.append({mon:(-1)**i*c for mon,c in determinant([row for j,row in enumerate(V)if j!=i]).items()})
 delta=lambdas[0]
 for j in range(4):
  result={}
  for i in range(5):
   for mon,c in multiply(lambdas[i],V[i][j]).items():add(result,mon,c)
  assert not result
 # FourtimesEuler avoids fractions: 4delta=sum_j x_j*d_j(delta).
 euler={}
 for mon,c in delta.items():
  for j in set(mon):
   reduced=list(mon);reduced.remove(j);add(euler,tuple(sorted((j,)+tuple(reduced))),mon.count(j)*c)
 assert euler=={mon:4*c for mon,c in delta.items()}
 square=multiply(delta,delta);scalar={}
 for mon,c in euler.items():
  for term,d in multiply({mon:c},delta).items():add(scalar,term,d)
 assert scalar=={mon:4*c for mon,c in square.items()}
 certificate=dict(universal_five_row_standard_module_toy=dict(parameter_variables=20,source_columns=100,source_image='m*R^5',pivot_columns=4,all4_lambdaV_identities_QQ_integer=True,delta_terms=len(delta),Euler_vector_and_square_membership_checked=True),
  conditional_kernel_identity='IfB=0,ker(Rprime(a)^5->E5)=m*Rprime(a)^5 because1,h,...,h4isKbasisE5.',
  exact_vector_criterion='B=0 iff delta*e0 inimagePhi',exact_scalar_criterion='B=0 iff delta^2 inJ_g',
  scalar_forward='Jg=lambda(imagePhi)=m*(lambda_0,...,lambda_4);delta=lambda_0 in m impliesdelta^2 inJg',
  vector_reverse='e0mapsunitA;delta*e0membershipimpliesdelta=0inA;deltaregularonB(a)forcesB=0',
  scalar_reverse='Jg subsetFitt0A subsetAnnA;deltaregularB(a)',
  formal_target_degrees=dict(vector=8,scalar=16),cofactor_degree_bound_not_claimed=True,actual_S61_membership_NOT_computed=True,old_cubic_power_criterion_remains_valid=True)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),*deps,ROOT/'research/notes/S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63',all_passed=True,universal5row_toy_lambda_and_Euler_verified=True,exact_delta_vector_and_square_scalar_exhaustion_criteria=True,no_actual_membership_or_unit_claim=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')})
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('DELTAVECTOR8/SQUARESCALAR16 EXACTCRITERIA PASS;MEMBERSHIPSUNCOMPUTED')
if __name__=='__main__':main()
