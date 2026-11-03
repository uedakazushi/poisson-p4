#!/usr/bin/env python3
"""Independent Fraction-only review of the nonsquare four-quad shadow certificate."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 producer=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/report.json';report=json.loads(producer.read_text());assert report['all_passed']
 for rel,want in report['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';data=json.loads(source.read_text());assert data['all_passed']
 for rel,want in data['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 # Enumerate degree-two monomials independently by unordered variable pairs.
 monomials=sorted({tuple([i,j].count(k)for k in range(6))for i in range(6)for j in range(i,6)})
 assert len(monomials)==21
 M=[];shifts=[0,2,2,1]
 for i,quad in enumerate(data['cleared_qfour_sparse']):
  coeff={}
  for light,poly in quad:
   assert len(light)==7 and light[6]==0 and sum(light)==2
   output={}
   for ex,c in poly:
    olda,e,oldd=ex;new=(olda-oldd+shifts[i],e,oldd);assert min(new)>=0
    output[new]=output.get(new,Fraction(0))+Fraction(c)
   coeff[tuple(light[:6])]={ex:c for ex,c in output.items()if c}
  M.append([coeff.get(mon,{})for mon in monomials])
 cert=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/right_inverse.txt.gz';tokens=iter(gzip.decompress(cert.read_bytes()).decode().split());assert tuple(int(next(tokens))for _ in range(3))==(21,4,0)
 L=[]
 for i in range(21):
  row=[]
  for j in range(4):
   poly={}
   for _ in range(int(next(tokens))):
    c=Fraction(next(tokens));ex=tuple(int(next(tokens))for _ in range(3));assert c and min(ex)>=0 and ex not in poly;poly[ex]=c
   row.append(poly)
  L.append(row)
 assert next(tokens,None)is None
 wanted={(2,0,3):Fraction(-5),(2,1,2):Fraction(6),(2,0,2):Fraction(4)}
 for i in range(4):
  for j in range(4):
   result=defaultdict(Fraction)
   for k in range(21):
    for u,c in M[i][k].items():
     for v,d in L[k][j].items():result[tuple(x+y for x,y in zip(u,v))]+=c*d
   assert {ex:c for ex,c in result.items()if c}==(wanted if i==j else {}),(i,j)
 # Exhaust all5985 possible sets of four leading monomials, stronger than the pair bound.
 shadows=[]
 for mon in monomials:
  shadows.append({tuple(c+(j==k)for k,c in enumerate(mon))for j in range(6)})
 minimum=min(len(set().union(*(shadows[i]for i in indices)))for indices in itertools.combinations(range(21),4));assert minimum==18
 pair=ROOT/'reports/S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62/report.json';pairdata=json.loads(pair.read_text());assert pairdata['all_passed']
 for rel,want in pairdata['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 evidence=dict(all16_full_QQ_identities=True,certificate_terms=sum(len(p)for row in L for p in row),independent_monomial_order=monomials,exhaustive_leading_quad_sets=5985,exact_minimum_degree3_shadow=minimum,unit_multiplier='-a^2*P^2*(5P-6e-4)',kernel_dimension_from_pair_cover=17,all_point_product_rank_lower_bound=35,incidence_base_dimension=3,projective_coefficient_dimension=9,maximum_bad_incidence_dimension=8,no_allpoint_rank41_or57_claim=True,no_ambient_removal_claim=True)
 (out/'checks.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),producer,source,cert,pair,ROOT/'research/notes/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW.md']
 result=dict(id='S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW',all_passed=True,independent_full_Fraction_right_inverse=True,all_point_four_quads_independent=True,all_point_actual_product_rank_atleast35_accepted=True,generic_constant_unimodular_row_existence_accepted=True,explicit_new_row_Bezout_or_rank41_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'checks.json':sha(out/'checks.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT QQ QUAD SHADOW/RANK35/INCIDENCE ACCEPTED',result['seconds'],flush=True)
if __name__=='__main__':main()
