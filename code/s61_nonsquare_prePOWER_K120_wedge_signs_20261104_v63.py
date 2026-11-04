#!/usr/bin/env python3
"""Independent finite formal-wedge proof; no parameter-polynomial arithmetic."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K120_WEDGE_SIGNS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 out.mkdir(parents=True,exist_ok=False)
 def mul(a,b):
  ans={}
  for m,c in a.items():
   for n,d in b.items():
    e=tuple(x+y for x,y in zip(m,n));ans[e]=ans.get(e,0)+c*d
  return {m:c for m,c in ans.items()if c}
 def add(a,b,s):
  ans=a.copy()
  for m,c in b.items():ans[m]=ans.get(m,0)+s*c
  return {m:c for m,c in ans.items()if c}
 # Formalweights order(g,h,h120); exteriorbasisorder(A,B,C).
 g={(1,0,0):1};h={(0,1,0):-1};hp={(0,0,1):-1};K1=[g,{},h];K2=[{},g,hp];actual={};discarded_diagonal=[]
 for i in range(3):
  for j in range(3):
   prod=mul(K1[i],K2[j])
   if i==j:
    if prod:discarded_diagonal.append(dict(exterior_index=i,formal_coeff=[[list(m),c]for m,c in prod.items()]))
   else:
    key=tuple(sorted([i,j]));actual[key]=add(actual.get(key,{}),prod,1 if i<j else-1)
 expected={(0,1):{(2,0,0):1},(0,2):{(1,0,1):-1},(1,2):{(1,1,0):1}}
 assert actual==expected
 coeffs={','.join(map(str,k)):[[list(m),c]for m,c in v.items()]for k,v in actual.items()}
 proof=dict(weight_order=['g','h','h120'],exterior_basis=['A','B','C'],computed_coefficients=coeffs,CC_diagonal_terms_discarded_by_exterior_identity=discarded_diagonal,
  identity='(gA-hC) wedge(gB-h120C)=g*(g A wedgeB-h120 A wedgeC-h C wedgeB)',
  all_coefficients_derived_from_the_two_literal_vectors=True,no_parameter_polynomial_coefficients=True)
 dump(out/'proof.json',proof);gate=ROOT/'reports/S61_NONSQUARE_PREPOWER_K120_T14_NESTED_SOURCE_PREFLIGHT_20261104_v63/report.json';note=ROOT/f'research/notes/{ID}.md'
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,nested_gate_report_sha256=sha(gate),formal_wedge_three_signs_verified=True,
  limited_check_clarification='The nestedgate compared literal alreadywrittenwedge coefficient dictionaries; this separatecheck derives them bilinearly fromK1,K2. All120 and5040 signchecks inthenestedgatewerealreadycomputational.',
  no_primary_file_changed=True,evidence_sha256={'proof.json':sha(out/'proof.json')},code_sha256=sha(Path(__file__)),note_sha256=sha(note))
 dump(out/'report.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
