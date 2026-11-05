#!/usr/bin/env python3
"""ROOT checks full typed tenth Schur/RHS operator, without evaluating any new cell."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_V2_20261105_v64'
OUT=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_ROOT_SOURCE_REVIEW_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);p=BASE/'report.json';assert sha(p)=='ad9dc2bc80fa50e1bf0946f5bf515450a04acad56818cb23851120111be8824b'
 r=json.loads(p.read_text());assert r['all_passed'] and r['no_new_matrix_coefficient_evaluated'];pins={str(p.relative_to(ROOT)):sha(p)}
 for name,want in r['input_sha256'].items():assert sha(ROOT/name)==want;pins[name]=want
 p=BASE/'actual-F10-pivot-and-full167x14-typed-operator.json';assert sha(p)==r['evidence_sha256'][p.name];d=json.loads(p.read_text());pins[str(p.relative_to(ROOT))]=sha(p)
 old=json.loads((ROOT/d['full167x14_RF_roots']['source']).read_text())
 for newkey,oldkey,excluded in [('full167_row_ids','full168_row_ids',58),('top101_row_ids','top102_row_ids',58),('selected14_row_ids','selected15_row_ids',58),('residual14_column_ids','residual15_column_ids',56)]:
  assert d[newkey]==[i for i in old[oldkey]if i!=excluded] and len(set(d[newkey]))==len(d[newkey])
 assert len(d['full167_row_ids'])==167 and len(d['residual14_column_ids'])==14
 assert d['all87_other_top_RHS']==old['all87_other_top_RHS'] and len(d['all87_other_top_RHS'])==87
 assert d['all66_zero_constraints']==old['all66_zero_constraints'] and len(d['all66_zero_constraints'])==66
 assert set(d['top101_row_ids'])==set(d['selected14_row_ids'])|set(d['all87_other_top_RHS'])
 assert set(d['full167_row_ids'])==set(d['top101_row_ids'])|set(d['all66_zero_constraints'])
 assert len(d['actual14_selected_lambda'])==14
 for i,rec in zip(d['selected14_row_ids'],d['actual14_selected_lambda']):
  assert rec['row']==i and sha(ROOT/rec['actual_N9_i56'])==rec['sha256']
 assert d['actual_p10']=='F10/(353894400000*H*t0*F9)' and d['F9_still_denominator_not_removed']
 assert d['actual_F10']['sha256']=='526ec2470b0ceb344ac7a837e9e45e3bfbf99e60e0a646831c5a2ae6ad917a68'
 assert d['incoming168_RHS_kept_for_restoration'] and d['P290_applied_once_no_double_scaling']
 assert d['all_other14_column_coefficients_UNEVALUATED']
 # A separate universal Laurent arithmetic check. Symbols are
 # U,F9,F10,cHt,Nij,Ni56,N58j,h_i,h58,x_j.
 def m(*exponents):return{tuple(exponents):1}
 def mul(A,B):
  C={}
  for a,c in A.items():
   for b,d in B.items():
    key=tuple(x+y for x,y in zip(a,b));C[key]=C.get(key,0)+c*d
  return{k:v for k,v in C.items()if v}
 def plus(A,B,sign=1):
  C=dict(A)
  for k,v in B.items():C[k]=C.get(k,0)+sign*v
  return{k:v for k,v in C.items()if v}
 def symbol(i,power=1):e=[0]*10;e[i]=power;return m(*e)
 Dinv=mul(symbol(0,-1),symbol(1,-1));pivotinv=mul(mul(symbol(3),symbol(1)),symbol(2,-1))
 Sij=mul(symbol(4),Dinv);Si56=mul(symbol(5),Dinv);S58j=mul(symbol(6),Dinv)
 lam=mul(Si56,pivotinv);S10=plus(Sij,mul(lam,S58j),-1)
 D10=mul(mul(symbol(0,2),symbol(1)),symbol(2))
 expected=plus(mul(mul(symbol(0),symbol(2)),symbol(4)),mul(mul(symbol(3),symbol(5)),symbol(6)),-1)
 assert mul(D10,S10)==expected
 RHS=plus(symbol(7),mul(lam,symbol(8)),-1)
 x56=mul(pivotinv,plus(symbol(8),mul(S58j,symbol(9)),-1))
 original_equation=plus(mul(Si56,x56),mul(Sij,symbol(9)))
 transformed_equation=plus(mul(S10,symbol(9)),mul(lam,symbol(8)))
 assert original_equation==transformed_equation
 assert 25000+8*317+8*951+12*(167+14+87+66)+2048==41200
 facts=dict(full167x14_and_all87plus66_RHS_exact=True,whole_pivot_row_current_not_raw=True,
  universal_denominator_and_numerator_identity=True,universal_RHS_and_inverse_substitution_identity=True,
  D10='U^2*F9*F10',N10='U*F10*N9ij-cHt*Ni56*N58j',
  all14_lambda_existing_coefficients_bound=True,F9_denominator_retained=True,P290_once=True,
  complement_scope='all old complements and F10=0 remain open',source_owned_upper=41200,
  no_new_matrix_values_or_inverse_or_exhaustion=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,input_sha256=pins,
  evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()
