#!/usr/bin/env python3
"""Pin accepted exact-gcd peer and source consequences; no polynomial arithmetic rerun."""
from pathlib import Path
import argparse,hashlib,json,struct
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(kind,out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p,expected=None):
  h=sha(p)
  if expected:assert h==expected,(str(p),h)
  pins[str(p.relative_to(ROOT))]=h;return p
 def report(st,expected=None):return json.loads(pin(ROOT/'reports'/st/'report.json',expected).read_text())
 peer=report('S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63','6d7f852148000830f26d4f9ed709a7e7a5bc205e8c2678e4613648ea82b23bc7')
 assert peer['all_passed'] and peer['exact_char0_gcd_proved_up_to_QQ_scalar']
 assert not peer['irreducibility_proved'] and not peer['Fitting_or_global_rank_or_exhaustion_claimed']
 pin(ROOT/'code/s61_two_tail_CRT_ROOT_dense_review_20261105_v63.py',peer['code_sha256'])
 for p,h in peer['input_sha256'].items():pin(ROOT/p,h)
 for p,h in peer['evidence_sha256'].items():pin(ROOT/'reports/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63'/p,h)
 pin(ROOT/'research/notes/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63.md')
 proof=dict(accepted_independent_ALL_two_integer_identities=True,accepted_exact_QQ_gcd_up_to_scalar=True,
  no_new_polynomial_arithmetic=True,no_irreducibility_Fitting_global_rank_or_exhaustion_claim=True)
 if kind=='source':
  lm=report('S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_ROOT_REVIEW_20261105_v63','761e34bef478d9272c2d36f92d762a28a8288aaf284d41b7f0cdd3d1f97e70b6')
  assert lm['all_passed'] and lm['all_six_prime_and_nonassociate_forms_verified']
  glob=report('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63','f56e3d125e18ff1aa1c3d8720dc574062ff5a54759ec0bf13d17fc17267b7499')
  report('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW','6640c4f2be0e409a8f4854fd458fbbd700b4e28e7d6a360c153b4f2f67bcd25b')
  u=ROOT/'reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63/unit-inventory.json'
  inv=json.loads(pin(u,'75c20c39aa839ffb233d395b413cd59a3221f4a1de03e42a6b92e6be3b9b3d05').read_text())
  assert inv['factor_order']==['a','d','F','B','C','J']
  named_lm=[max(tuple(e) for e,c in p) for p in inv['named_factors']]
  assert all(e[1]==0 for e in named_lm)
  report('S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_SOURCE_20261105_v63','f0e6eadef71bc6ab661c5a5597bfc8c9f26527410b7f78d477d8127ba443a7de')
  pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_SOURCE_20261105_v63/N2-coprime-to-ALL-inverted-factors.json','b29d9f655fc01b3bf3db839e31a5c5e2238b7da6f1dd841c408ec135e976ca4c')
  union=report('S61_NONSQUARE_PREPOWER_TAIL1_FULL57_UNION_SOURCE_20261105_v63','3dbea3f3095fff06ffde0a809dd870d7370940d3567b149eccf2d650934b194c')
  assert union['all_passed'] and union['no_strict_cover_enlargement_or_bad_ideal_equality']
  report('S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW','c5ebf88f59f93ab7a4a4158ed10f50c89618fa9d1ffb9f8d3e7fab5ef38d5709')
  base=ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63'
  supports={}
  for n,magic,which,count in [('factor.bin',b'TTFAC001',2,94),('cofactor-N0.bin',b'TTCOF001',0,1795),('cofactor-N2.bin',b'TTCOF001',1,724)]:
   p=pin(base/n);support=[]
   with p.open('rb') as f:
    assert f.read(8)==magic and int.from_bytes(f.read(4),'little')==which
    assert int.from_bytes(f.read(4),'little')==count
    for _ in range(count):
     key=int.from_bytes(f.read(8),'little');sign=f.read(1);assert sign in (b'\0',b'\1')
     size=int.from_bytes(f.read(4),'little');assert 0<size<2048
     magnitude=f.read(size);assert len(magnitude)==size and magnitude[-1]
     k=key%(1<<26);e=(k>>20,(k>>10)&1023,k&1023);assert key//(1<<26)==2*e[0]+e[1]
     support.append(e)
    assert not f.read(1)
   supports[n]=dict(terms=count,lex_LM=max(support),coordinate_degrees=[max(e[j] for e in support) for j in range(3)])
  assert [supports[n]['lex_LM'] for n in supports]==[(6,1,4),(10,9,14),(6,10,10)]
  proof.update(complete_ring=inv['ring'],named_factor_lex_LM=named_lm,saved_support_only=supports,
   H_coprime_to_ALL_six='H divides N2; accepted N2 coprime to each inverted prime',
   H_A0_A2_nonunits='all have positive lex e exponent; every polynomial unit in this localization is QQstar times products of the six primes, with lex e exponent zero',
   exact_selected_tail_ideal='(Q0,Q2)=H*(A0,A2)',cofactor_ideal_if_proper_height_at_least_two=True,
   source_sufficient_open='D(D12*H) intersect (D(A0) union D(A2))',
   H_vs_D12_divisibility_UNRESOLVED=True,no_D12_leadingterm_assumption=True,
   common_H_divisor_NOT_actual_bad_Fitting_locus=True,no_new_factor_GCD_division_products_values_or_prime=True)
  pin(ROOT/'research/notes/S61_TWO_TAIL_COMMON_FACTOR_SOURCE_20261105_v63.md')
 pin(Path(__file__))
 target=out/'acceptance-and-source-proof.json';target.write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n')
 result=dict(id=out.name,all_passed=True,frozen=True,publication_and_source_only=True,kind=kind,
  independent_exact_gcd_peer_accepted=True,input_sha256=pins,evidence_sha256={target.name:sha(target)},
  no_new_polynomial_arithmetic_or_actual_bad_rank_claim=True)
 (out/'report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['peer','source'],required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.kind,a.output.resolve())
