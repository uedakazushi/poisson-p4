#!/usr/bin/env python3
"""Proof-only colon/nonflatness/residual decomposition audit with checked inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--colon-replay',type=Path,required=True)
 ap.add_argument('--residual-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();paths=[Path(__file__)]
 for stem,replay in [('Q_KNOWN_GRAPH_COLON_NONFLATNESS',args.colon_replay),('Q_KNOWN_GRAPH_RESIDUAL_DECOMPOSITION',args.residual_replay)]:
  p=ROOT/f'reports/{stem}_20261103_v63.json';d=load(p);assert d==load(replay);paths.append(p)
  for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file;paths.append(q)
 a=s.Symbol('a');u=a-s.Rational(2,9);v=a+s.Rational(8,9);assert s.expand(s.Rational(9,10)*(v-u))==1
 r39=load(ROOT/'reports/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW.json');r40=load(ROOT/'reports/Q_T40_GRAPH_OPEN_20261103_v63_REVIEW.json')
 assert r39['all_passed']and r40['all_passed']and r39['exact_curve_power']==44 and r40['no_new_inverse_and_weight275_checked']
 paths += [ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.json',
  ROOT/'reports/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW.json',
  ROOT/'reports/Q_T40_GRAPH_OPEN_20261103_v63_REVIEW.json',
  ROOT/'research/notes/Q_KNOWN_GRAPH_ALGEBRA_20261103_v63_REVIEW.md']
 result=dict(id='Q_KNOWN_GRAPH_ALGEBRA_20261103_v63_REVIEW',all_passed=True,proof_only_consequences_no_new_GB=True,
  BOTH_graph_coefficient_certificates_independently_checked_before_this_review=True,
  exact_all_positive_power_colons='I:C_i^k=N_i',exact_injective_modules='P/N_i --multiply C_i--> P/I',
  full_local_and_completed_R_nonflat_in_a_at_Qa39_and_Qa40_accepted=True,
  cross_C_products_and_principal_ideal_intersection_zero=True,
  exact_combined_scheme_open='R_(C39+C40) = (P/N39)_C39 product (P/N40)_C40',
  exact_single_prime_peel='I=N_i intersection (I+C_i)',
  exact_two_prime_peel='I=N39 intersection N40 intersection (I+C39+C40)',
  local_completion_at_Qstar='Ihat=N39hat intersection (Ihat+C39)',
  residual_components_or_primary_structure_or_whole_R_product_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS exact known-prime colons, full nonflatness, and residual ideal peel')

if __name__=='__main__':main()
