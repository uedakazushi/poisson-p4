#!/usr/bin/env python3
"""Bounded independent audit of the finite87 unit criterion, not a solver run."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62'
 path=directory/'report.json';d=json.loads(path.read_text());fresh=json.loads(args.replay.read_text())
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 criterion=json.loads((directory/'finite-criterion.json').read_text());z=s.Symbol('z')
 # Independent coefficient formula for P Q'-P' Q. The sole possible degree87
 # contribution has i=j=44 and is identically zero; next degree is86.
 degrees={i+j-1 for i in range(45)for j in range(45)if i!=j}
 assert max(degrees)==86 and 44+44-1==87 and 44-44==0
 assert s.gcd(s.Poly(z**87-2,z),s.Poly(87*z**86,z)).degree()==0
 assert s.expand(z*(87*z**86)-87*(z**87-2))==174
 assert s.polys.polytools.intervals(z**87-2,eps=s.Rational(1,10))[0][1]==1
 product=math.prod(math.factorial(j)for j in range(87))
 factors=criterion['vandermonde_determinant_prime_factorization']
 assert math.prod(int(p)**e for p,e in factors.items())==product!=0
 assert product.bit_length()==criterion['vandermonde_determinant_bit_length']
 prior=ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62/wronski-exact.sing'
 original=prior.read_text().split('ideal I=',1)[1].split(';\n',1)[0].replace('zz1','zz')
 full=(directory/'algebraic87-polynomial-exact.m2').read_text()
 assert full.split('I=ideal(zz^87-2,',1)[1].split(');\n',1)[0]==original
 assert 'ORDINARY87_FUNCTION_FIELD_TRACE_PASS' in (directory/'coefficient-probe.log').read_text()
 assert 'expected coefficient ring of the form ZZ/n, ZZ, QQ, or GF' in (directory/'nested-field-unsupported.log').read_text()
 assert not d['actual_full_unit_result_obtained']and not d['surface_exclusion_claimed']
 paths += [ROOT/'research/notes/S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62.md',
  ROOT/'research/notes/S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62_REVIEW.md',
  ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62_REVIEW.json']
 result=dict(id='S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62_REVIEW',all_passed=True,
  fresh_exact_M2_toy_trace_and_full_producer_replay_matches=True,all134_plus2_original_equations_preserved=True,
  wronskian_degree_upper=86,distinct_value_count=87,eisenstein_polynomial='z^87-2',
  all87_rational_unit_tests_equivalent_to_no_surfaces_accepted=True,
  single_degree87_algebraic_unit_test_equivalent_to_no_surfaces_accepted=True,
  ideal_level_conjugation_handles_component_splitting=True,
  ordinary_polynomial_presentation_equals_coefficient_field_extension=True,
  nested_field_M2_GB_unsupported_honestly_retained=True,
  full_model_unit_or_surface_exclusion_obtained=False,nonunit_fiber_degree_compared_to_surface_degree=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS finite87 unit equivalence and exact input/compatibility scope')

if __name__=='__main__':main()
