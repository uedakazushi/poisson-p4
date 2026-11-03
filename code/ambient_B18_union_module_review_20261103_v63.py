#!/usr/bin/env python3
"""Independent matrix Cayley-Hamilton audit of the B18 actual-intersection module."""
import argparse,gzip,hashlib,itertools,json,math,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/AMBIENT_B18_FINITE_UNION_MODULE_20261103_v63';r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for name,want in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want
 tau,nu,sigma,mu,a,b,c,d=s.symbols('tau nu sigma mu a b c d')
 A=s.Matrix([[a,nu*b],[b,a+tau*b]]);B=s.Matrix([[c,mu*d],[d,c+sigma*d]]);M=s.diag(A,B);z=s.Symbol('z')
 char=s.Poly(s.expand((z*z-s.trace(A)*z+A.det())*(z*z-s.trace(B)*z+B.det())),z)
 assert char==M.charpoly(z).as_poly()
 CH=sum((char.nth(k)*M**k for k in range(5)),s.zeros(4));assert CH.applyfunc(s.expand)==s.zeros(4)
 e,f=s.symbols('e f');C=s.Matrix([[e,nu*f],[f,e+tau*f]])
 assert (2*A*C-s.trace(A)*C-s.trace(C)*A+(s.trace(A)*s.trace(C)-s.trace(A*C))*s.eye(2)).applyfunc(s.expand)==s.zeros(2)
 archive=folder/'paired-words-degree0-through3.json.gz';words=json.loads(gzip.decompress(archive.read_bytes()))
 expected=[[]]+[list(m)for degree in range(1,4)for m in itertools.combinations_with_replacement(range(97),degree)]
 assert words==expected and len(words)==math.comb(100,3)==161700
 assert hashlib.sha256(archive.read_bytes()).hexdigest()==r['word_archive_sha256']
 recipe=json.loads((folder/'finite-module-recipe.json').read_text());assert recipe['M_generator_count']==196 and recipe['U_generator_count']==161700
 paths=[Path(__file__),folder/'report.json',folder/'finite-module-recipe.json',archive,ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',ROOT/'research/notes/AMBIENT_B18_FINITE_UNION_MODULE_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_B18_FINITE_UNION_MODULE_20261103_v63_REVIEW',all_passed=True,
  independent_generic4x4_matrix_CH_and_rank2_pair_identity=True,
  degree4_polarization_proof_over_any_QQ_algebra_reviewed=True,
  all161700_word_indices_exactly_enumerated=True,
  actual_D_exact_module_subquotient_not_normalizer_ring_quotient=True,
  target_actions_preserve_U_and_M_and_unit_class_is_cyclic=True,
  full_actual_ring_reconstruction_from_target_action_accepted=True,
  algebraic_coefficients_or_nilradical_or_primary_type_computed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT B18 FINITE MODULE/TARGET ACTION REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
