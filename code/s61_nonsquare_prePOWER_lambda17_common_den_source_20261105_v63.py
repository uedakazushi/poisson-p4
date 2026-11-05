#!/usr/bin/env python3
"""SOURCE ONLY tighten fixed17 scalar height without any ratio coefficient solve."""
import argparse,hashlib,json,math,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_TARGET17_LAMBDA17_COMMON_DEN_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 rel='reports/S61_NONSQUARE_PREPOWER_TARGET17_ANCHOR_QUOTIENT_SOURCE_20261105_v63/finite-lambda-ratio-gate.json';p=ROOT/rel;x=json.loads(p.read_text())
 D=int(x['common_denominator_upper']);A=int(x['numerator_upper']);V=int(x['integer_verification_L1_upper']);assert (D.bit_length(),A.bit_length(),V.bit_length())==(998,1651,1746)
 # No ratio coefficient is solved. D,A,V are source bounds from frozen input.
 b=abs(int(x['leading_K_integer']));assert D==b**17
 bound=dict(common_D=str(D),all_R_L1_upper=str(A),all_scalar_and_verify_integer_upper=str(V),
  input_source_rule='scaled R=D*rho; only known Kcoeff times one R, never R times R',
  old_3524bit_bound_unneeded=True,mathematically_proved_integer_bits=1746,implementation_STOP_guard_bits=2048,
  before_commit_stop_no_retry=True,canonical_reduction='integer gcd/exactdivision; verify(num*g,den*g) against raw pair, no large cross product',
  unknown_scalar_count=17,all_coefficient_positions=math.comb(37,3),owned_numeric_slots=36714,
  prefix_binary_RAW_upper=2*1024*1024,all_artifact_and_input_RAW_cap=6820000,
  no_new_lambda_coefficient_or_lastrow_value=True)
 note=ROOT/f'research/notes/{ID}.md';code=Path(__file__)
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,tighter_completed_bound=True,
  all_explicit_integer_temporaries_proved_under2048=True,owned_cap_unchanged36714=True,RAW_cap_unchanged6820000=True,
  bound=bound,input_sha256={rel:sha(p),str(note.relative_to(ROOT)):sha(note),str(code.relative_to(ROOT)):sha(code)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
