#!/usr/bin/env python3
"""SOURCE ONLY: forced scalar normalization and next fixed-F Hasse budget."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T1_C4_GAUSS_F_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 q=load('S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_ENTRY_20261105_v63');assert q['all_passed']and q['result']['full_C4_multiplyback_zero']
 p=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_ENTRY_20261105_v63/quotient-integer-common-D.txt')
 with p.open()as f:
  assert f.readline().strip()=='T1C4QV001';Dq=int(f.readline());n=int(f.readline());assert n==4811
  L=0;mx=0;previous=None
  for _ in range(n):
   key,c=map(int,f.readline().split());assert c and(previous is None or previous<key);previous=key
   L+=abs(c);mx=max(mx,abs(c).bit_length())
  assert not f.read().strip()
 assert L.bit_length()==433 and mx==428
 p=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63/native-input.txt')
 with p.open()as f:f.readline();D1=int(f.readline())
 c=100**20;assert Dq==D1*c and D1.bit_length()==197 and c.bit_length()==133
 normalized_l1_bits=433-132;assert normalized_l1_bits==301
 owned=2*(4811+289+1254)+16000;assert owned==28708
 inner=4811*289+2*4811*33;scaled=2*4811*19+768+16+100+4812
 bits=301+16+(500**16).bit_length()+2;assert bits==463
 d=dict(source_only=True,normalization_coefficients_NOT_written_or_evaluated=True,
  quotient_peer_REQUIRED_before_any_actual_normalization_or_F_job=True,
  source_identity='C^4*Qint=100^20*I, I=(T1/a^3)_integer in ZZ[a,e,d]',
  C_integer_coefficients=[100,18,-9,-2],primitive_identity='-1*(-9)+4*(-2)=1',
  Gauss_proof='content(C^4)=1; hence content(C^4*Qint)=content(Qint)=100^20*content(I); every Qint coefficient is divisible by100^20',
  fixed_scalar_c='100^20',no_integer_content_search_or_polynomial_GCD=True,
  same_Q_representation='K=Qint/(100^20) in ZZ[a,e,d], Q=K/D1; D1 ORIGINAL197bits',
  future_exact_scalar_checks='for each existing integer coefficient x, assert x%c==0; y=x/c; assert y*c==x; and(Dq/c)*c==Dq. Never use two rational cross-products.',
  raw_Qint_coefficient_bits428=428,raw_Qint_l1_bits433=433,
  normalized_K_coefficient_bits_upper296=296,normalized_K_l1_bits_upper301=301,
  F='500*a+90*e^2+3*e-24',F_a0='(24-3*e-90*e^2)/500',F_root_coefficients_low_to_high=[24,-3,-90],F_root_l1=117,
  F_Hasse_formula='D1*500^16*[z^k]Q(a0+z,e,d)=sum K_ijl*choose(i,k)*500^(16-i+k)*(24-3e-90e^2)^(i-k)*e^j*d^l',
  orders='0..16 stopFIRSTnonzero; only F, no further unit or tail',
  source_height_proof='117<500, binom<=2^16, L1K<2^301 => all integer partials <2^463',
  hard_all_integer_bound518_UNCHANGED=518,normalization_backcheck_bits_upper429=429,
  Hasse_integer_bits_upper=bits,common_denominator_bits_upper342=342,
  input4811=4811,root_cache289=289,one_accumulator1254=1254,
  support_proof='e<=33+32=65,d<=18 =>66*19=1254',
  owned28708=owned,ownership='2*(inputcoeff/key +rootcoeff/index +accumulatorcoeff/key)+allmetadata16000',
  inner_products_upper=inner,ALL_scaled_root_and_polynomial_products_upper=inner+scaled,
  exact_known_scalar_divisibility_checks_upper4812=4812,
  RAW_upper8000000=8000000,soft45_hard60_RAM2GB=True,
  no_F_power_Hasse_normalization_or_other_coefficient_job_run=True,
  geometry='same exact Q polynomial, knownQQ scalar only; no bad-locus elimination or strict cover claim',
  new_parent_budget_REQUIRED=True)
 (out/'Gauss-and-F-preflight.json').write_text(json.dumps(d,indent=2)+'\n')
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,
  forced_scalar_normalization_by_Gauss_proved=True,normalized_coefficient_file_NOT_created=True,
  next_F_known_unit_job_ONLY_preflight=True,quotient_independent_peer_pending=True,
  future_hard_integer518=518,HH_bits463=463,owned28708=28708,ALL_products_upper=inner+scaled,RAW_upper8000000=8000000,
  input_sha256=pins,evidence_sha256={'Gauss-and-F-preflight.json':sha(out/'Gauss-and-F-preflight.json')},seconds=time.monotonic()-t)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
