"""Check the remaining d=J=0 coefficient in the six-quadric rank proof."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
from s61_nonsquare_universal_quad22_review_20261103_v63 import decode,linear_coefficient

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CONSTANT_QUAD6_20261104_v63'

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
 out=parser.parse_args().output.resolve();out.mkdir(parents=True,exist_ok=False)
 start=time.monotonic()
 source=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json'
 ell=decode(json.loads(source.read_text())['ell'])
 x0=linear_coefficient(ell,0);special={}
 for mon,c in x0.items():
  assert not any(mon[:7])
  if mon[9]:continue
  n=mon[7];special[n]=special.get(n,Fraction(0))+c*Fraction(-2,3)**mon[8]
 special={n:c for n,c in special.items()if c}
 assert special=={1:Fraction(3200,243)}
 paths=[Path(__file__),ROOT/'research/notes'/f'{ID}.md',source,
  ROOT/'code/s61_nonsquare_universal_quad22_review_20261103_v63.py',
  ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json',
  ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW/report.json']
 report=dict(id=ID,all_passed=True,d_zero_J_zero_forces_e='-2/3',
  ell_x0_on_d_zero_J_zero='3200*a/243',six_quadrics_pointwise_constant_rank=6,
  raw_open='a*C!=0',actual_weighted_chart='a*F*B*C!=0',
  weighted_degree2_quotient_dimension_on_actual_chart=23,
  kernel22_independence_extended_off_FULLU=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in ['id','all_passed','seconds']}))

if __name__=='__main__':main()
