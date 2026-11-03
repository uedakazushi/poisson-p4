#!/usr/bin/env python3
"""Verify the fourth auxiliary five-quadric pivot; no Groebner basis search.

Extract the report archive at repository root and source cloud activation first.
The third core native variable d denotes P=ad. This is not main pivot 35.
"""
from pathlib import Path
import hashlib,json,subprocess,shlex,sys,time
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/five-fourth-unit-replay-v60'
OUT.mkdir(parents=True,exist_ok=True);t=time.monotonic()
W=ROOT/'work/nonsquare-v62/five-quadric-degree23-v60'
V=ROOT/'work/nonsquare-v62/five-rank-last-row-v61'
report=json.loads((ROOT/'reports/S61_NONSQUARE_FIVE_QUADRIC_FOURTH_UNIT_20261103_v63/report.json').read_text())
for path,sha in report['input_sha256'].items():
 assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
def build(stem):
 exe=OUT/stem
 subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'code'/(stem+'.cpp')),'-o',str(exe),*flags],check=True)
 return exe
check=build('polynomial_exact_scalar_unit_check_20261103_v62')
units=ROOT/'work/nonsquare-v62/raw-fourquad-structure-v49/aPJ-units.txt'
with(OUT/'scalar-identity.log').open('w')as f:
 subprocess.run([str(check),str(W/'source-row0.txt'),str(W/'D2-row0.witness'),str(units)],stdout=f,check=True)
text=(OUT/'scalar-identity.log').read_text();assert 'EXPONENTS 3 3 2' in text and 'SCALAR_UNIT_IDENTITY_PASS' in text
pivot=build('polynomial_unit_span_pivots_20261103_v62')
with(OUT/'pivot.log').open('w')as f:
 subprocess.run([str(pivot),str(ROOT/'work/nonsquare-v62/five-quadric-minors-v50/after-primitive.txt'),str(OUT/'last-row.txt'),str(units),'1',str(W/'D2-row0.witness'),'0'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
assert hashlib.sha256((OUT/'last-row.txt').read_bytes()).hexdigest()==hashlib.sha256((W/'last-row.txt').read_bytes()).hexdigest()
primitive=build('s61_nonsquare_five_quadric_scalar_primitive_20261103_v63')
with(OUT/'scalings.json').open('w')as f,(OUT/'primitive.log').open('w')as err:
 subprocess.run([str(primitive),str(OUT/'last-row.txt'),str(units),str(OUT/'primitive-row.txt')],stdout=f,stderr=err,check=True)
assert json.loads((OUT/'scalings.json').read_text())==json.loads((V/'column-scalings.json').read_text())
assert hashlib.sha256((OUT/'primitive-row.txt').read_bytes()).hexdigest()==hashlib.sha256((V/'primitive-row.txt').read_bytes()).hexdigest()
r=dict(all_exact_replay_gates_passed=True,auxiliary_five_quadric_pivots=4,remaining_shape=[1,18],known_unit='a^3*P^3*J^2',only_known_units=['a','P','J'],full_rank5_not_proved=True,main_pivot34_count_unchanged=True,raw_rank24_or_product41_or_ambient_not_proved=True,seconds=time.monotonic()-t)
(OUT/'replay-report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
