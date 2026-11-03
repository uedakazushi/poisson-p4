#!/usr/bin/env python3
"""Certify the reduced support of the full old-quadratic central S over F_p.

This includes rho.  It is a special-fiber support computation, not an ambient
nilpotence assertion and not a characteristic-zero radical calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

from s61_marked_sigma_old_plane_central_20261103_v62 import program
from s61_marked_sigma_pencil_rank_20261103 import PRIME, ROOT, read_rows


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--root', type=int, choices=(3083, 28920), default=3083)
    parser.add_argument('--timeout', type=int, default=180)
    args = parser.parse_args()
    expected = json.loads((ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json').read_text())['row_output_sha256']
    assert digest(args.rows) == expected
    assert (args.root**2+2) % PRIME == 0
    args.output.mkdir(parents=True, exist_ok=False)
    text = program(read_rows(args.rows), args.root, include_rho=True)
    text = text[:text.index('print "INPUT READY";')]
    text += f'''print "INPUT READY"; flush stdio;
G=gb I; print "GB READY"; flush stdio;
print("HS="|toString hilbertSeries(R/I,Reduce=>true));
print("GB_GENERATORS="|toString numColumns gens G); flush stdio;
L=ideal(v0,v1,v2,v3,v4,v5,v6,v7,v8,v9-{args.root}*v10,v11,v13*v10,v13*v12);
assert(isSubset(I,L)); print "QUOTIENT CHECKED"; flush stdio;
scan(flatten entries gens L,f->(n:=1; while n<=8 and f^n%G!=0 do n=n+1;
  print("NIL_WITNESS="|toString f|":"|toString n); flush stdio; assert(n<=8)));
print "FULL SUPPORT AXIS UNION PLANE PASS"; flush stdio; exit 0;
'''
    script = args.output/'support.m2'
    script.write_text(text)
    log = args.output/'support.log'
    started = time.monotonic()
    timed_out = False
    with log.open('w') as stream:
        try:
            result = subprocess.run([shutil.which('M2'), '--script', str(script.resolve())],
                                    stdout=stream, stderr=subprocess.STDOUT,
                                    timeout=args.timeout)
            returncode = result.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            returncode = None
    contents = log.read_text()
    passed = returncode == 0 and 'FULL SUPPORT AXIS UNION PLANE PASS' in contents
    sources = [Path(__file__), ROOT/'code/s61_marked_sigma_old_plane_central_20261103_v62.py',
               ROOT/'code/s61_marked_sigma_pencil_rank_20261103.py']
    report = dict(id='S61_MARKED_SIGMA_OLD_FULL_SUPPORT_20261103_V62',
                  field=f'GF({PRIME})', root=args.root, includes_rho=True,
                  execution_completed=returncode == 0, returncode=returncode,
                  timed_out=timed_out, special_fiber_reduced_support_proved=passed,
                  reduced_support='rho-axis union weighted (2,3) plane',
                  ideal_generators=['v0', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8',
                                    f'v9-{args.root}*v10', 'v11', 'rho*v10', 'rho*v12'],
                  characteristic_zero_radical_proved=False,
                  full_ambient_removal_proved=False,
                  global_rows_sha256=expected,
                  input_sha256={str(p.resolve().relative_to(ROOT)): digest(p) for p in sources},
                  output_sha256={'support.m2': digest(script), 'support.log': digest(log)},
                  seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
