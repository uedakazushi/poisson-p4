#!/usr/bin/env python3
"""Exact finite-field J10/base quotient and all-r parameter-adic cutoff."""
import argparse
import base64
import gzip
import hashlib
import json
import os
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement as cwr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def verify(binary):
    started = time.time()
    source = 'reports/S61_MINIMAL_FREE_FLAG_INITIAL_PROBE_20261101.json'
    record = json.loads((ROOT / source).read_text())
    raw = gzip.decompress(base64.b64decode(record['initial_rows_compressed']))
    assert hashlib.sha256(raw).hexdigest() == record['initial_rows_sha256']
    assert record['prime'] == 32003 and record['maximum_degree'] == 10
    monomials = [m for d in range(11) for m in cwr(range(11), d)]
    equations = []
    terms = 0
    for line in raw.decode().splitlines()[1:]:
        row = list(map(int, line.split()))
        assert len(row) == 1 + 2 * row[0]
        parts = []
        for index, coefficient in zip(row[1::2], row[2::2]):
            monomial = monomials[index]
            if any(j >= 8 for j in monomial):
                continue
            parts.append(str(coefficient) + '*' + '*'.join('x_' + str(j) for j in monomial))
        terms += len(parts)
        if parts:
            equations.append('+'.join(parts))
    assert len(equations) == 434 and terms == 7322
    expected = [1, 8, 24, 38, 38, 23, 8, 3, 1, 0, 0, 0]
    script = 'R=ZZ/32003[x_0..x_7];J=ideal(' + ','.join(equations) + ');\n'
    script += 'M=R/J;h=apply(toList(0..11),d->hilbertFunction(d,M));\n'
    script += 'assert(h=={' + ','.join(map(str, expected)) + '});\n'
    script += 'assert(dim M==0 and degree M==144);print("CUTOFF",dim M,degree M,h);\n'
    script += 'print("ALL_CUTOFF_CHECKS_PASS");stdio<<flush;\n'
    env = os.environ.copy()
    local_libraries = ROOT / 'work/m2-runtime/lib'
    if local_libraries.is_dir():
        env['DYLD_LIBRARY_PATH'] = str(local_libraries) + ':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    with tempfile.TemporaryDirectory(prefix='poisson-free-M0-cutoff-') as temporary:
        path = Path(temporary) / 'cutoff.m2'
        path.write_text(script)
        result = subprocess.run([str(binary), '--script', str(path)], env=env,
                                capture_output=True, text=True, timeout=60, check=True)
    assert 'ALL_CUTOFF_CHECKS_PASS' in result.stdout
    paths = [source, str(Path(__file__).relative_to(ROOT)),
             'research/notes/S61_MINIMAL_FREE_FLAG_R_ADIC_CUTOFF_20261102.md',
             'research/notes/S61_MINIMAL_FREE_FLAG_INITIAL_PROBE_20261101.md',
             'reports/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.json']
    return dict(id='S61_MINIMAL_FREE_FLAG_R_ADIC_CUTOFF_20261102', all_passed=True,
                prime=32003, base_variables=[8, 9, 10], specialized_generators=434,
                specialized_terms=7322, Hilbert_values=expected, quotient_length=144,
                highest_nonzero_degree=8, polynomial_module_generator_degree_upper=8,
                all_r_cutoff='m^(r+8) is contained in (v8,v9,v10)^r S for every r>=1',
                actual_transverse_length_from_previous_record=141,
                M2_script=script, M2_stdout=result.stdout,
                limitations=['One specified finite-field point only.',
                             'J10/base has length144, not the actual transverse length141.',
                             'No equality of J10 with the full initial ideal, no t-nilpotence or relative removal.'],
                input_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m2-binary', type=Path,
                        default=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify(args.m2_binary)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('FREE M0 R-ADIC CUTOFF PASS', report['seconds'])
