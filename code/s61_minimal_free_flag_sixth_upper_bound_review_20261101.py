#!/usr/bin/env python3
"""Independent QQ(a1) lower-ideal Hilbert check and sixth-bound proof review."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s
from s61_minimal_free_flag_sixth_upper_bound_20261101 import run as primary_run

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.time()
    files = {degree: 'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_'+degree+'_20261101.json'
             for degree in ['QUADRATIC', 'CUBIC', 'QUARTIC', 'QUINTIC']}
    saved = {degree: json.loads((ROOT/path).read_text()) for degree, path in files.items()}
    def polynomial(row, monomials):
        return '+'.join('('+c.replace('**', '^').replace('a1', 'a')+')*'
                        +'*'.join('v_'+str(i) for i in mon)
                        for c, mon in zip(row, monomials) if c != '0') or '0'
    polys = [polynomial(row, saved['QUADRATIC']['quadratic_monomials'])
             for row in saved['QUADRATIC']['constant_quadratic_rref']]
    v = s.symbols('v0:11')
    for expression in saved['CUBIC']['constant_cubic_polynomials']:
        poly = s.Poly(s.sympify(expression), *v)
        polys.append('+'.join('('+str(c)+')*'+'*'.join('v_'+str(i) for i, e in enumerate(mon) for _ in range(e))
                              for mon, c in poly.terms()))
    for label, rows, monomial_key in [('QUARTIC', 'new_quartic_rref', 'monomials'),
                                       ('QUINTIC', 'constant_fifth_rref', 'fifth_monomials')]:
        data = saved[label]
        monomials = [data[monomial_key][j] for j in data['quotient_free_columns']]
        polys.extend(polynomial(row, monomials) for row in data[rows])
    program = 'B=QQ[a];K=frac B;R=K[v_0..v_10];I=ideal('+','.join(polys)+');\n'
    program += 'assert(apply(toList(0..6),d->hilbertFunction(d,R/I))=={1,11,54,168,391,740,1237});print "INDEPENDENT KNOWN SIXTH 1237 PASS";\n'
    with tempfile.TemporaryDirectory(prefix='free-flag-sixth-review-') as directory:
        source = Path(directory)/'known.m2'
        source.write_text(program)
        env = os.environ.copy()
        env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        process = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(source)],
                                 env=env, capture_output=True, text=True, check=True, timeout=45)
        assert 'INDEPENDENT KNOWN SIXTH 1237 PASS' in process.stdout
    primary_path = 'reports/S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_20261101.json'
    original = json.loads((ROOT/primary_path).read_text())
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
               for path, digest in original['input_sha256'].items())
    replay = primary_run()
    assert replay['all_passed'] and replay['exact_modular_ideal_rank'] == 9808
    assert replay['modular_input_sha256'] == original['modular_input_sha256']
    assert replay['filtered_trace_sha256'] == original['filtered_trace_sha256']
    assert 1237-replay['generic_sixth_hilbert_upper_bound'] == 34
    paths = list(files.values())+[primary_path, str(Path(__file__).relative_to(ROOT)),
                                  'code/s61_minimal_free_flag_sixth_upper_bound_20261101.py',
                                  'research/notes/S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_20261101.md']
    return dict(id='S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_REVIEW_20261101', all_passed=True,
                independent_field='QQ(a1)', independent_known_ideal_hilbert_function=[1, 11, 54, 168, 391, 740, 1237],
                accepted_generic_HF6_upper_bound=1203, accepted_new_sixth_rank_lower_bound=34,
                primary_producer_reproduced=True, modular_input_and_trace_hashes_match=True,
                independent_raw6_reimplementation=False, m2_log=process.stdout+process.stderr,
                limitations=['The independent calculation concerns the known initial ideal, not raw sixth coefficients.',
                             'The raw sixth jet and both modular closures are producer reproduction.',
                             'No generic equality, full initial ideal, nilpotence, or stratum removal.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_REVIEW_20261101.md',
                input_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as output:
        json.dump(result, output, indent=2)
        output.write('\n')
    print('FREE FLAG SIXTH UPPER REVIEW PASS', result['seconds'])
