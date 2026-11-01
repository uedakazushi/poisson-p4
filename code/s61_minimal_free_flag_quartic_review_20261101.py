#!/usr/bin/env python3
"""Independent constant lower-ideal Hilbert check and proof-scope audit of true quartics."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def run():
    started = time.time()
    stems = ['UNIVERSAL_QUADRATIC', 'UNIVERSAL_CUBIC', 'UNIVERSAL_QUARTIC']
    paths = [ROOT / f'reports/S61_MINIMAL_FREE_FLAG_{stem}_20261101.json' for stem in stems]
    q, c, r = [json.loads(path.read_text()) for path in paths]
    for path, expected in r['input_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    v = s.symbols('v0:11')
    quadratics = [sum(s.Rational(coefficient) * s.prod(v[j] for j in mon)
                      for coefficient, mon in zip(row, q['quadratic_monomials']))
                  for row in q['constant_quadratic_rref']]
    cubics = list(map(s.sympify, c['constant_cubic_polynomials']))

    def m2_expression(polynomial):
        terms = []
        for exponent, coefficient in s.Poly(polynomial, *v).terms():
            terms.append('(' + str(coefficient) + ')'
                         + ''.join('*v_' + str(j) + '^' + str(d)
                                   for j, d in enumerate(exponent) if d))
        return '+'.join(terms) or '0'

    program = 'R=QQ[v_0..v_10];I=ideal(' + ','.join(map(m2_expression, quadratics + cubics)) + ');\n'
    program += 'assert(hilbertFunction(4,R/I)==406);print("INDEPENDENT LOWER RANK595 PASS");\n'
    env = os.environ.copy()
    env['DYLD_LIBRARY_PATH'] = str(ROOT / 'work/m2-runtime/lib') + ':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    with tempfile.TemporaryDirectory(prefix='free-flag-quartic-review-') as directory:
        script = Path(directory) / 'verify.m2'
        script.write_text(program)
        result = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(script)],
                                env=env, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0 and 'INDEPENDENT LOWER RANK595 PASS' in result.stdout, result.stdout + result.stderr
    reduced = s.Matrix(r['new_quartic_rref']).applyfunc(s.sympify)
    assert reduced.rows == 15 and reduced[:, r['new_quartic_pivots']] == s.eye(15)
    parameters = s.symbols('a0:5')
    entries = [(i, j, value) for i, row in enumerate(reduced.tolist()) for j, value in enumerate(row)
               if value.free_symbols]
    assert len(entries) == 1
    i, j, value = entries[0]
    assert value.free_symbols == {parameters[1]} and s.degree(value, parameters[1]) == 1
    assert r['monomials'][r['quotient_free_columns'][j]] == [3, 3, 3, 10]
    coefficients = s.Matrix(r['new_quartic_coefficients']).applyfunc(s.sympify)
    witness = coefficients[r['new_quartic_selected_rows'], :].subs(dict(zip(parameters, range(1, 6)))).det()
    assert witness != 0
    sources = paths + [Path(__file__), ROOT / 'code/s61_minimal_free_flag_universal_quartic_20261101.py',
                       ROOT / 'research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.md',
                       ROOT / 'research/notes/S61_MINIMAL_FREE_FLAG_LOW_REVIEW_20261101.md']
    return dict(id='S61_MINIMAL_FREE_FLAG_QUARTIC_REVIEW_20261101', all_passed=True,
                independent_lower_quartic_rank=595, stored_extra_basis_rank=15,
                selected_coefficient_minor_at_sample=str(witness), single_parameter_entry_verified=True,
                accepted_true_hilbert_function=[1, 11, 54, 168, 391],
                scope='Independent lower-ideal Hilbert and stored-basis checks; conceptual review of filtered completeness.',
                limitations=['No independent raw fourth-jet contraction or main producer replay.',
                             'No fifth or higher initial equality, flatness, t-nilpotence, or endpoint removal.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_QUARTIC_REVIEW_20261101.md',
                m2_log=result.stdout + result.stderr,
                input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sources}, seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=report['seconds'])))
