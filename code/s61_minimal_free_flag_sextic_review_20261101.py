#!/usr/bin/env python3
"""Independent full polynomial-ideal HF check and sextic proof/source review."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def independent_hilbert(data):
    saved = {degree: json.loads((ROOT/('reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_'+degree+'_20261101.json')).read_text())
             for degree in ['QUADRATIC', 'CUBIC', 'QUARTIC', 'QUINTIC']}
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
        record = saved[label]
        monomials = [record[monomial_key][j] for j in record['quotient_free_columns']]
        polys.extend(polynomial(row, monomials) for row in record[rows])
    monomials = [data['sixth_monomials'][j] for j in data['quotient_free_columns']]
    sixth = [polynomial([str(coefficient) for column, coefficient in row],
                        [monomials[column] for column, coefficient in row])
             for row in data['constant_sixth_rref_sparse']]
    program = 'B=QQ[a];K=frac B;R=K[v_0..v_10];I=ideal('+','.join(polys)+');\n'
    program += 'assert(hilbertFunction(6,R/I)==1237);J=I+ideal('+','.join(sixth)+');\n'
    program += 'assert(apply(toList(0..6),d->hilbertFunction(d,R/J))=={1,11,54,168,391,740,1203});print "INDEPENDENT SEXTIC HF PASS";\n'
    with tempfile.TemporaryDirectory(prefix='free-flag-sixth-review-') as directory:
        source = Path(directory)/'review.m2'
        source.write_text(program)
        env = os.environ.copy()
        env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        process = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(source)],
                                 env=env, capture_output=True, text=True, check=True, timeout=60)
    assert 'INDEPENDENT SEXTIC HF PASS' in process.stdout
    return process.stdout+process.stderr


def verify():
    start = time.time()
    primary = 'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_SEXTIC_20261101.json'
    data = json.loads((ROOT/primary).read_text())
    assert data['all_passed']
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
               for path, digest in data['input_sha256'].items())
    assert data['relation_rows'] == 3263
    assert data['old_residual_rows']+data['independent_product_syzygies'] == 1180+2083 == 3263
    assert data['new_sixth_rank'] == 34 and data['nonzero_canonical_coefficients'] == 179
    assert data['full_fifth_leading_compatibility_checked'] is True
    assert data['selected_minor_determinant'] % data['selected_minor_prime'] != 0
    log = independent_hilbert(data)
    paths = list(data['input_sha256'])+[primary,
        'research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_SEXTIC_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_MINIMAL_FREE_FLAG_SEXTIC_REVIEW_20261101.md']
    return dict(id='S61_MINIMAL_FREE_FLAG_SEXTIC_REVIEW_20261101', all_passed=True,
                independent_full_eleven_variable_ideal_hilbert_check=True,
                independent_field='QQ(a1)', known_sixth_quotient_before_new_relations=1237,
                accepted_hilbert_function=[1, 11, 54, 168, 391, 740, 1203],
                proof_and_source_review_accepted=True, primary_input_hashes_verified=True,
                full_fifth_leading_compatibility_assert_source_audited=True,
                independent_raw6_or_3263_tail_reimplementation=False, m2_log=log,
                limitations=['The independent computation is on the reported homogeneous generators.',
                             'Raw135 contraction, full filtered tails, fifth compatibility and universal membership use the source-audited primary producer.',
                             'The result is over the five-parameter function field, not every core value.',
                             'No full associated-graded equality, ambient nilpotence or saturation removal.'],
                input_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('FREE FLAG SEXTIC REVIEW PASS', result['seconds'])
