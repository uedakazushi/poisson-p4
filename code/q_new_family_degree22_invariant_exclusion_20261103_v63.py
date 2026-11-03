#!/usr/bin/env python3
"""One conditional exact QQ21minor after the rank21 degree22 product probe."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import resource
import shlex
import subprocess
import sys
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

from q_new_family_degree16_invariant_probe_20261103_v63 import ROOT, PRIME, read, sha, matrix, mod, necklaces
from weight1_structure import sp_basis

STEM = 'Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63'


def cap():
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (2_200_000_000, 2_200_000_000))


def main(out, probe):
    started = time.monotonic()
    sys.set_int_max_str_digits(0)
    out.mkdir(parents=True, exist_ok=False)
    pp = probe / 'report.json'
    data = read(pp)
    assert data['all_passed'] and data['actual_restriction_rank'] == data['full_character_dimension22'] == 21
    assert data['existing_point_count'] == 36 and data['new_beta_points_generated'] == 0
    for path, digest in data['input_sha256'].items():
        assert sha(ROOT / path) == digest, path
    invpath = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    old18 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    old20 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/report.json'
    old14 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/report.json'
    invreview = ROOT / 'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
    assert read(old20)['all_passed'] and read(old14)['all_passed']
    assert read(invreview)['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
    info = read(invpath)
    O = matrix(info['symplectic_form16'])
    reps = [matrix(A) for A in info['all_ten_representation_matrices']]
    gs = sp_basis()
    Gi = s.Matrix([[s.trace(A * B) for B in gs] for A in gs]).inv()
    selected_rows, selected_columns = data['selected_rows'], data['selected_columns']
    assert len(selected_rows) == len(selected_columns) == 21
    source = [data['same_accepted_QQ_samples'][i] for i in selected_rows]
    for row, sample in zip(selected_rows, source):
        assert [mod(x) for x in sample['beta16']] == data['actual_beta16_modp'][row]
        assert Fraction(sample['value_pivot']) and Fraction(sample['value_pfaffian'])
    native = ROOT / 'code/q_new_family_degree22_invariant_exact_20261103_v63.cpp'
    np = out / 'native-input.txt'
    with np.open('w') as f:
        def write(M):
            f.write(f'{M.rows} {M.cols}\n')
            for row in M.tolist():
                f.write(' '.join(map(str, row)) + '\n')
        f.write('21 24\n')
        for M in [O, O.inv(), Gi, *reps]:
            write(M)
        for word in necklaces(4)[:5]:
            f.write(''.join(map(str, word)) + '\n')
        for word in read(old18)['selected_words']:
            f.write(''.join(map(str, word)) + '\n')
        for pair in [[1, 2], [1, 3], [2, 3]]:
            f.write(' '.join(map(str, pair)) + '\n')
        f.write(' '.join(map(str, selected_columns)) + '\n')
        for sample in source:
            write(s.Matrix([s.Rational(x) for x in sample['beta16']]))
    binary = ROOT / 'work/new-family-degree22-exact-native-v63'
    flags = shlex.split(subprocess.check_output(['pkg-config', '--cflags', '--libs', 'flint'], text=True))
    subprocess.run(['c++', '-std=c++17', '-O2', str(native), *flags, '-lgmp', '-o', str(binary)], check=True)
    evidence = out / 'exact-evaluation.json'
    with (out / 'native.log').open('w') as f:
        subprocess.run([str(binary), str(np), str(evidence)], stdout=f, stderr=subprocess.STDOUT,
                       check=True, timeout=65, preexec_fn=cap)
    exact = read(evidence)
    assert exact['all_passed'] and len(exact['evaluation']) == 21
    for i, row in enumerate(selected_rows):
        assert [mod(x) for x in exact['evaluation'][i]] == data['actual_evaluation'][row]
        expected = data['factor_values'][row]
        assert [mod(x) for x in exact['factors_F4_Inv8five_Inv14three_Inv18nine'][i]] == \
            [expected['F4'], *expected['Inv8'], *expected['Inv14'], *expected['Inv18']]
    det = Fraction(exact['exact_selected21_minor'])
    assert det
    modular = s.Matrix([[data['actual_evaluation'][r][c] for c in selected_columns] for r in selected_rows])
    md = int(DomainMatrix.from_Matrix(modular).convert_to(s.GF(PRIME)).det()) % PRIME
    assert md == mod(str(det)) != 0
    samples = dict(probe_rows=selected_rows, selected_columns=selected_columns,
                   selected_definitions=[data['definitions'][j] for j in selected_columns],
                   existing_accepted_QQ_source_points=source)
    (out / 'samples.json').write_text(json.dumps(samples, indent=2) + '\n')
    files = [Path(__file__), native, pp, invpath, old18, old20, old14, invreview,
             ROOT / 'code/weight1_structure.py', ROOT / 'code/q_new_family_degree16_invariant_probe_20261103_v63.py']
    result = dict(id=STEM, all_passed=True, full_degree22_invariant_dimension=21,
                  universal_product_candidate_count=24, exact_QQ_restriction_rank=21,
                  exact_selected21_minor_nonzero=True, exact21_minor_file='exact-evaluation.json',
                  exact21_minor_numerator_bits=abs(det.numerator).bit_length(),
                  exact21_minor_denominator_bits=det.denominator.bit_length(), exact21_minor_mod32003=md,
                  all504_QQ_product_values_verified_modp=True, all378_QQ_factor_values_verified_modp=True,
                  native_degree14_gradients_use_direct_dtheta_dH=True,
                  existing_accepted_QQ_records_reused=True, no_new_source_or_beta_points=True,
                  complete_degree22_space_spanned_and_restriction_injective=True,
                  defining_homogeneous_invariant_degree_lower_bound=24,
                  exact_native_run_count=1, no_degree24_probe_or_equation_claimed=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in files},
                  evidence_sha256={p.name: sha(p) for p in [np, evidence, out / 'native.log', out / 'samples.json']},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'exact21_minor_numerator_bits',
          'exact21_minor_denominator_bits', 'exact21_minor_mod32003', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True)
    a = p.parse_args()
    main(a.output.resolve(), a.probe.resolve())
