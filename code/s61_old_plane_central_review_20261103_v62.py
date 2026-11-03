#!/usr/bin/env python3
"""Check the integral plane quotient and degreewise lifting arithmetic.

The saved M2 finite-field full Groebner calculation is input evidence; this
review does not reimplement that engine or claim ambient t-nilpotence.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003
WEIGHTS = [1]*7+[2]*4+[3]*2+[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(rows, output):
    model = ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json'
    assert sha(rows) == json.loads(model.read_text())['row_output_sha256']
    rational_terms = 0
    plane_terms = 0
    with rows.open() as stream:
        assert list(map(int, stream.readline().split())) == [14, 136]
        for row in range(136):
            weight, count = map(int, stream.readline().split())
            restriction = defaultdict(Fraction)
            for _ in range(count):
                scalar, *powers = stream.readline().split()
                scalar = Fraction(scalar)
                hp, *exponents = map(int, powers)
                assert len(exponents) == 14
                assert scalar.denominator % PRIME
                assert sum(a*b for a, b in zip(exponents, WEIGHTS)) == weight
                rational_terms += 1
                # rho=0, all v0..v8=v11=0, v9=h*u, v10=u, v12=v.
                if row == 135 or any(exponents[i] for i in [*range(9), 11, 13]):
                    continue
                hp += exponents[9]
                restriction[(exponents[9]+exponents[10], exponents[12], hp % 2)] += scalar*(-2)**(hp//2)
                plane_terms += 1
            if row < 135:
                assert not any(restriction.values())
        assert not stream.read().strip()
    assert (3083**2+2) % PRIME == 0 and (2*3083) % PRIME
    T = s.Symbol('T')
    numerator = 1+7*T+22*T**2+27*T**3-26*T**4-56*T**5-30*T**6+19*T**7+34*T**8+3*T**9
    kernel = s.cancel((numerator-1)/((1-T**2)*(1-T**3)))
    assert s.expand(kernel) == 7*T+22*T**2+34*T**3+3*T**4
    central_dir = ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_CENTRAL_20261103_v62'
    flat_dir = ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_QQ_FLAT_20261103_v62'
    for folder in [central_dir, flat_dir]:
        report = json.loads((folder/'report.json').read_text())
        assert report['execution_completed'] and not report['full_ambient_removal_proved']
        assert report['global_rows_sha256'] == sha(rows)
        for p, digest in report['input_sha256'].items():
            assert sha(ROOT/p) == digest
    log = (central_dir/'central.log').read_text()
    assert 'GB_GENERATORS=171' in log and 'RADICAL_PROVED=true' in log
    assert 'T_POWER=(2,3,4,5,6,7,8,9,10,11,12)' in log
    # Independently recompute the Tor-derived bounds; these are upper bounds
    # because Tor_1(A,k) may map nontrivially into Tor_1(C,k).
    section_results = []
    for name in ['S61_AFFINE_ENLARGED_SECTION_20261103_v62', 'S61_AFFINE_ENLARGED_SECTION_P28_20261103_v62']:
        folder = ROOT/'reports'/name
        report = json.loads((folder/'report.json').read_text())
        for p, digest in report['input_sha256'].items():
            assert sha(ROOT/p) == digest
        for p, digest in report['evidence_sha256'].items():
            assert sha(folder/p) == digest
        r = report['retained_normal_count']
        central = report['hilbert_function_mod_32003']
        c = [1, 4, 5, 5, 5]
        bounds = [r*c[d-1]+central[d]-c[d] for d in range(1, 5)]
        assert bounds == report['finite_normal_module_generator_bounds_by_degree']
        assert sum(bounds) == report['finite_normal_module_generator_bound']
        section_results.append(dict(section=report['projective_section_dimension'], base_variables=r,
                                    generator_degree_bounds=bounds, total_bound=sum(bounds)))
    paths = [Path(__file__), model,
             ROOT/'research/notes/S61_MARKED_SIGMA_OLD_PLANE_CENTRAL_20261103_v62.md',
             ROOT/'research/notes/S61_AFFINE_ENLARGED_SECTION_20261103_v62.md']
    for folder in [central_dir, flat_dir]:
        paths.extend(folder/f for f in ['report.json', 'central.m2', 'central.log'])
    for name in ['S61_AFFINE_ENLARGED_SECTION_20261103_v62', 'S61_AFFINE_ENLARGED_SECTION_P28_20261103_v62']:
        paths.extend(ROOT/'reports'/name/f for f in ['report.json', 'review.json'])
    result = dict(id='S61_OLD_PLANE_CENTRAL_AND_ENLARGED_SECTION_20261103_v62_REVIEW', all_passed=True,
                  original_rows_sha256=sha(rows), rational_coefficients_checked=rational_terms,
                  all_rational_denominators_p_units=True, original_weight_homogeneity_checked=True,
                  all135_integral_plane_restrictions_zero=True, surviving_plane_terms=plane_terms,
                  split_simple_root_and_Hensel_embedding=True,
                  special_kernel_hilbert=[0,7,22,34,3], characteristic_zero_kernel_length_upper=66,
                  degreewise_Nakayama_proof_accepted=True, ambient_nilpotence_asserted=False,
                  source_finite_field_Groebner_engine_independently_reimplemented=False,
                  enlarged_section_Tor_and_generation_proof_accepted=True, enlarged_sections=section_results,
                  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')
    print('PASS', rational_terms, 'rational coefficients;', plane_terms, 'plane terms;', section_results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.rows, args.output)
