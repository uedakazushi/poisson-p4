#!/usr/bin/env python3
"""Independent rooted-flow/weighted-determinant audit of actual42 Euler delta28."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63_REVIEW'
P = 32003


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def residue(value):
    a, sep, b = str(value).partition('/')
    den = int(b) % P if sep else 1
    assert den
    return int(a) % P * pow(den, -1, P) % P


def accepted(path):
    d = load(path)
    assert d['all_passed'] is True
    for name, digest in d.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in d.get('evidence_sha256', {}).items():
        assert sha(path.parent / name) == digest, name
    return d


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    output.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63'
    accepted(producer / 'report.json')
    saved = load(producer / 'certificate.json.gz')
    dependency_paths = [
        ROOT / 'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW.json',
        ROOT / 'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json',
        ROOT / 'reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63/report.json',
        ROOT / 'reports/S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63_REVIEW.json']
    for path in dependency_paths:
        accepted(path)
    modelpath = ROOT / 'reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63/coefficient-model-mod32003.json.gz'
    model = load(modelpath)
    assert model['prime'] == P
    presentationpath = ROOT / 'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz'
    presentation = load(presentationpath)
    specpath = ROOT / 'code/data/s61_affine_normal_global_quintic_20261103_v63.json'
    spec = load(specpath)

    # Independently stream all exact QQ reduction coefficients.  The accepted
    # integer receipt supplies the two inverse identities, not a modular lift.
    finitepath = ROOT / 'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/certificate.json'
    finite = load(finitepath)
    native = ROOT / 'work/normal-finite5-public-v63/cross-output.txt'
    assert sha(native) == finite['native_output_sha256']
    assert native.stat().st_size == finite['native_output_bytes']
    reduced_P = []
    with native.open() as f:
        assert f.readline().split() == ['126', '945']
        assert list(map(int, f.readline().split())) == finite['selected_cross_rows']
        for _ in range(126):
            row = f.readline().split()
            assert len(row) == 126
            for value in row:
                residue(value)  # all selected QQ inverse denominators are p-units
        for _ in range(126):
            row = f.readline().split()
            assert len(row) == 945
            reduced_P.append([residue(value) for value in row])
        assert not f.read().strip()
    assert reduced_P == model['P_coefficients']

    # Check the full actual normal weight bound for every remainder coefficient.
    weights = [1] * 42 + [4, 4, 5]
    coefficient_bounds = [0] * 5
    for remainder in presentation['all134_remainders']:
        for degree, terms in remainder:
            for mon, coefficient in terms:
                residue(coefficient)
                coefficient_bounds[degree] = max(coefficient_bounds[degree], sum(weights[j] for j in mon))
                if degree == 4:
                    assert len(mon) <= 1 and 44 not in mon
    assert coefficient_bounds[0] <= 8 and max(coefficient_bounds[1:4]) <= 5 and coefficient_bounds[4] <= 4

    # Universal weighted companion determinant, using all 24 permutation terms.
    H, Q = s.symbols('H Q')
    L = s.symbols('L0:5')
    g = s.symbols('g0:5')
    T = s.zeros(5)
    for i in range(4):
        T[i + 1, i] = 1
    T[:, 4] = s.Matrix([H - g[0], -g[1], -g[2], -g[3], -g[4]])
    vectors = [s.Matrix([L[0] + Q, *L[1:]])]
    for _ in range(3):
        vectors.append(T * vectors[-1])
    bottom = s.Matrix.hstack(*vectors)[1:, :]
    determinant = 0
    for perm in itertools.permutations(range(4)):
        inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4))
        determinant += (-1) ** inversions * s.prod(bottom[i, perm[i]] for i in range(4))
    variables = (H, Q, *L, *g)
    polynomial = s.Poly(s.expand(determinant), *variables, domain=s.ZZ)
    variable_weights = (5, 8, 5, 5, 5, 5, 4, 0, 0, 0, 0, 0)
    top = 0
    max_weight = 0
    for exponents, c in polynomial.terms():
        weight = sum(a * b for a, b in zip(exponents, variable_weights))
        max_weight = max(max_weight, weight)
        if weight == 28:
            top += c * s.prod(v ** e for v, e in zip(variables, exponents))
    expected_top = H ** 2 * L[3] ** 2 * L[4] ** 2 - H * L[3] ** 3 * Q - L[4] * Q ** 3
    assert max_weight == 28 and s.expand(top - expected_top) == 0
    assert s.expand(s.sympify(saved['exact_universal_degree28_top'], locals={str(v): v for v in variables}) - expected_top) == 0

    # Bind all plane cofactor identities and check their exact QQ denominators.
    z = s.symbols('z0:3')
    planes = [s.Poly(s.sympify(p), *z, domain=s.QQ) for p in spec['plane']]
    relations = sorted(spec['relations'], key=lambda r: r['target_index'])
    highest = []
    for index, relation in enumerate(relations):
        cofactors = [s.Poly(s.sympify(p), *z, domain=s.QQ) for p in relation['cofactors']]
        assert s.Poly(sum((a.as_expr() * b.as_expr() for a, b in zip(cofactors, planes)), s.Integer(0)) - s.sympify(relation['target']), *z).is_zero
        degree = 2 if index < 2 else 3
        high = []
        for p in cofactors:
            row = []
            for ex, coefficient in p.terms():
                reduced = residue(coefficient)
                if sum(ex) == degree and coefficient:
                    row.append((tuple(j for j, count in enumerate(ex) for _ in range(count)), reduced))
            high.append(row)
        highest.append(high)
    for p in planes:
        for coefficient in p.coeffs():
            residue(coefficient)

    # Rooted composition, independent of the author's Hessian/Jacobian formula.
    # For nilpotent u,v,w (each square zero), the mixed coefficient of
    # f(Flow_c(w, Flow_b(v, Flow_a(u,x)))) is D_a D_b D_c f(x).
    words = model['normal_monomials']
    qcols = [i for i, word in enumerate(words) if len(word) == 2]
    qa = np.array([words[i][0] for i in qcols], dtype=np.int64)
    qb = np.array([words[i][1] for i in qcols], dtype=np.int64)
    vectors_P = np.array(reduced_P, dtype=np.int64)[:, qcols].reshape(3, 42, len(qcols))
    forms_K = np.array(model['K_coefficients'], dtype=np.int64)[:, qcols]
    subsets = [[u for u in range(8) if u & mask == u] for mask in range(8)]

    def quadratic_jet(coefficients, point_jet):
        monomial = np.zeros((8, len(qcols)), dtype=np.int64)
        for mask in range(8):
            for sub in subsets[mask]:
                monomial[mask] += point_jet[sub, qa] * point_jet[mask ^ sub, qb]
            monomial[mask] %= P
        return (coefficients @ monomial.T % P).T

    def rooted(word):
        point = np.zeros((8, 42), dtype=np.int64)
        point[0, 0] = 1
        for j, vector_index in enumerate(word):
            field = quadratic_jet(vectors_P[vector_index], point)
            bit = 1 << j
            updated = point.copy()
            for mask in range(8):
                if not mask & bit:
                    updated[mask | bit] = (updated[mask | bit] + field[mask]) % P
            point = updated
        return quadratic_jet(forms_K, point)[(1 << len(word)) - 1]

    second = {word: rooted(word) for word in itertools.combinations_with_replacement(range(3), 2)}
    third = {word: rooted(word) for word in itertools.combinations_with_replacement(range(3), 3)}
    assert {','.join(map(str, k)): v.tolist() for k, v in second.items()} == saved['all_ordered_second_derivative_values']
    assert {','.join(map(str, k)): v.tolist() for k, v in third.items()} == saved['all_ordered_third_derivative_values']
    S4 = [-sum(coefficient * int(second[word][i]) for i, row in enumerate(cofs) for word, coefficient in row) * pow(6, -1, P) % P for cofs in highest[:2]]
    H5 = -sum(coefficient * int(third[word][i]) for i, row in enumerate(highest[2]) for word, coefficient in row) * pow(24, -1, P) % P
    assert S4 == [15338, 11727] == saved['S1_S2_degree4']
    assert H5 == 31338 == saved['H_degree5']

    # Independent polynomial substitution in r0's actual weighted highest parts.
    replacements = {42: S4[0], 43: S4[1], 44: H5}
    values = []
    first = dict(presentation['all134_remainders'][0])
    for hdegree, weight in ((0, 8), (3, 5), (4, 4)):
        value = 0
        for mon, coefficient in first[hdegree]:
            if sum(weights[j] for j in mon) != weight:
                continue
            term = residue(coefficient)
            for j in mon:
                term = term * replacements.get(j, int(j == 0)) % P
            value = (value + term) % P
        values.append(value)
    Q8, L3, L4 = values
    assert values == [30427, 13746, 2943]
    nonzero = (H5 ** 2 * L3 ** 2 * L4 ** 2 - H5 * L3 ** 3 * Q8 - L4 * Q8 ** 3) % P
    assert nonzero == 31492 == saved['actual_delta_degree28_witness']
    assert 3 + 2583 + 210 + 42 == 2838
    assert saved['scalar_generator_count'] == 2838 and saved['scalar_generator_degree_bound'] == 4
    evidence = dict(
        all_exact_cross_inverse_denominators_p_units=126 * 126,
        all_exact_cross_P_coefficients_p_units_and_reductions_bound=126 * 945,
        all_plane_cofactor_QQ_identities_and_p_unit_denominators=True,
        all134_remainder_actual_coefficient_weight_bounds=coefficient_bounds,
        universal_determinant_24_permutation_terms_checked=True,
        weighted_degree=28,
        weighted_top=str(s.expand(top)),
        independent_square_zero_flow_second={','.join(map(str, k)): v.tolist() for k, v in second.items()},
        independent_square_zero_flow_third={','.join(map(str, k)): v.tolist() for k, v in third.items()},
        S1_S2_degree4=S4, H_degree5=H5, Q_degree8=Q8, L3_degree5=L3, L4_degree4=L4,
        nonzero_mod32003=nonzero,
        actual_QQ_nonzero_via_p_integral_positive_witness=True,
        Euler_actual_algebra_and_scalar_criterion_proof_accepted=True,
        actual_scalar_membership_proved=False,
        normal_degree_cofactor_necessary_lower_bound=24,
        cofactor_upper_bound_or_nonmembership_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
             ROOT / 'code/s61_affine_actual42_euler_delta28_20261103_v63.py',
             ROOT / 'research/notes/S61_AFFINE_ACTUAL42_EULER_DELTA28_20261103_v63.md',
             producer / 'report.json', producer / 'certificate.json.gz',
             modelpath, specpath, presentationpath, finitepath, *dependency_paths]
    report = dict(id=STEM, all_passed=True, producer_rerun=False,
                  actual_degree28_and_cofactor_lower24_accepted=True,
                  actual_scalar_criterion_accepted=True, actual_membership_or_exhaustion_proved=False,
                  independent_method='24 permutation determinant; square-zero flow rooted derivatives; exact QQ cofactor identities and streaming p-unit cross-coefficient reduction.',
                  accepted_native_QQ_source_sha256=sha(native),
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
                  evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
                  seconds=time.monotonic() - start)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)
