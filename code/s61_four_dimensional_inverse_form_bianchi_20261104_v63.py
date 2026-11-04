#!/usr/bin/env python3
"""Exact universal 4D inverse-form/Bianchi comparison; no S61 expansion."""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_FOUR_DIMENSIONAL_INVERSE_FORM_BIANCHI_20261104_v63'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(60)
    output.mkdir(parents=True, exist_ok=False)
    pairs = list(it.combinations(range(4), 2))
    jet2 = list(it.combinations_with_replacement(range(4), 2))
    names = ([f'p{i}{j}' for i, j in pairs]
             + [f'p{i}{j}_{k}' for i, j in pairs for k in range(4)]
             + [f'p{i}{j}_{k}{l}' for i, j in pairs for k, l in jet2]
             + [f'h{i}' for i in range(4)]
             + [f'h{j}_{i}' for i in range(4) for j in range(4)])
    R, *gens = ring(','.join(names), s.QQ)
    z = dict(zip(names, gens))
    zero, one = R.zero, R.one
    P = [[zero for _ in range(4)] for _ in range(4)]
    p1, p2 = {}, {}
    for i, j in pairs:
        P[i][j], P[j][i] = z[f'p{i}{j}'], -z[f'p{i}{j}']
        for k in range(4):
            p1[i, j, k] = z[f'p{i}{j}_{k}']
        for k, l in jet2:
            p2[i, j, k, l] = z[f'p{i}{j}_{k}{l}']

    def D(poly, k):
        answer = zero
        for i, j in pairs:
            answer += poly.diff(P[i][j]) * p1[i, j, k]
            for l in range(4):
                answer += poly.diff(p1[i, j, l]) * p2[i, j, *sorted((k, l))]
        return answer

    def sign(indices):
        return (-1) ** sum(indices[i] > indices[j]
                          for i in range(len(indices)) for j in range(i + 1, len(indices)))

    def exterior_derivative(form):
        answer = {}
        for indices, poly in form.items():
            for k in range(4):
                if k not in indices:
                    key = tuple(sorted((k, *indices)))
                    answer[key] = answer.get(key, zero) + sign((k, *indices)) * D(poly, k)
        return {key: poly for key, poly in answer.items() if poly}

    def wedge(left, right):
        answer = {}
        for I, u in left.items():
            for J, v in right.items():
                if set(I).isdisjoint(J):
                    key = tuple(sorted((*I, *J)))
                    answer[key] = answer.get(key, zero) + sign((*I, *J)) * u * v
        return {key: poly for key, poly in answer.items() if poly}

    def multiply(form, poly):
        return {key: poly * value for key, value in form.items() if poly * value}

    def subtract(left, right):
        answer = dict(left)
        for key, poly in right.items():
            answer[key] = answer.get(key, zero) - poly
        return {key: poly for key, poly in answer.items() if poly}

    a, b, c, d, e, f = [P[i][j] for i, j in pairs]
    pf = a * f - b * e + c * d
    A = [[zero, -f, e, -d], [f, zero, -c, b],
         [-e, c, zero, -a], [d, -b, a, zero]]
    checks = {}
    checks['P_A_equals_A_P_equals_pf_identity'] = all(
        sum((P[i][k] * A[k][j] for k in range(4)), zero) == pf * int(i == j)
        and sum((A[i][k] * P[k][j] for k in range(4)), zero) == pf * int(i == j)
        for i in range(4) for j in range(4))
    dp = [D(pf, i) for i in range(4)]
    X = [sum((D(P[i][j], i) for i in range(4)), zero) for j in range(4)]
    J = {}
    for i, j, k in it.combinations(range(4), 3):
        J[i, j, k] = sum((P[i][l] * D(P[j][k], l)
                            + P[j][l] * D(P[k][i], l)
                            + P[k][l] * D(P[i][j], l) for l in range(4)), zero)
    dual = [(-1) ** i * J[tuple(j for j in range(4) if j != i)] for i in range(4)]
    checks['actual_J_dual_equals_minus_dp_minus_A_X'] = all(
        dual[i] == -dp[i] - sum((A[i][j] * X[j] for j in range(4)), zero)
        for i in range(4))
    formA = {(i, j): A[i][j] for i, j in pairs}
    formdp = {(i,): dp[i] for i in range(4)}
    dA = exterior_derivative(formA)
    E = subtract(multiply(dA, pf), wedge(formdp, formA))
    edual = [(-1) ** i * E[tuple(j for j in range(4) if j != i)] for i in range(4)]
    checks['E_dual_equals_minus_P_J_dual'] = all(
        edual[i] == -sum((P[i][j] * dual[j] for j in range(4)), zero)
        for i in range(4))
    checks['dE_equals_2_dp_wedge_dA'] = not subtract(
        exterior_derivative(E), multiply(wedge(formdp, dA), R(2)))
    checks['pf_dE_minus_2_dp_wedge_E_zero_in_all_second_jets'] = not subtract(
        multiply(exterior_derivative(E), pf), multiply(wedge(formdp, E), R(2)))
    checks['actual_Bianchi_zero_in_all_second_jets'] = not (
        sum((X[j] * dual[j] for j in range(4)), zero)
        - sum((P[i][j] * D(dual[j], i) for i in range(4) for j in range(4)), zero))

    h = [z[f'h{i}'] for i in range(4)]
    dh = [[z[f'h{j}_{i}'] for j in range(4)] for i in range(4)]
    B = sum((X[j] * h[j] for j in range(4)), zero) - sum(
        (P[i][j] * dh[i][j] for i in range(4) for j in range(4)), zero)
    # Independently build the same Schouten coefficient by alternating indices,
    # matching the exact convention of the frozen full135 operator.
    schouten = zero
    for i in range(4):
        for l in range(4):
            schouten += P[i][l] * dh[l][i]
    for i, k in pairs:
        other = tuple(j for j in range(4) if j not in (i, k))
        for l in range(4):
            if l not in other:
                tri = tuple(sorted((l, *other)))
                omitted = next(j for j in range(4) if j not in tri)
                schouten += ((-1) ** (i + k + omitted) * sign((l, *other))
                              * D(P[i][k], l) * h[omitted])
    checks['Schouten_convention_equals_X_dot_h_minus_P_dh'] = schouten == B
    div_minus_Ph = -sum((X[j] * h[j] for j in range(4)), zero) - sum(
        (P[i][j] * dh[i][j] for i in range(4) for j in range(4)), zero)
    deRham = pf * div_minus_Ph + 2 * sum(
        (dp[i] * P[i][j] * h[j] for i in range(4) for j in range(4)), zero)
    koszul = -2 * sum((dual[i] * P[i][j] * h[j]
                      for i in range(4) for j in range(4)), zero)
    checks['operator_identity_deRham_equals_pf_Bianchi_plus_Koszul'] = deRham == pf * B + koszul
    checks['Koszul_after_h_equals_actual_J_is_skew_zero'] = not sum(
        (dual[i] * P[i][j] * dual[j] for i in range(4) for j in range(4)), zero)
    assert all(checks.values()), checks
    certificate = {
        'jet_variable_order': names,
        'pfaffian': str(pf), 'P': [[str(v) for v in row] for row in P],
        'inverse_numerator_A': [[str(v) for v in row] for row in A],
        'dp': list(map(str, dp)), 'divergence_X': list(map(str, X)),
        'Jacobi_dual': list(map(str, dual)), 'E_dual': list(map(str, edual)),
        'formal_Bianchi_operator': str(B), 'formal_deRham_operator': str(deRham),
        'difference_Koszul_operator': str(koszul),
        'all_identity_normal_forms_zero': True,
    }
    (output / 'certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
    inputs = [
        'research/notes/S61_NONSQUARE_BIANCHI_WEIGHT4_CONTAINMENT_20261104_v63.md',
        'research/notes/S61_NONSQUARE_BIANCHI_GENERIC_COMBINATION_DIVISION_20261104_v63.md',
        'code/s61_nonsquare_full_bianchi_operator_20261104_v63.py',
        'research/notes/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW.md',
    ]
    report = {
        'id': STEM, 'all_passed': True,
        'method': 'one exact universal 4x4 first/second-jet polynomial test over QQ',
        'checks': checks, 'jet_variables': len(names),
        'no_new_source_points_fields_GB_or_full135_expansion': True,
        'new_primitive_syzygy_from_this_operation': False,
        'general_nonmonomial_saturation_decided': False,
        'actual_weight4_transport_or_full239_rank_claimed': False,
        'input_sha256': {name: sha(ROOT / name) for name in inputs},
        'source_sha256': sha(Path(__file__)),
        'evidence_sha256': {'certificate.json': sha(output / 'certificate.json')},
        'elapsed_seconds': time.monotonic() - start,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'all_passed': True, 'checks': len(checks),
                      'seconds': report['elapsed_seconds'], 'max_rss_kib': report['max_rss_kib']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    main(args.out)
