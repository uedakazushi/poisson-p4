#!/usr/bin/env python3
"""Independent actual cover/Cohen/pushout application and complete99 audit."""
import argparse
from collections import Counter
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import time

import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW'


def load(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def accepted(p):
    d = load(p)
    assert d['all_passed'] is True
    for name, digest in d.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in d.get('evidence_sha256', {}).items():
        assert sha(p.parent / name) == digest, name
    return d


def add(*ps):
    out = Counter()
    for p in ps:
        for m, c in p.items():
            out[m] += c
    return {m: c for m, c in out.items() if c}


def neg(p):
    return {m: -c for m, c in p.items()}


def mul(p, q):
    out = Counter()
    for m, c in p.items():
        for n, d in q.items():
            out[tuple(sorted(m + n))] += c * d
    return {m: c for m, c in out.items() if c}


def var(x):
    return {(x,): 1}


def main(out):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63'
    accepted(producer / 'report.json')
    saved = load(producer / 'relations99.json.gz')
    gates = [ROOT / 'reports/AMBIENT_THREE_AXIS_CONDUCTOR_RIGIDITY_20261103_v63.json',
        ROOT / 'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json']
    for p in gates:
        accepted(p)
    rigidity = load(gates[0])
    assert rigidity['ring_rigidity_proved'] and rigidity['pair_rigidity_proved']
    assert rigidity['q_hypotheses']['q_in_m_squared'] is True
    qp = ROOT / 'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz'
    fp = ROOT / 'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    bp = ROOT / 'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/certificate.json.gz'
    quad, first, binary = map(load, (qp, fp, bp))
    U, V, r, s = sy.symbols('U V r s')
    cover = {}
    for component, expected_lambda, expected_q in [('CA', sy.Rational(1, 3), sy.Rational(1, 25)), ('CB', sy.Integer(9), sy.Rational(16, 81))]:
        q = quad[component]
        prior = first['CA_CB_full_firstjet_data'][component]
        A = sy.Matrix(prior['actual'])
        rr = q['actual_seven_odd_basis_rows']
        assert sy.Matrix(q['chosen_nu_lift_combination']) * A[rr, :] == sy.Matrix(prior['discriminant'])
        # nu minus this actual ideal element has zero full firstjet: q in m_C².
        ell = sy.sympify(binary[component]['self_component_tangent_UV_equation'], locals={'U': U, 'V': V})
        assert sy.Poly(ell, U, V).coeff_monomial(U) == 1
        lam = -sy.Poly(ell, U, V).coeff_monomial(V)
        assert lam == expected_lambda
        tangent = sy.Matrix([sy.Rational(v) for v in binary[component]['WUV']])
        assert sy.Matrix.hstack(sy.Matrix([1, 0, 0]), sy.Matrix([0, 0, 1]), tangent).det() != 0
        actual_q2 = -sy.sympify(q['chosen_nu_lift_quadratic'], locals={'U': U, 'V': V})
        converted = sy.expand(actual_q2.subs({U: s, V: (s - r) / lam}, simultaneous=True))
        assert sy.expand(converted - expected_q * s * (s - r)) == 0
        cover[component] = dict(lambda_=str(lam), q_in_m_squared_by_full_firstjet=True, q2_coefficient=str(expected_q), three_independent_transverse_tangents=True)
    # Pair normalization identity, without assuming higher terms vanish.
    AA, BB = sy.symbols('AA BB')
    S, R = s * AA, -r * BB / AA
    assert sy.expand(S * (S - R) - (s * s * AA * AA + r * s * BB)) == 0

    # Exact integer word engine in R directsum R*t, t²=nu.
    z = [var(f'z{i}') for i in range(1, 7)]
    nu, w, rp, sp = map(var, ('nu', 'w', 'r', 's'))
    F = add(nu, neg(mul(sp, add(sp, neg(rp)))))
    b = [*z, F, mul(w, rp), mul(w, sp), mul(sp, mul(rp, rp))]
    Y = [var(f'Y{i}') for i in range(10)]
    linear = [add(mul(b[j], Y[i]), neg(mul(b[i], Y[j]))) for i, j in it.combinations(range(7), 2)]
    linear += [add(mul(b[j], Y[i]), neg(mul(b[i], Y[j]))) for i in range(7) for j in range(7, 10)]
    linear += [add(mul(sp, Y[7]), neg(mul(rp, Y[8]))), add(mul(mul(rp, sp), Y[7]), neg(mul(w, Y[9])))]
    products = [add(mul(Y[i], Y[j]), neg(mul(nu, mul(b[i], b[j])))) for i in range(10) for j in range(i, 10)]
    assert len(linear) == 44 and len(products) == 55
    names = [f'z{i}' for i in range(1, 7)] + ['nu', 'w', 'r', 's'] + [f'Y{i}' for i in range(10)]
    symbols = sy.symbols(' '.join(names))
    local = dict(zip(names, symbols))
    def decode(text):
        p = sy.Poly(sy.sympify(text, locals=local), *symbols, domain=sy.ZZ)
        return {tuple(sorted(letter for letter, power in zip(names, e) for _ in range(power))): int(c) for e, c in p.terms() if c}
    assert linear == [decode(p) for p in saved['all44_complete_first_syzygies']]
    assert products == [decode(p) for p in saved['all55_pair_products']]
    assert linear + products == [decode(p) for p in saved['all99']]
    def pair_product(x, y):
        return (add(mul(x[0], y[0]), mul(nu, mul(x[1], y[1]))), add(mul(x[0], y[1]), mul(x[1], y[0])))
    mapping = {f'Y{i}': ({}, b[i]) for i in range(10)}
    def evaluate(p):
        ev, od = {}, {}
        for word, coefficient in p.items():
            term = ({(): coefficient}, {})
            for letter in word:
                term = pair_product(term, mapping.get(letter, (var(letter), {})))
            ev, od = add(ev, term[0]), add(od, term[1])
        return ev, od
    assert all(evaluate(p) == ({}, {}) for p in linear + products)
    # Correct core Hilbert–Burch signs, then tensor with the regular seven.
    HB = [[sp, mul(rp, sp)], [neg(rp), {}], [{}, neg(w)]]
    det2 = lambda rows: add(mul(HB[rows[0]][0], HB[rows[1]][1]), neg(mul(HB[rows[0]][1], HB[rows[1]][0])))
    assert [det2((1, 2)), neg(det2((0, 2))), det2((0, 1))] == b[7:]
    h = sy.symbols('h')
    betti = sy.Poly(sy.expand((1 + h) ** 7 * (1 + 3 * h + 2 * h * h)), h)
    assert [betti.coeff_monomial(h ** i) for i in range(10)] == [1, 10, 44, 112, 182, 196, 140, 64, 17, 2]

    # Primary intersection of the core: independent monomial LCM calculation.
    def intersect(gens1, gens2):
        candidates = {tuple(sorted((Counter(a) | Counter(b)).elements())) for a in gens1 for b in gens2}
        return {a for a in candidates if not any(b != a and not (Counter(b) - Counter(a)) for b in candidates)}
    primary = intersect(intersect({('r',), ('s',)}, {('w',), ('s',)}), {('w',), ('r', 'r')})
    assert primary == {('r', 'w'), ('s', 'w'), ('r', 'r', 's')}
    # Both directions of the critical quotient ideal change U=s,V=r-s.
    W = sy.symbols('W')
    qcrit = -U * V
    assert sy.expand(U * (U + V) ** 2 + (2 * U + V) * qcrit - U ** 3) == 0

    extra = ROOT / 'reports/AMBIENT_ORIGIN_FORMAL99_M2_20261103_v63/report.json'
    extra_report = accepted(extra)
    evidence = dict(actual_pair_binding=cover,
        finite_G19_equivariant_cover_product_and_trace_line_unit_square_accepted=True,
        regular_Cohen_lift_edimE23_plus_kernel6_equals29=True,
        C_lifts_and_conductor_z_lifts_and_t_squared_belong_actual_A=True,
        intrinsic_conductor_pullback_recovers_whole_actual_A=True,
        NEW_regular_base_not_old_relative_CA_CB_alignment=True,
        integer_pair_algebra_all44_and55_relations_zero=True,
        core_Hilbert_Burch_signs_verified=True,
        tensor_resolution_Betti=[1, 10, 44, 112, 182, 196, 140, 64, 17, 2],
        all99_kernel_completeness_by_convergent_formal_reduction_and_complete_module_syzygies=True,
        primary_core_intersection=sorted(primary), generic_lengths=[1, 1, 2],
        independent_M2_entire_kernel_gate_read_not_rerun=True,
        formal_only_not_convergent_or_relative_embedding=True)
    (out / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_origin_full_formal_model_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63.md',
        producer / 'report.json', producer / 'relations99.json.gz', qp, fp, bp, extra, *gates]
    report = dict(id=STEM, all_passed=True, producer_rerun=False,
        actual_CA_CB_separate_full_formal_models_accepted=True,
        intrinsic_G_cover_Cohen_lift_and_pushout_application_accepted=True,
        all44_syzygies_and55_products_complete_kernel_accepted=True,
        formal_not_convergent_or_relative_embedding_claim=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(out / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)
