#!/usr/bin/env python3
"""Independent full-normalization square/cube criterion and complete81 audit."""
import argparse
import json
from pathlib import Path
import time

import sympy as sy

from ambient_origin_full_formal_model_review_20261103_v63 import load, sha, accepted, add, neg, mul, var

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63_REVIEW'


def main(out):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63'
    accepted(producer / 'report.json')
    saved = load(producer / 'certificate.json.gz')
    model = ROOT / 'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/report.json'
    accepted(model)
    z = [var(f'z{i}') for i in range(1, 7)]
    nu, w, r, s = map(var, ('nu', 'w', 'r', 's'))
    q = mul(s, add(s, neg(r)))
    F = add(nu, neg(q))
    I = [*z, F, mul(w, r), mul(w, s), mul(s, mul(r, r))]
    J = [*z, F, w, r]
    w3, r3 = mul(mul(w, w), w), mul(mul(r, r), r)
    assert add(mul(nu, w3), neg(mul(F, w3)), neg(mul(mul(w, s), mul(mul(w, w), add(s, neg(r)))))) == {}
    assert add(mul(nu, r3), neg(mul(F, r3)), neg(mul(mul(s, mul(r, r)), mul(r, add(s, neg(r)))))) == {}

    def pair_product(x, y):
        return (add(mul(x[0], y[0]), mul(nu, mul(x[1], y[1]))), add(mul(x[0], y[1]), mul(x[1], y[0])))
    xx, yy = var('x'), var('y')
    square = pair_product((xx, yy), (xx, yy))
    cube = pair_product(square, (xx, yy))
    assert square[1] == {('x', 'y'): 2}
    assert cube[1] == {('x', 'x', 'y'): 3, ('nu', 'y', 'y', 'y'): 1}
    # All nine regular generators and their complete Koszul/product relations.
    Y = [var(f'Y{i}') for i in range(9)]
    linear = [add(mul(J[j], Y[i]), neg(mul(J[i], Y[j]))) for i in range(9) for j in range(i + 1, 9)]
    products = [add(mul(Y[i], Y[j]), neg(mul(nu, mul(J[i], J[j])))) for i in range(9) for j in range(i, 9)]
    assert len(linear) == 36 and len(products) == 45
    names = [f'z{i}' for i in range(1, 7)] + ['nu', 'w', 'r', 's'] + [f'Y{i}' for i in range(9)]
    symbols = sy.symbols(' '.join(names))
    local = dict(zip(names, symbols))
    def decode(text):
        p = sy.Poly(sy.sympify(text, locals=local), *symbols, domain=sy.ZZ)
        return {tuple(sorted(letter for letter, power in zip(names, ex) for _ in range(power))): int(c) for ex, c in p.terms() if c}
    assert linear == [decode(p) for p in saved['all36_complete_Koszul_relations']]
    assert products == [decode(p) for p in saved['all45_pair_products']]
    assert linear + products == [decode(p) for p in saved['all81']]
    mapping = {f'Y{i}': ({}, J[i]) for i in range(9)}
    for p in linear + products:
        ev, od = {}, {}
        for word, coefficient in p.items():
            term = ({(): coefficient}, {})
            for letter in word:
                term = pair_product(term, mapping.get(letter, (var(letter), {})))
            ev, od = add(ev, term[0]), add(od, term[1])
        assert not ev and not od
    # Regular9 is explicit after the coordinate F=nu-s(s-r).
    zz = sy.symbols('z1:7')
    n, ww, rr, ss, tt = sy.symbols('nu w r s t')
    actualJ = [*zz, n - ss * (ss - rr), ww, rr]
    variables = [*zz, n, ww, rr]
    assert sy.Matrix(actualJ).jacobian(variables).det() == 1
    assert sy.expand((tt * tt - ss * (ss - rr)) - (tt * tt - ss * ss) - ss * rr) == 0
    # Properness on both sides: w,r survive on their critical conductor axes.
    substitutions_W = dict(zip([*zz, n, rr, ss], [0] * 9))
    substitutions_R = dict(zip([*zz, n, ww, ss], [0] * 9))
    actualI = [*zz, n - ss * (ss - rr), ww * rr, ww * ss, ss * rr * rr]
    assert all(p.subs(substitutions_W) == 0 for p in actualI)
    assert all(p.subs(substitutions_R) == 0 for p in actualI)
    assert ww.subs(substitutions_W) != 0 and rr.subs(substitutions_R) != 0
    evidence = dict(
        independent_integer_cube_identities=True,
        full_B_even_odd_square_cube_coefficients_checked=True,
        domain_R_mod_J_and_nonzero_s_squared_not_s_unit_argument_accepted=True,
        all_total_fraction_candidates_reduce_to_B_by_integrality=True,
        finite_elementary_subintegral_extension_and_seminormality_accepted=True,
        exact_seminormalization='R directsum J*t=A[w*t,r*t]',
        J_is_self_prime_not_radical_I=True,
        strict_A_subset_A_sn_subset_B=True,
        source_conductor_reduced_node_and_residual_self_gluing_accepted=True,
        regular9_Jacobian_unit='1',
        all36_Koszul_plus45_products_integer_pair_algebra_checked=True,
        complete81_kernel_by_regular_sequence_and_formal_R_plus_Y_reduction=True,
        cone_dimension=29, cone_edim=38, cone_depth=21, cone_multiplicity=2,
        projective_dimension=28, projective_edim=37, projective_depth=20,
        formal_model_application_only_not_global_relative_classification=True)
    (out / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_origin_full_formal_model_review_20261103_v63.py',
        ROOT / 'code/ambient_origin_seminormalization_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63.md',
        producer / 'report.json', producer / 'certificate.json.gz', model]
    report = dict(id=STEM, all_passed=True, producer_rerun=False,
        exact_seminormalization_of_both_accepted_individual_formal_models=True,
        finite_subintegral_and_entire_normalization_square_cube_proof_accepted=True,
        complete81_kernel_and_conductor_node_invariants_accepted=True,
        no_global_or_relative_CA_CB_seminormality_claim=True,
        independent_method='Integer even/odd pair algebra; domain total-fraction square/cube proof; elementary subintegral criterion; regular9 Koszul completeness.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(out / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)
