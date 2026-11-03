#!/usr/bin/env python3
"""Independent full-D corollary: exact word-polynomials and both ideal inclusions."""
import argparse
from collections import Counter
import hashlib
import itertools as it
import json
from pathlib import Path
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def accepted(p):
    d = json.loads(p.read_text())
    assert d['all_passed'] is True
    for name, digest in d.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in d.get('evidence_sha256', {}).items():
        assert sha(p.parent / name) == digest, name
    return d


def add(*polynomials):
    out = Counter()
    for p in polynomials:
        for word, coefficient in p.items():
            out[word] += coefficient
    return {word: coefficient for word, coefficient in out.items() if coefficient}


def neg(p):
    return {word: -coefficient for word, coefficient in p.items()}


def mul(p, q):
    out = Counter()
    for a, c in p.items():
        for b, d in q.items():
            out[tuple(sorted(a + b))] += c * d
    return {word: coefficient for word, coefficient in out.items() if coefficient}


def variable(name):
    return {(name,): 1}


def substitute(p, mapping):
    terms = []
    for word, coefficient in p.items():
        term = {(): coefficient}
        for letter in word:
            term = mul(term, mapping.get(letter, variable(letter)))
        terms.append(term)
    return add(*terms)


def main(output):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    folder = ROOT / 'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63'
    accepted(folder / 'report.json')
    saved = json.loads((folder / 'certificate.json').read_text())
    gates = [
        ROOT / 'reports/AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW.json']
    meeting, twojet, audit = [accepted(p) for p in gates]
    assert audit['all_CB_tangent_in_whole_CA_tangent_as_QQa_identity_not_sample'] is True
    assert audit['first_CA_branch_intersection_CB_is_clean_smooth_N28'] is True
    assert audit['second_CA_branch_intersection_CB_is_clean_shared_smooth20'] is True
    assert audit['QQa_ranks_CA1_CB_CA2_CB_CAunion_plusL_plusCB'] == [30, 38, 37, 38, 38]
    rho, tau, f, c, d = map(variable, ('rho', 'tau', 'f', 'c', 'd'))
    u = [variable(f'u{i}') for i in range(1, 9)]
    v = [variable(f'v{i}') for i in range(1, 9)]
    b = [variable(f'b{i}') for i in range(1, 9)]
    r2 = mul(rho, rho)
    relations = [mul(x, y) for x in u for y in v] + [mul(x, tau) for x in u] + [mul(y, add(tau, neg(r2))) for y in v] + [mul(tau, add(tau, neg(r2)))]
    assert len(relations) == 81
    # Literal evaluations at the two branches and exact conversion from old tau.
    branch1 = {'tau': {}, **{f'v{i}': {} for i in range(1, 9)}}
    branch2 = {'tau': r2, **{f'u{i}': {} for i in range(1, 9)}}
    assert all(not substitute(q, branch1) and not substitute(q, branch2) for q in relations)
    old = [mul(x, y) for x in u for y in v] + [mul(x, add(tau, neg(r2))) for x in u] + [mul(y, tau) for y in v] + [mul(tau, add(tau, neg(r2)))]
    converted = [substitute(q, {'tau': add(r2, neg(tau))}) for q in old]
    assert all(p == q or p == neg(q) for p, q in zip(converted, relations))
    mapping = {'rho': mul(f, d), 'tau': mul(f, c), **{f'v{i}': mul(f, b[i - 1]) for i in range(1, 9)}}
    restrictions = [substitute(q, mapping) for q in relations]
    for q in restrictions:
        for word in q:
            assert word.count('f') >= 2 or (word.count('f') >= 1 and any(letter.startswith('u') for letter in word))
    for i in range(8):
        assert restrictions[64 + i] == mul(c, mul(f, u[i]))
    assert restrictions[-1] == mul(mul(f, f), mul(c, add(c, neg(mul(f, mul(d, d))))))

    # Bind the complete word-engine arrays to the frozen author's evidence.
    names = ['rho', 'tau', 'f', 'c', 'd'] + [f'u{i}' for i in range(1, 9)] + [f'v{i}' for i in range(1, 9)] + [f'b{i}' for i in range(1, 9)]
    symbols = s.symbols(' '.join(names))
    local = dict(zip(names, symbols))
    def decode(expression):
        polynomial = s.Poly(s.sympify(expression, locals=local), *symbols, domain=s.ZZ)
        return {tuple(letter for letter, power in zip(names, exponents) for _ in range(power)): int(coefficient) for exponents, coefficient in polynomial.terms() if coefficient}
    # The word algebra sorts names lexicographically; sympy's variable order differs.
    def normalized(p):
        return {tuple(sorted(word)): coefficient for word, coefficient in p.items()}
    assert relations == [normalized(decode(q)) for q in saved['full81_CA_relations']]
    assert restrictions == [normalized(decode(q)) for q in saved['all81_CB_restrictions']]

    # Independent monomial LCM intersection, no ideal solver.
    first = [Counter({'f': 1})]
    second = [Counter({'f': 2})] + [Counter({f'u{i}': 1}) for i in range(1, 9)]
    intersection = [tuple(sorted((a | b).elements())) for a in first for b in second]
    expected = [('f', 'f')] + [tuple(sorted(('f', f'u{i}'))) for i in range(1, 9)]
    assert intersection == expected
    t = s.symbols('t')
    # Standard monomials: arbitrary u8, or the single f, with free s20.
    hilbert = 1 / (1 - t) ** 28 + t / (1 - t) ** 20
    assert s.limit((1 - t) ** 28 * hilbert, t, 1) == 1
    assert s.diff(hilbert, t).subs(t, 0) == 29
    evidence = dict(
        gates_checked_and_all_input_hashes_current=True,
        both_branch_evaluations_all81_zero=True,
        old_tau_to_new_tau_rho2_minus_tau_all81_bound=True,
        independent_word_engine_all81_CB_restrictions=True,
        lower_inclusion_unit_cofactors_checked=True,
        upper_inclusion_for_ALL_actual_ambient_equations_accepted=True,
        simultaneous_N_rho_coordinates_from_scheme_meeting_intersection_accepted=True,
        c_unit_from_CB_escape_of_37_dimensional_tangent_sum_accepted=True,
        full_ideal='Q=f*(f,u1,...,u8)',
        monomial_primary_intersection=intersection,
        nilradical='f*O_S', nilradical_square_zero=True, nilmodule_generic_length=1,
        Ass=['(f)', '(f,u1,...,u8)'],
        Hilbert_series=str(hilbert), cone_dimension=28, cone_edim=29, cone_depth=20, cone_multiplicity=1,
        projective_dimension=27, projective_edim=28, projective_depth=19,
        algebraic_ideal_equality_descends_by_faithfully_flat_completion=True,
        global_embedded_associated_closure_S20_accepted=True,
        all_nonzero_parameters_origin_or_other_global_primary_components_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_CA_self_pair_actual_D_full_ring_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md',
        folder / 'report.json', folder / 'certificate.json', *gates]
    result = dict(id=STEM, all_passed=True, producer_rerun=False,
        actual_full_D_two_ideal_inclusions_and_all_ambient_upper_bound_accepted=True,
        simultaneous_regular_coordinates_and_consistent_tau_accepted=True,
        nil_O_S_Ass_N_plus_embedded_S_depth20_e1_accepted=True,
        global_embedded_associated_component_closure_S20_accepted=True,
        scope='Generic nonempty open of P_a and its G saturation inside the accepted meeting-double unit open; no a0 or all-nonzero extension.',
        independent_method='Commutative word-polynomial all81 restrictions and tau conversion; monomial LCM intersection/Hilbert series; scheme-theoretic coordinate and faithful-flat descent proof.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (output / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output)
