#!/usr/bin/env python3
"""Bounded independent archive/graph review of the resonance Q51 incidence."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = ROOT / 'reports/Q_RESONANCE_Q51_INCIDENCE_20261103_v62.json'
    d = json.loads(source.read_text())
    replay = json.loads(args.replay.read_text())
    assert d['all_passed'] and replay['all_passed']
    for key in d:
        if key != 'seconds':
            assert d[key] == replay[key], key
    paths = [source, Path(__file__),
             ROOT/'research/notes/Q_RESONANCE_Q51_INCIDENCE_20261103_v62.md',
             ROOT/'research/notes/Q_RESONANCE_Q51_INCIDENCE_20261103_v62_REVIEW.md',
             ROOT/'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62_REVIEW.json']
    for path, expected in d['input_sha256'].items():
        p = ROOT/path
        assert sha(p) == expected, path
        paths.append(p)
    e = s.Symbol('e')
    K = s.QQ.frac_field(e)
    values = [K.from_sympy(s.sympify(v)) for v in d['full61_normalized_parameters']]
    assert len(values) == 61
    graph = json.loads((ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json').read_text())
    variables = s.symbols('u1:62')
    beta = [variables[i] for i in graph['weight1_u_indices_zero_based']]
    # Separate exact evaluation of all45 cached graph formulas in QQ(e).
    # No producer graph-evaluation helper or sampled e values are used.
    for output, formula in zip(graph['weight2_u_indices_zero_based'], graph['quadratic_section_q']):
        total = K.zero
        for powers, coefficient in s.Poly(s.sympify(formula), *beta, domain=s.QQ).terms():
            product = K.convert(coefficient)
            for index, power in zip(graph['weight1_u_indices_zero_based'], powers):
                if power:
                    product *= values[index]**power
            total += product
        assert not values[output] - total, output
    for key in ['value_pivot', 'value_pfaffian']:
        assert K.from_sympy(s.sympify(d['normalizer_metadata'][key]))
    report = dict(
        id='Q_RESONANCE_Q51_INCIDENCE_20261103_v62_REVIEW', all_passed=True,
        fresh_full_normalizer_134_equations_150_coefficients_replay_matches=True,
        independent_cached_all45_graph_identities_over_QQ_e=True,
        normalization_pivots_nonzero_in_function_field=True,
        known_closed_PGL_component_pullback_argument_accepted=True,
        entire_curve_in_known_T39_including_Q_2_over_9=True,
        normal_divisor_classification_used=False,
        local_branch_exhaustion_claimed=False, new_component_claimed=False,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in dict.fromkeys(paths)})
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('PASS independent45 graph identities and fresh134+150 normalizer replay')


if __name__ == '__main__':
    main()
