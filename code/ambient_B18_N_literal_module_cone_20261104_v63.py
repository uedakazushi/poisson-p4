#!/usr/bin/env python3
"""Literal actual N/R0 cone, reusing the frozen coefficient-Hensel binding.

This constructs a resolution, not the N binding producer or new source jets.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import resource
import signal
import time

from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG, Quadratic
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors
from ambient_B18_T3_syzygy_bockstein_preflight_20261104_v63 import add, mul, scale, ONE, ZERO, var

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for q in iter(lambda: f.read(1048576), b''):
            h.update(q)
    return h.hexdigest()


def sparse_decode(col):
    return {int(i): {tuple(mon): int(c) for c, mon in enc} for i, enc in col}


def compose(left, col):
    out = {}
    for j, c in col.items():
        for i, p in left[j].items():
            out[i] = add(out.get(i, ZERO), mul(c, p))
    return {i: p for i, p in out.items() if p}


def pair_add(a, b):
    return (add(a[0], b[0]), add(a[1], b[1]))


def pair_mul(a, b):
    return (add(mul(a[0], b[0]), mul(var(8), mul(a[1], b[1]))),
            add(mul(a[0], b[1]), mul(a[1], b[0])))


def pair_poly(p):
    """Substitute v=z1*tau, tau^2=nu in universal coefficient polys."""
    variables = [(var(i), ZERO) for i in range(9)] + [(ZERO, var(0))]
    ans = (ZERO, ZERO)
    for mon, c in p.items():
        q = ({(): c}, ZERO)
        for j in mon:
            q = pair_mul(q, variables[j])
        ans = pair_add(ans, q)
    return ans


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    np = ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    nr = np.parent/'report.json'
    up = ROOT/'reports/AMBIENT_B18_T3_SYZYGY_BOCKSTEIN_PREFLIGHT_20261104_v63/universal-syzygies-and-method.json'
    ur = up.parent/'report.json'
    d, universal = read(np), read(up)
    assert read(nr)['all_passed'] and read(ur)['all_passed']
    nodes = [tuple(n) for n in d['DAG_scalar_nodes']]
    D = DAG()
    D.nodes = nodes[:]
    D.cache = {n: i for i, n in enumerate(nodes)}
    assert len(D.cache) == len(nodes)
    D.values = []
    old_inverse_origins = []
    for i, n in enumerate(nodes):
        o = n[0]
        if o == 'c': q = F(n[1], n[2])
        elif o == 'v': q = F(0)
        elif o == 'a': q = D.values[n[1]]+D.values[n[2]]
        elif o == 'n': q = -D.values[n[1]]
        elif o == 'm': q = D.values[n[1]]*D.values[n[2]]
        elif o == 'i':
            assert D.values[n[1]]
            q = 1/D.values[n[1]]
            old_inverse_origins.append([i, n[1], str(D.values[n[1]])])
        else: raise ValueError(o)
        D.values.append(q)
    assert D.zero == 0 and D.one == 1

    # Recover all even metadata from the existing exact quadratic source, not
    # from a jet.  The 56 coefficients, including B-centers, are already fixed.
    Q = Quadratic(D, 'N_EXACT_FOLD')
    assert [Q.S, Q.T] == d['quadratic_trace_S_T_nodes']
    unknown = d['normalization_source28_variable_nodes']
    pnode = d['oldN_source11_parameter_values']['p']
    others = [i for i in unknown if i != pnode]
    us = d['all56_coefficient_unknown_nodes']
    assert len(others) == 27 and len(us) == 56
    kappa_position = d['mixed_jet_monomials'].index([7])
    kernel = {i: F(q[kappa_position]) for i, q in d['all28_mixed_source2jets']}
    pairs = {pnode: Q.p}
    for j, i in enumerate(others):
        pairs[i] = (us[j], D.add(us[27+j], D.const(kernel[i])))
    roots = d['physical150_N_source_nodes'] + d['actual_same_R0_readback_nodes9']
    qmap = {}
    for i in sorted(ancestors(nodes, roots)):
        n = nodes[i]
        o = n[0]
        if o == 'c': q = Q.scalar(D.const(F(n[1], n[2])))
        elif o == 'v':
            assert i in pairs, (i, n)
            q = pairs[i]
        elif o == 'a': q = Q.add(qmap[n[1]], qmap[n[2]])
        elif o == 'n': q = Q.neg(qmap[n[1]])
        elif o == 'm': q = Q.mul(qmap[n[1]], qmap[n[2]])
        elif o == 'i': q = Q.inv(qmap[n[1]])
        else: raise ValueError(o)
        qmap[i] = q
    assert [list(qmap[i]) for i in d['actual_same_R0_readback_nodes9']] == d['all9_R0_quadratic_p_readback_pairs']
    p150 = [qmap[i] for i in d['physical150_N_source_nodes']]
    trace150 = [Q.tracepair(q) for q in p150]
    selected = d['chosen8_odd_target_physical_rows']
    assert selected == [61, 18, 19, 28, 29, 45, 46, 60]
    z = [trace150[i][1] for i in selected]
    assert z == d['chosen8_odd_coefficient_nodes']
    aa = trace150[61][0]
    assert aa == d['ninth_actual_even_a_node']
    nu = d['nu_tracefree_square_node']
    assert nu == D.add(Q.T, D.mul(D.const(F(1, 4)), D.mul(Q.S, Q.S)))
    assert all(D.values[b] == 0 for a, b in trace150)

    # The ninth actual R0 formal coordinate remains independent in the base.
    rho61 = D.var('N_LITERAL_R0_rho61')
    vv = D.add(rho61, D.neg(aa))
    actual = z + [nu, vv]
    u = universal['conditional_small_N_cone']
    n1 = [sparse_decode(c) for c in u['small_d1_8_by36']]
    n2 = [sparse_decode(c) for c in u['small_d2_36_by84']]
    assert len(n1) == 36 and len(n2) == 84
    assert all(not compose(n1, c) for c in n2)
    generic_generators = [(ONE, ZERO)] + [(ZERO, var(i)) for i in range(1, 8)]
    for col in n1:
        value = (ZERO, ZERO)
        for i, p in col.items():
            value = pair_add(value, pair_mul(pair_poly(p), generic_generators[i]))
        assert value == (ZERO, ZERO)

    poly_memo = {}
    def at(p):
        key = tuple(sorted(p.items()))
        if key in poly_memo: return poly_memo[key]
        terms = []
        for mon, c in sorted(p.items()):
            value = D.const(c)
            for j in mon: value = D.mul(value, actual[j])
            terms.append(value)
        value = D.total(terms)
        poly_memo[key] = value
        return value
    def at_matrix(cols):
        return [{str(i): at(p) for i, p in col.items()} for col in cols]
    d1 = at_matrix(n1)
    d2 = at_matrix(n2)
    assert all(D.values[i] == 0 for c in d1+d2 for i in c.values())

    rel_labels = []
    for i in range(8):
        for j in range(i+1, 8): rel_labels.append(f'Koszul_z{i+1}_z{j+1}_after_E1_equals_v')
    rel_labels += [f'v_action_on_original_z{i+1}_tau' for i in range(8)]
    assert len(rel_labels) == 36
    even_rows = [108, 109, 112, 113, 117, 118, 119, 27]
    cohen = d['actual_Cohen_R8_parameter_nodes']
    R0map = {str(row): node for row, node in zip(even_rows, cohen)}
    R0map['61'] = rho61
    data = {
        'DAG_scalar_nodes': D.nodes,
        'actual_N_same_R0_binding_input': str(np.relative_to(ROOT)),
        'same_actual_R0_rows9': d['actual_same_R0_rows9'],
        'same_actual_R0_coordinate_nodes_by_row': R0map,
        'actual_R8_physical_rows8': even_rows,
        'actual_R8_parameter_nodes8': cohen,
        'all56_coefficient_unknown_nodes': us,
        'all56_actual_coefficient_Hensel_equations': d['all56_actual_coefficient_Hensel_equations'],
        'unit56_minor': d['unit56_minor'],
        'all150_actual_quadratic_p_pairs': p150,
        'all150_actual_tracefree_pairs_even_odd': trace150,
        'actual_tau_definition': {'p_source_node': pnode, 'S_node': Q.S, 'formula': 'tau=p-S/2'},
        'nu_actual_node': nu,
        'z1_through_z8_actual_odd_coefficient_nodes': z,
        'z1_through_z8_physical_target_rows': selected,
        'unit8_odd_minor': d['chosen8_odd_linear_determinant'],
        'actual_a61_node': aa,
        'actual_independent_R0_rho61_node': rho61,
        'actual_v_node': vv,
        'eight_N_module_generators': [
            {'index': 0, 'lift_in_any_actual_target_algebra': '1', 'N_trace_pair': [D.one, D.zero]}
        ] + [
            {'index': j, 'physical_target_row': row, 'actual_even_part_node': trace150[row][0],
             'lift_in_any_actual_target_algebra': f'rho{row}-a{row}(R8)',
             'N_trace_pair': [D.zero, z[j]]}
            for j, row in enumerate(selected[1:], 1)
        ],
        'R0_sparse_d1_8_by36_columns': d1,
        'R0_d1_column_labels': rel_labels,
        'R0_sparse_d2_36_by84_columns': d2,
        'D2_column_order': 'first56 lex triples of original z1..z8; last28 lex pairs in cone shifted Koszul2',
        'coefficient_polynomial_variables_in_order': ['z1','z2','z3','z4','z5','z6','z7','z8','nu','v'],
        'd1_universal_sparse_polynomial_columns': u['small_d1_8_by36'],
        'd2_universal_sparse_polynomial_columns': u['small_d2_36_by84'],
        'N_full_minimal_free_resolution_ranks': [8,36,84,126,126,84,36,9,1],
        'source_and_target_traces_not_identified_with_CA_or_CB': True,
        'D4_extension_not_yet_constructed': True,
        'T3_nonzero_depth_or_splitting_not_claimed': True,
    }
    evidence = out/'actual-N-literal-cone-DAG.json'
    evidence.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    assert evidence.stat().st_size < 32*1024*1024
    paths = [Path(__file__), ROOT/'research/notes'/f'{STEM}.md', np, nr, up, ur,
             ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py',
             ROOT/'code/ambient_B18_top_graph_telescoping_20261104_v63.py',
             ROOT/'code/ambient_B18_T3_syzygy_bockstein_preflight_20261104_v63.py']
    report = {
        'id': STEM, 'all_passed': True,
        'status': 'actual_N_same_R0_literal_8_36_84_resolution_constructed',
        'actual_N_binding_reused_not_independently_reproved_in_this_constructor': True,
        'all150_exact_quadratic_source_pairs_recovered': True,
        'all9_stored_R0_quadratic_pairs_and_all8_selected_b_nodes_match_exactly': True,
        'all36_relations_vanish_in_actual_bound_quadratic_algebra_by_polynomial_identity': True,
        'all84_d2_columns_are_syzygies_by_full_polynomial_coefficient_identity': True,
        'full_resolution_exact_by_Koszul_and_monic_mapping_cone': True,
        'module_generators': 8, 'first_relations': 36, 'second_syzygy_generators': 84,
        'D4_nil_extension_horseshoe_or_T3_presentation_computed': False,
        'no_new_jets_points_GB_primes_or_coefficient_expansion': True,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'evidence_sha256': {evidence.name: sha(evidence)},
        'DAG_nodes': len(D.nodes),
        'seconds': time.monotonic()-started,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'certificate_bytes': evidence.stat().st_size,
    }
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['all_passed','module_generators','first_relations','second_syzygy_generators','seconds','max_rss_kib','certificate_bytes']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)
