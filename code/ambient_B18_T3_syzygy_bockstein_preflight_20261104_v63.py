#!/usr/bin/env python3
"""Read-only, universal exact syzygies for the CB99/Bockstein route.

No source evaluation, jet extension, standard basis, or module membership job.
The 882 identities use an independent sparse integer polynomial implementation.
"""
import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_T3_SYZYGY_BOCKSTEIN_PREFLIGHT_20261104_v63'
NAMES = [f'z{i+1}' for i in range(7)] + ['r', 's', 'w', 'A', 'B', 'nu'] + [f'Z{i}' for i in range(10)]
ZERO = {}
ONE = {(): 1}


def var(i):
    return {(i,): 1}


def add(*args):
    out = Counter()
    for p in args:
        for k, v in p.items():
            out[k] += v
    return {k: v for k, v in out.items() if v}


def scale(p, c):
    return {k: v*c for k, v in p.items() if v*c}


def mul(a, b):
    out = Counter()
    for ka, va in a.items():
        for kb, vb in b.items():
            out[tuple(sorted(ka+kb))] += va*vb
    return {k: v for k, v in out.items() if v}


def put(out, i, p):
    out[i] = add(out.get(i, ZERO), p)
    if not out[i]:
        del out[i]


def enc(p):
    return [[c, list(mon)] for mon, c in sorted(p.items())]


def stream_sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    z = [var(i) for i in range(7)]
    r, s, w, aa, bb, nu = [var(i) for i in range(7, 13)]
    zz = [var(i) for i in range(13, 23)]
    g = add(mul(r, aa), mul(s, bb))
    core_b = [g, mul(w, r), mul(w, s)]
    b = z + core_b
    core_d2 = [{0: w, 1: scale(aa, -1), 2: scale(bb, -1)},
               {1: s, 2: scale(r, -1)}]

    # Total resolution Koszul(z7) tensor the Hilbert--Burch core resolution.
    def basis(n):
        ans = []
        for q, multiplicity in [(0, 1), (1, 3), (2, 2)]:
            p = n-q
            if 0 <= p <= 7:
                ans += [(ks, q, a) for ks in combinations(range(7), p) for a in range(multiplicity)]
        return ans

    bases = [basis(n) for n in range(10)]
    indices = [{x: i for i, x in enumerate(q)} for q in bases]

    def differential(n, label):
        ks, q, alpha = label
        p = len(ks)
        col = {}
        for t, j in enumerate(ks):
            put(col, indices[n-1][(ks[:t]+ks[t+1:], q, alpha)], scale(z[j], (-1)**t))
        if q == 1:
            put(col, indices[n-1][(ks, 0, 0)], scale(core_b[alpha], (-1)**p))
        if q == 2:
            for j, coeff in core_d2[alpha].items():
                put(col, indices[n-1][(ks, 1, j)], scale(coeff, (-1)**p))
        return col

    # basis(1) is exactly z1..z7,g,wr,ws.
    assert bases[1] == [((i,), 0, 0) for i in range(7)] + [((), 1, i) for i in range(3)]
    syz_b = [differential(2, q) for q in bases[2]]
    s_polys = [add(*(mul(c, zz[i]) for i, c in col.items())) for col in syz_b]
    assert len(syz_b) == 44
    assert all(not add(*(mul(c, b[i]) for i, c in col.items())) for col in syz_b)

    pairs = list(combinations_with_replacement(range(10), 2))
    pair_index = {p: 44+i for i, p in enumerate(pairs)}
    p_polys = [add(mul(zz[i], zz[j]), scale(mul(nu, mul(b[i], b[j])), -1)) for i, j in pairs]
    equations = s_polys + p_polys
    labels = [f'linear_{i}' for i in range(44)] + [f'product_{i}_{j}' for i, j in pairs]

    def pair_lift(i, j):
        """Coefficients in the 44 S's for b_i Z_j-b_j Z_i."""
        assert i < j
        if i < 7:
            label = ((i, j), 0, 0) if j < 7 else ((i,), 1, j-7)
            return {indices[2][label]: ONE}
        h0 = indices[2][((), 2, 0)]
        h1 = indices[2][((), 2, 1)]
        if (i, j) == (7, 8):
            return {h0: scale(r, -1), h1: bb}
        if (i, j) == (7, 9):
            return {h0: scale(s, -1), h1: scale(aa, -1)}
        assert (i, j) == (8, 9)
        return {h1: scale(w, -1)}

    pair_lifts = {}
    for i, j in combinations(range(10), 2):
        q = pair_lift(i, j)
        assert add(*(mul(c, s_polys[k]) for k, c in q.items())) == add(mul(b[i], zz[j]), scale(mul(b[j], zz[i]), -1))
        pair_lifts[f'{i},{j}'] = q

    columns = []
    kinds = []
    for label in bases[3]:
        columns.append(differential(3, label))
        kinds.append('Koszul_HB_second_syzygy')
    assert len(columns) == 112

    for alpha, rel in enumerate(syz_b):
        for k in range(10):
            col = {alpha: zz[k]}
            for i, coeff in rel.items():
                put(col, pair_index[tuple(sorted((i, k)))], scale(coeff, -1))
            columns.append(col)
            kinds.append('mixed_linear_product')

    for i, j, k in combinations_with_replacement(range(10), 3):
        if i == k:
            continue
        anchor = ((i, j), k)
        others = [((i, k), j, i, j, k), ((j, k), i, j, i, k)]
        seen = {anchor}
        for pair, rem, common, lo, hi in others:
            if (pair, rem) in seen:
                continue
            seen.add((pair, rem))
            col = {}
            put(col, pair_index[(i, j)], zz[k])
            put(col, pair_index[pair], scale(zz[rem], -1))
            for index, coeff in pair_lift(lo, hi).items():
                put(col, index, mul(mul(nu, b[common]), coeff))
            columns.append(col)
            kinds.append('triple_monomial_overlap')

    assert Counter(kinds) == {'Koszul_HB_second_syzygy': 112,
                              'mixed_linear_product': 440,
                              'triple_monomial_overlap': 330}
    assert len(columns) == 882
    for col in columns:
        assert not add(*(mul(coeff, equations[i]) for i, coeff in col.items()))

    # Homogeneity check in the actual relative weight (w and Z8,Z9 weight 1).
    weights = [0]*23
    for j in [9, 21, 22]:
        weights[j] = 1

    def weight(p):
        ws = {sum(weights[j] for j in mon) for mon in p}
        assert len(ws) == 1
        return next(iter(ws))

    eq_weights = [weight(p) for p in equations]
    col_weights = []
    for col in columns:
        ws = {weight(coeff)+eq_weights[i] for i, coeff in col.items()}
        assert len(ws) == 1
        col_weights.append(next(iter(ws)))
    assert Counter(eq_weights) == {0: 64, 1: 32, 2: 3}
    assert Counter(col_weights) == {0: 448, 1: 368, 2: 64, 3: 2}

    # Conditional SMALL N-module cone.  These are universal identities only;
    # the actual same-R0 Cohen-coordinate binding is a separate, pending gate.
    # Reuse independent polynomial indices for z1..z8,nu,v; their meanings here
    # are recorded separately from the CB99 variables above.
    nz = [var(i) for i in range(8)]
    nnu, vv = var(8), var(9)
    nk2 = list(combinations(range(8), 2))
    nk3 = list(combinations(range(8), 3))
    k2index = {q: i for i, q in enumerate(nk2)}
    nd1 = [{j+1: nz[i], i+1: scale(nz[j], -1)} for i, j in nk2]
    nd2 = []
    for ks in nk3:
        col = {}
        for t, j in enumerate(ks):
            put(col, k2index[ks[:t]+ks[t+1:]], scale(nz[j], (-1)**t))
        nd2.append(col)
    action0 = [{1: ONE}] + [{0: mul(nnu, mul(nz[0], nz[i]))} for i in range(8)]
    h0 = []
    for i, col in enumerate(action0):
        q = {i: vv}
        for j, coeff in col.items():
            put(q, j, scale(coeff, -1))
        h0.append(q)
    cone_d1 = nd1 + h0
    cone_d2 = [dict(c) for c in nd2]
    for i, rel in enumerate(nd1):
        col = {i: vv}
        for j, coeff in rel.items():
            put(col, 28+j, scale(coeff, -1))
        cone_d2.append(col)
    assert len(cone_d1) == 37 and len(cone_d2) == 84

    def compose_cols(left, right):
        ans = {}
        for j, coeff in right.items():
            for i, p in left[j].items():
                put(ans, i, mul(coeff, p))
        return ans

    assert all(not compose_cols(cone_d1, c) for c in cone_d2)
    # Unit -1 in row E1 / column shifted E0; exact Schur cancellation.
    pivot_position = 28
    pivot = cone_d1[pivot_position]
    assert pivot[1] == scale(ONE, -1)
    keep_rows = [i for i in range(9) if i != 1]
    keep_columns = [i for i in range(37) if i != pivot_position]
    def nproject(col):
        return {i: p for i, row in enumerate(keep_rows)
                if (p := add(col.get(row, ZERO), mul(pivot.get(row, ZERO), col.get(1, ZERO))))}
    nsmall_d1 = [nproject(cone_d1[j]) for j in keep_columns]
    # The transformed pivot row of d2 vanishes by the exact chain identity.
    for col in cone_d2:
        pivot_new = add(col.get(pivot_position, ZERO),
                        *(scale(mul(cone_d1[j].get(1, ZERO), col.get(j, ZERO)), -1)
                          for j in keep_columns))
        assert not pivot_new
    nsmall_d2 = [{i: col[j] for i, j in enumerate(keep_columns) if j in col} for col in cone_d2]
    assert all(not compose_cols(nsmall_d1, c) for c in nsmall_d2)

    source_paths = [
        ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',
        ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/report.json',
    ]
    actual = json.loads(source_paths[0].read_text())
    assert actual['CB_HB_matrix'] == '[[w,0],[-A,s],[-B,-r]], minors(g,-wr,ws)'
    assert Counter(q['weight'] for q in actual['all99_CB_kernel_and87_graph_generators_evaluated_on_CA']) == {0: 68, 1: 85, 2: 27, 3: 6}
    # Only input labels/counts are compared here: actual arithmetic was already reviewed.
    for path in source_paths[1:]:
        assert json.loads(path.read_text())['all_passed']
    hp = ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    honest = json.loads(hp.read_text())
    variable_names = [n[1] for n in honest['DAG_scalar_nodes'] if n[0] == 'v']
    source_counts = {
        'physical_base_parameters': 9,
        'original_unit_source_unknowns': 96,
        'pointed_unit57_instances': len(actual['all16_pointed_unit57_inverse_instances']) + len(honest['all8_pure_CA_prefix57_unit_systems']),
        'scalar_variable_nodes_including_unused_aliases': len(variable_names),
        'distinct_scalar_variable_names': len(set(variable_names)),
        'regular_coordinate_quotient_nodes': sum(n[0] == 'r' for n in honest['DAG_scalar_nodes']),
        'regular_noncoordinate_g_quotients': sum(n[0] == 'g' for n in honest['DAG_scalar_nodes']),
    }
    assert source_counts['pointed_unit57_instances'] == 24
    assert source_counts['distinct_scalar_variable_names'] == 9+96+24*57 == 1473
    assert source_counts['regular_coordinate_quotient_nodes'] == 1576
    assert source_counts['regular_noncoordinate_g_quotients'] == 4

    ev = {
        'variables': NAMES,
        'polynomial_encoding': 'integer coefficient, sorted variable-index multiset',
        'generic_b10': [enc(p) for p in b],
        'generic99_labels': labels,
        'generic99_equations': [enc(p) for p in equations],
        'all44_first_syzygies_b10': [[[i, enc(p)] for i, p in sorted(c.items())] for c in syz_b],
        'all45_pair_syzygy_lifts': {k: [[i, enc(p)] for i, p in sorted(c.items())] for k, c in pair_lifts.items()},
        'all882_syzygies_of_99': [{'kind': kind, 'weight': wt,
                                 'entries': [[i, enc(p)] for i, p in sorted(col.items())]}
                                for kind, wt, col in zip(kinds, col_weights, columns)],
        'R_over_I_resolution_ranks': [len(q) for q in bases],
        'ideal_I_resolution_ranks': [len(q) for q in bases[1:]],
        'Bockstein_formula': 'T3=coker[H1(C)_4 -> H1(C/w)_4], C=F_CB tensor_P A_CA',
        'connecting_map': 'cycle c mod w -> class of d1(c)/w in D3',
        'H1_of_evaluated_resolution_not_assumed_zero': True,
        'conditional_small_N_cone': {
            'variables_here': ['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'nu', 'v'],
            'hypothesis_not_yet_actual_binding': 'same R0=R8[[v]], N=R8+m8*tau, tau^2=nu, v=z1*tau',
            'T0_columns': [[[i, enc(p)] for i, p in sorted(c.items())] for c in action0],
            'positive_resolution_actions_Tj': 'zero for all j>=1',
            'cone_d1_9_by37': [[[i, enc(p)] for i, p in sorted(c.items())] for c in cone_d1],
            'cone_d2_37_by84': [[[i, enc(p)] for i, p in sorted(c.items())] for c in cone_d2],
            'unit_pivot': {'row': 1, 'column': 28, 'value': -1},
            'small_d1_8_by36': [[[i, enc(p)] for i, p in sorted(c.items())] for c in nsmall_d1],
            'small_d2_36_by84': [[[i, enc(p)] for i, p in sorted(c.items())] for c in nsmall_d2],
            'actual_N_Cohen_binding_computed': False,
            'with_actual_regular_nilcurve_horseshoe_free_ranks_through2': [9, 44, 112],
        },
        'actual_coefficient_Hensel_support_counts': source_counts,
        'no_claim_actual_T3_presentation_or_nonzero': True,
    }
    evidence = out/'universal-syzygies-and-method.json'
    evidence.write_text(json.dumps(ev, separators=(',', ':'))+'\n')
    assert evidence.stat().st_size < 32*1024*1024
    report = {
        'id': STEM, 'all_passed': True,
        'status': 'constructive_complete_component_syzygy_prefix_and_exact_Bockstein_method_not_T3_solution',
        'independent_sparse_integer_polynomial_identities': 882,
        'component_first_syzygies_complete_by_filtered_reduction_proof': True,
        'preserves_all_z7_syzygies': True,
        'complete_first_syzygy_counts': dict(Counter(kinds)),
        'component_d2_weight_counts': dict(Counter(col_weights)),
        'component_plus87graphs_second_free_rank_if_expanded': 882+99*87+87*86//2,
        'graph_Koszul_factor_should_remain_unexpanded': True,
        'conditional_N_cone_exact_unit_cancel_and_chain_identities_verified': True,
        'conditional_D4_horseshoe_free_ranks_through2': [9, 44, 112],
        'actual_N_same_R0_coordinate_gate_still_pending': True,
        'actual_T3_homology_computed': False,
        'actual_T3_nonzero_or_depth_claimed': False,
        'unit_reduction_review_reused_without_rerun': True,
        'HONEST9_LOW2_read_only_inconclusive_not_independently_replayed': True,
        'no_GB_new_jets_new_points_new_primes_or_source_rerun': True,
        'input_sha256': {str(p.relative_to(ROOT)): stream_sha(p) for p in [
            Path(__file__), ROOT/'research/notes'/f'{STEM}.md', *source_paths, hp]},
        'evidence_sha256': {evidence.name: stream_sha(evidence)},
        'seconds': time.monotonic()-started,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'certificate_bytes': evidence.stat().st_size,
    }
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['all_passed', 'independent_sparse_integer_polynomial_identities', 'seconds', 'max_rss_kib', 'certificate_bytes']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)
