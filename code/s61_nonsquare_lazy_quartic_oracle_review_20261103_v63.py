#!/usr/bin/env python3
"""Independent original135/JB/DAG binding and two exact H² recursions."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.rings import ring
from s61_model import S61

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63_REVIEW'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63'
    packet, report = folder/'oracle.json', folder/'report.json'
    data, ar = read(packet), read(report)
    assert ar['all_exact_source_and_graph_checks_passed'] and sha(packet) == ar['oracle_sha256']
    chart_path = ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
    chart = read(chart_path)
    extra_path = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    heavy_path = ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json'
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    frac = lambda x: K.from_expr(s.sympify(x))
    record = data['arithmetic']
    atoms = []
    known = [a, 500*a+90*e*e+3*e-24,
             (50*a+18*e*e-9*e-2)*d+60*e-40, 100*a+18*e*e-9*e-2]
    for num, den in record['atoms']:
        load = lambda ts: R.from_dict({tuple(ex): s.QQ(c) for ex, c in ts})
        c = K(load(num))/K(load(den))
        dd = c.denom
        for unit in known:
            while dd.rem(unit) == 0:
                dd = dd.exquo(unit)
        assert dd and dd.is_ground
        atoms.append(c)
    assert len(atoms) == 3484
    atomid = {str(c): i for i, c in enumerate(atoms)}
    assert len(atomid) == len(atoms)
    T = [[atoms[i] for i in row] for row in data['T_atom_ids']]
    U = [[atoms[i] for i in row] for row in data['U_atom_ids']]
    selected = data['selected34_rows']
    local_weights, row_weights = data['local_coordinate_weights'], data['equation_weights']
    vw = record['variable_weights']
    assert vw == [1, 1, 1, 1, 1, 2, 1, 1]
    assert selected == chart['selected34_rows']
    assert all(T[i] == [frac(row[j]) for j in chart['normal34']] for i, row in enumerate(chart['T34']))
    assert U == [[frac(x) for x in row] for row in chart['U34']]
    assert all(not any(T[i]) and not any(U[i]) for i in [19, 23])
    assert all(1 <= row_weights[i] <= 4 for i in selected)
    # Rebuild the homogenized135 directly from the original61 monomial data.
    # Introduce t, make every q polynomial degree2, then set q_fixed=1.
    model = S61()
    weights = model.reduced['parameter_weights']
    two = [i for i, w in enumerate(weights) if w == 2]
    qindex = {i: j for j, i in enumerate(two)}
    beta = {int(k[1:])-1: frac(v) for k, v in chart['beta'].items()}
    fixed, cols = data['fixed_original_coordinate'], data['other_original_coordinate_order']
    assert cols == [i for i in range(45) if i != fixed]
    local_index = {q: j for j, q in enumerate(cols)}
    local_index[45] = 44
    equations = model.reduced['equations']+[read(extra_path)['universal_new_equation']]
    assert len(equations) == 135
    J, B = [], []
    for equation in equations:
        homogeneous = {}
        for mon, numerator, denominator in equation:
            c = K(s.QQ(numerator, denominator))
            qmon = []
            for i in mon:
                if weights[i] == 1:
                    c *= beta[i]
                else:
                    qmon.append(qindex[i])
            assert len(qmon) <= 2
            key = tuple(sorted(qmon+[45]*(2-len(qmon))))
            homogeneous[key] = homogeneous.get(key, K.zero)+c
        jj, bb = {}, {}
        for mon, c in homogeneous.items():
            if not c:
                continue
            nfixed = mon.count(fixed)
            assert nfixed != 2
            if nfixed == 1:
                q = mon[0] if mon[1] == fixed else mon[1]
                j = local_index[q]
                jj[j] = jj.get(j, K.zero)+c
            else:
                key = tuple(sorted(local_index[q] for q in mon))
                bb[key] = bb.get(key, K.zero)+c
        J.append({j: c for j, c in jj.items() if c})
        B.append({jk: c for jk, c in bb.items() if c})
    for i in range(135):
        assert J[i] == {j: atoms[c] for j, c in data['J_sparse_atom_ids'][i]}
        assert B[i] == {tuple(sorted([j, k])): atoms[c] for j, k, c in data['B_sparse_atom_ids'][i]}
        assert all(local_weights[j] == row_weights[i] for j in J[i])
        assert all(local_weights[j]+local_weights[k] == row_weights[i] for j, k in B[i])
        assert all(sum((c*T[j][k] for j, c in J[i].items()), K.zero) == 0 for k in range(8))
    P = [[sum((c*U[j][k] for j, c in row.items()), K.zero) for k in range(34)] for row in J]
    assert P == [[atoms[i] for i in row] for row in data['P_JU_atom_ids']]
    for r, i in enumerate(selected):
        assert P[i] == [K.one if r == k else K.zero for k in range(34)]
    for i in range(45):
        assert all(not c or local_weights[i] == vw[j] for j, c in enumerate(T[i]))
        assert all(not c or local_weights[i] == row_weights[selected[j]] for j, c in enumerate(U[i]))
    nodes = record['nodes']
    nw, nd = record['node_weights'], record['node_ordinary_degrees']
    assert len(nodes) == 1373 and nodes[0] == ['z']
    nodeid = {}
    for i, n in enumerate(nodes):
        typ = n[0]
        if typ == 'z':
            key, w, degree = ('z',), -1, -1
        elif typ == 'v':
            key, w, degree = tuple(n), vw[n[1]], 1
        elif typ == 'c':
            key, w, degree = tuple(n), 0, 0
        elif typ == 'm':
            assert n[1] < i and n[2] < i
            key, w = tuple(n), nw[n[1]]+nw[n[2]]
            degree = None if nd[n[1]] is None or nd[n[2]] is None else nd[n[1]]+nd[n[2]]
            assert w <= 4
        else:
            assert typ == 'l' and all(j < i and atoms[c] for c, j in n[1])
            key = ('l', tuple(map(tuple, n[1])))
            assert len({nw[j] for c, j in n[1]}) == 1
            w = nw[n[1][0][1]]
            ds = {nd[j] for c, j in n[1]}
            degree = next(iter(ds)) if len(ds) == 1 else None
        assert w == nw[i] and degree == nd[i]
        assert key not in nodeid
        nodeid[key] = i
    def lin(terms):
        c = {}
        for a, n in terms:
            if n and a:
                c[n] = c.get(n, K.zero)+a
        c = {n: a for n, a in c.items() if a}
        if not c:
            return 0
        if len(c) == 1:
            n, a = next(iter(c.items()))
            if a == 1:
                return n
        return nodeid[('l', tuple((atomid[str(a)], n) for n, a in sorted(c.items())))]
    def product(i, j):
        if not i or not j or nw[i]+nw[j] > 4:
            return 0
        return nodeid[('m', min(i, j), max(i, j))]
    E = [lin(zip(row, [nodeid[('v', j)] for j in range(8)])) for row in T]
    assert E == data['linear_graph_roots']
    def bilinear(X, Y):
        out = []
        for i, row in enumerate(B):
            terms = []
            if row_weights[i] <= 4:
                for (j, k), c in row.items():
                    terms += ([(c, product(X[j], Y[k]))] if j == k else
                              [(c/2, product(X[j], Y[k])), (c/2, product(X[k], Y[j]))])
            out.append(lin(terms))
        return out
    lift = lambda w: [lin(zip(row, w)) for row in U]
    EE = bilinear(E, E)
    w2 = [lin([(-K.one, EE[i])]) for i in selected]
    Z2 = lift(w2)
    EZ2 = bilinear(E, Z2)
    w3 = [lin([(-2*K.one, EZ2[i])]) for i in selected]
    Z3 = lift(w3)
    ZZ, EZ3 = bilinear(Z2, Z2), bilinear(E, Z3)
    b4 = [lin([(K.one, x), (2*K.one, y)]) for x, y in zip(ZZ, EZ3)]
    w4 = [lin([(-K.one, b4[i])]) for i in selected]
    Z4 = lift(w4)
    for key, value in [('implicit_w2_roots', w2), ('implicit_w3_roots', w3), ('implicit_w4_roots', w4),
                       ('local_Z2_roots', Z2), ('local_Z3_roots', Z3), ('local_Z4_roots', Z4)]:
        assert value == data[key]
    project = lambda b: [lin([(K.one, b[i])]+[(-c, b[selected[k]]) for k, c in enumerate(P[i]) if c]) for i in range(135)]
    ff = [project(EE), project([lin([(2*K.one, n)]) for n in EZ2]), project(b4)]
    assert ff == data['residual_by_ordinary_degree_roots']
    roots = [lin([(K.one, f[i]) for f in ff]) for i in range(135)]
    assert roots == data['residual_weight_le4_roots'] and all(roots[i] == 0 for i in selected)
    genuine = [i for i, w in enumerate(row_weights) if w == 4 and i not in selected]
    assert len(genuine) == 33 and genuine == data['genuine_quartic_original_rows']
    assert [roots[i] for i in genuine] == data['genuine_quartic_roots']
    # Independent one-variable H recursion: only selected rows and the two
    # requested actual rows are expanded in the core coefficient field.
    H = [row[5] for row in T]
    def HH(i):
        return sum((c*H[j]*H[k] for (j, k), c in B[i].items()), K.zero)
    selected_HH = [HH(i) for i in selected]
    x2 = [-sum((c*v for c, v in zip(row, selected_HH)), K.zero) for row in U]
    cache = {0: (K.zero, K.zero, K.zero)}
    def evalH(i):
        if i in cache:
            return cache[i]
        n = nodes[i]
        if n[0] == 'v':
            v = (K.zero, K.one if n[1] == 5 else K.zero, K.zero)
        elif n[0] == 'c':
            v = (atoms[n[1]], K.zero, K.zero)
        elif n[0] == 'l':
            v = tuple(sum((atoms[c]*evalH(j)[k] for c, j in n[1]), K.zero) for k in range(3))
        else:
            a0, b0 = evalH(n[1]), evalH(n[2])
            v = tuple(sum((a0[j]*b0[k-j] for j in range(k+1)), K.zero) for k in range(3))
        cache[i] = v
        return v
    heavy = read(heavy_path)
    checks = []
    for row, oldformula in zip([77, 131], heavy['coefficient_formulas']):
        actual = HH(row)+sum((c*x2[j] for j, c in J[row].items()), K.zero)
        assert actual == evalH(roots[row])[2] == frac(oldformula)
        assert all(evalH(f[row])[2] == 0 for f in ff[1:])
        checks.append(dict(row=row, exact_H2_coefficient=str(actual)))
    # Exercise the public coefficient API only after both independent paths.
    from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle
    api, api_data = Oracle.load(packet)
    before = len(api.nodes)
    for check in checks:
        index = api_data['genuine_quartic_original_rows'].index(check['row'])
        node = api.coefficient(api_data['genuine_quartic_roots'][index], (0, 0, 0, 0, 0, 2, 0, 0))
        assert frac(api.core_value(node).as_expr()) == frac(check['exact_H2_coefficient'])
    after = len(api.nodes)
    assert before == 1373 and after == 1407 and len(nodes) == 1373
    assert sha(packet) == ar['oracle_sha256']
    clarification = folder/'API_COUNT_CLARIFICATION.json'
    cr = read(clarification)
    assert cr['original_serialized_DAG_nodes'] == before and cr['DAG_nodes_after_two_coefficient_queries'] == after
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', report, packet, chart_path,
             extra_path, heavy_path, clarification, folder/'API_H2_CHECK.json',
             ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py', ROOT/'code/s61_model.py']
    result = dict(id=STEM, all_passed=True,
        independent_original135_homogenized_source_J_B_including_G6=True,
        all135_JT_zero_and_selected_JU_identity_exact=True,
        all_core3484_atoms_only_known_a_F_B_C_denominators=True,
        all_source_J_B_weights_and_T_U_weights_verified=True,
        all1373_node_weights_degrees_and_exact_graph_recurrences_verified=True,
        selected34_graph_exact_by_positive_weight_not_jet_extrapolation=True,
        all33_original_quartic_roots_bound_to_full_source=True,
        independent_H2_formal_recursion_and_independent_DAG_coefficient_evaluation=checks,
        public_API_two_queries_match_both_independent_paths=True,
        serialized_node_count=before, post_query_in_memory_count=after,
        alias_count_clarification_accepted_and_frozen_file_unchanged=True,
        no_new_core_point_field_full_quartic_matrix_rank64_or_ambient_claim=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths}, seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'serialized_node_count', 'post_query_in_memory_count', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())
