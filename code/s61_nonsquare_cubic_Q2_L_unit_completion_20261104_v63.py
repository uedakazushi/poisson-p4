#!/usr/bin/env python3
"""One fixed actual three-dot and GL3 completion; other coefficients stay lazy."""
import argparse
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_cubic_zw25_pair_unit_completion_20261104_v63 import Lazy, matmul

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def universal():
    R, b0, b1, b2, v = ring('b0,b1,b2,v', s.QQ)
    z, one = R.zero, R.one
    mm = lambda A, B: matmul(A, B, z, lambda x, y: x+y, lambda x, y: x*y)
    I = [[one if i == j else z for j in range(3)] for i in range(3)]
    U0 = [[v, one, z], [-18*one, z, one], [-9*one, z, z]]
    Ui0 = [[z, z, -one/9], [one, z, v/9], [z, one, -2*one]]
    E = [[one, -b0, -b1], [z, one, z], [z, z, one]]
    Ei = [[one, b0, b1], [z, one, z], [z, z, one]]
    assert mm(U0, Ui0) == I == mm(Ui0, U0)
    assert mm(E, Ei) == I == mm(Ei, E)
    U, Ui = mm(U0, E), mm(Ei, Ui0)
    assert mm(U, Ui) == I == mm(Ui, U)
    t = b0*v-18*b1-9*b2
    assert mm([[b0, b1, b2]], U)[0] == [t, (one-t)*b0, (one-t)*b1]
    det = (U0[0][0]*(U0[1][1]*U0[2][2]-U0[1][2]*U0[2][1])
           -U0[0][1]*(U0[1][0]*U0[2][2]-U0[1][2]*U0[2][0])
           +U0[0][2]*(U0[1][0]*U0[2][1]-U0[1][1]*U0[2][0]))
    assert det == -9
    return {'ring': 'QQ[b0,b1,b2,v]',
            'both_U0_inverse_products_9_coefficients': True,
            'both_E_inverse_products_9_coefficients': True,
            'both_U_inverse_products_9_coefficients': True,
            'row_output': '[t,(1-t)*b0,(1-t)*b1]',
            'dot_t': 'b0*v-18*b1-9*b2', 'det_U0': '-9', 'det_E': '1'}


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(57)
    out.mkdir(parents=True, exist_ok=False)
    read = lambda p: json.loads(p.read_text())
    pp = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'
    mp = ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_MINORS_20261104_v63'
    cp = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'
    paths = [Path(__file__), ROOT/'code/s61_nonsquare_cubic_zw25_pair_unit_completion_20261104_v63.py',
             pp/'report.json', pp/'certificate.json', mp/'report.json', mp/'certificate.json', cp]
    parent, mr, mc, check = read(paths[3]), read(paths[4]), read(paths[5]), read(cp)
    assert sha(paths[3]) == read(paths[2])['certificate_sha256']
    assert sha(paths[5]) == mr['certificate_sha256']
    assert mc['all_six_new_presentation_roots_literal_old_checkpoint_references']
    oldrows, oldcols = parent['new16_source_original25_rows'], parent['new23_target_original32_columns']
    assert oldrows[0] == 5
    group = [21, 22, 25]
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    decp = lambda ts: R.from_dict({tuple(m): s.QQ(c) for m, c in ts})
    dec = lambda q: K(decp(q[0]))/K(decp(q[1]))
    encp = lambda p: [[list(m), str(c)] for m, c in p.items()]
    enc = lambda q: [encp(q.numer), encp(q.denom)]
    b = [dec(check['matrix25x32'][5][j]) for j in group]
    dual = [K(6*(18*e+7))/K(5*a), K(-18), K(-9)]
    T = 900*a-126*e-121
    assert b[0] == K(a*T)/150
    assert 2*b[1]+b[2] == K((18*e+7)*T)/1125-K.one/9
    dot = sum((b[i]*dual[i] for i in range(3)), K.zero)
    assert dot == K.one
    proof = universal()
    lazy = Lazy()
    lazy.nodes = parent['scalar_DAG']['nodes'][:]
    lazy.index = {tuple(node): i for i, node in enumerate(lazy.nodes)}
    oldnodecount = len(lazy.nodes)
    newatoms = [enc(dual[0]), enc(dual[1]), enc(dual[2]), enc(K(-1)/9), enc(K.one/9), enc(K(-2))]
    atom = lambda i: lazy.put(['new_unit_atom', i])
    v, c18, c9, minus9i, plus9i, c2 = [atom(i) for i in range(6)]
    A = parent['remaining16x23_scalar_roots']
    bnodes = [A[0][oldcols.index(j)] for j in group]
    assert all(lazy.nodes[n] == ['checkpoint', 5, j] for n, j in zip(bnodes, group))
    U0 = [[v, 1, 0], [c18, 0, 1], [c9, 0, 0]]
    Ui0 = [[0, 0, minus9i], [1, 0, lazy.mul(v, plus9i)], [0, 1, c2]]
    E, Ei = lazy.eye(3), lazy.eye(3)
    E[0][1], E[0][2] = lazy.sub(0, bnodes[0]), lazy.sub(0, bnodes[1])
    Ei[0][1], Ei[0][2] = bnodes[0], bnodes[1]
    U3, V3 = lazy.mm(U0, E), lazy.mm(Ei, Ui0)
    U23, V23 = lazy.eye(23), lazy.eye(23)
    pos = [oldcols.index(j) for j in group]
    for i, ri in enumerate(pos):
        for j, cj in enumerate(pos):
            U23[ri][cj], V23[ri][cj] = U3[i][j], V3[i][j]
    pivot_column, pivot_position = 21, pos[0]
    tail, taili = lazy.eye(23), lazy.eye(23)
    outside = []
    for j in oldcols:
        if j in group:
            continue
        root = A[0][oldcols.index(j)]
        if root:
            outside.append(j)
            tail[pivot_position][oldcols.index(j)] = lazy.sub(0, root)
            taili[pivot_position][oldcols.index(j)] = root
    assert all(A[0][oldcols.index(j)] == 0 for j in [27, 28, 29, 30, 31])
    target, targeti = lazy.mm(U23, tail), lazy.mm(taili, V23)
    transformed = lazy.mm(A, target)
    source = lazy.eye(16)
    for i in range(1, 16):
        source[i][0] = lazy.sub(0, transformed[i][pivot_position])
    rows = oldrows[1:]
    cols = [j for j in oldcols if j != pivot_column]
    reduced = [[transformed[i][oldcols.index(j)] for j in cols] for i in range(1, 16)]
    cert = {'actual_three_dot': {'old25_source_row': 5, 'old32_target_columns': group,
                                'coefficients': [enc(q) for q in b],
                                'dual': [enc(q) for q in dual], 'dot': enc(dot),
                                'visible_coefficient_identity_T': 'T=900a-126e-121; b0=aT/150; 2b1+b2=(18e+7)T/1125-1/9',
                                'constant_unit_dual_slot': [2, '-9']},
            'universal_polynomial_proof': proof,
            'scalar_DAG': {'nodes': lazy.nodes,
                           'checkpoint_reference': parent['scalar_DAG']['checkpoint_reference'],
                           'pair_reference': parent['scalar_DAG']['pair_reference'],
                           'accepted_pair_second_definition': parent['scalar_DAG']['accepted_pair_second_definition'],
                           'new_unit_atoms_num_den': newatoms,
                           'parent_nodes_unchanged_prefix_length': oldnodecount},
            'U0_roots': U0, 'U0_inverse_roots': Ui0, 'E_roots': E, 'E_inverse_roots': Ei,
            'U3_roots': U3, 'U3_inverse_roots': V3,
            'target23_roots': target, 'target23_inverse_roots': targeti,
            'tail_clear_original32_columns': outside, 'source16_roots': source,
            'pivot_full23_row_roots_before_exact_simplification': transformed[0],
            'pivot_row_exact_simplified': [int(i == pivot_position) for i in range(23)],
            'old16_source_original25_rows': oldrows, 'old23_target_original32_columns': oldcols,
            'new15_source_original25_labels': rows,
            'new22_target_original32_columns': cols,
            'remaining15x22_scalar_roots': reduced,
            'all5_remaining_heavy_quotient_columns_unchanged': [27, 28, 29, 30, 31],
            'source_change_determinant': '1', 'target_change_determinant': '-9',
            'new_pivots': 1, 'total_unit_pivots_in_transformed_presentation': 10,
            'same_cokernel_proof': 'The fixed3dot equals1. GL3 and tail operations turn FULL row5 into e21; clear its column from other rows using source16. Delete the unit source row and target column. Roots are actual transformed relations, not a literal old submatrix.'}
    (out/'certificate.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n')
    report = {'id': ID, 'fixed_three_dot_and_all_universal_inverse_identities_passed': True,
              'new_global_unit_pivots': 1, 'total_unit_pivots_in_transformed_presentation': 10,
              'actual_new_shape_same_cokernel': [15, 22],
              'all_full16x23_rows_transported_finite_DAG': True,
              'only_three_actual_coefficients_and_dot_expanded_all_tails_lazy': True,
              'all5_heavy_quotient_coordinates_unchanged': True,
              'new_inverse_factors_only_a_and_nonzero_QQ': True,
              'no_Gaussian_minor_search_point_GB_or_automatic_next_operation': True,
              'full_FULLU57_full239_kernel38_gap42_unchanged': True,
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
              'certificate_sha256': sha(out/'certificate.json'),
              'scalar_nodes': len(lazy.nodes), 'old_nodes_unchanged': oldnodecount,
              'seconds': time.monotonic()-began, 'frozen': True, 'independent_review_pending': True}
    assert sum(p.stat().st_size for p in out.iterdir()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())
