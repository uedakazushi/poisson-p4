#!/usr/bin/env python3
"""One fixed six-dot identity and global GL6 completion, without Gaussian."""
import argparse
import hashlib
import itertools
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def matmul(A, B, zero, add, mul):
    return [[sumprod([(A[i][k], B[k][j]) for k in range(len(B))], zero, add, mul)
             for j in range(len(B[0]))] for i in range(len(A))]


def sumprod(terms, zero, add, mul):
    value = zero
    for a, b in terms:
        value = add(value, mul(a, b))
    return value


class Lazy:
    def __init__(self):
        self.nodes = [['z'], ['one']]
        self.index = {('z',): 0, ('one',): 1}

    def put(self, node):
        key = tuple(node)
        if key not in self.index:
            assert len(self.nodes) < 200000
            self.index[key] = len(self.nodes)
            self.nodes.append(node)
        return self.index[key]

    def add(self, a, b):
        if not a:
            return b
        if not b:
            return a
        return self.put(['a', a, b])

    def sub(self, a, b):
        if not b:
            return a
        if a == b:
            return 0
        return self.put(['s', a, b])

    def mul(self, a, b):
        if not a or not b:
            return 0
        if a == 1:
            return b
        if b == 1:
            return a
        return self.put(['m', a, b])

    def mm(self, A, B):
        return matmul(A, B, 0, self.add, self.mul)

    def eye(self, n):
        return [[int(i == j) for j in range(n)] for i in range(n)]


def universal_completion():
    """Integer polynomial identities in 11 indeterminates; no actual minors."""
    rr = ring(','.join([f'b{i}' for i in range(5)]+[f'v{i}' for i in range(5)]+['k']), s.ZZ)
    R, xx = rr[0], rr[1:]
    b, v, k = list(xx[:5]), list(xx[5:10]), xx[10]
    zero, one = R.zero, R.one
    mm = lambda A, B: matmul(A, B, zero, lambda x, y: x+y, lambda x, y: x*y)
    eye = lambda n: [[one if i == j else zero for j in range(n)] for i in range(n)]
    t = sum((b[i]*v[i] for i in range(5)), zero)
    U0 = [[(one if i == j else zero)-v[i]*b[j] for j in range(5)]+[v[i]] for i in range(5)]
    U0 += [b+[zero]]
    E, Ei = eye(6), eye(6)
    for j in range(5):
        E[5][j], Ei[5][j] = -k*b[j], k*b[j]
    W = [[v[i]*b[j] for j in range(5)]+[-v[i]] for i in range(5)]
    W += [[-x for x in b]+[one]]
    U, Ui = mm(U0, E), mm(Ei, U0)
    U02, UV, VU, EiWE = mm(U0, U0), mm(U, Ui), mm(Ui, U), mm(mm(Ei, W), E)
    assert mm(E, Ei) == eye(6) == mm(Ei, E)
    for i in range(6):
        for j in range(6):
            ij = one if i == j else zero
            assert U02[i][j]-ij == (t-one)*W[i][j]
            assert UV[i][j]-ij == (t-one)*W[i][j]
            assert VU[i][j]-ij == (t-one)*EiWE[i][j]
    row = mm([b+[k]], U)[0]
    assert row == [(one-t)*(one+k)*x for x in b]+[t]
    # A fixed universal 6x6 Leibniz identity, not a core-parameter minor.
    determinant = zero
    for p in itertools.permutations(range(6)):
        term = one
        for i in range(6):
            term *= U0[i][p[i]]
        if sum(p[i] > p[j] for i in range(6) for j in range(i+1, 6)) % 2:
            determinant -= term
        else:
            determinant += term
    assert determinant == -t
    return {'integer_polynomial_base_variables': [str(x) for x in xx],
            'dot_symbol': 't=sum(b_i*v_i)',
            'all_36_U0_square_remainders': 'U0^2-I=(t-1)*[[v*b,-v],[-b,1]]',
            'all_36_right_inverse_remainders': 'U*Ui-I=(t-1)*W',
            'all_36_left_inverse_remainders': 'Ui*U-I=(t-1)*Ei*W*E',
            'both_E_inverse_products_exact': True,
            'six_row_output_exact': '[(1-t)*(1+k)*b,t]',
            'universal_720_term_determinant_identity': 'det(U0)=-t',
            'at_actual_t1': 'U0 is an involution, Ui=Ei*U0, det(U)=-1, (b,k)*U=(0,1)'}


def main(out):
    began = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    signal.alarm(57)
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'certificate.json').exists() and not (out/'report.json').exists()
    assert all(p.name == 'failed-attempt' for p in out.iterdir())
    read = lambda p: json.loads(p.read_text())
    old = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    pd = ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63'
    paths = {'checkpoint': old/'checkpoint.json', 'unit_report': old/'report.json',
             'unit_review': ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json',
             'pair_certificate': pd/'certificate.json', 'pair_report': pd/'report.json',
             'input': ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json',
             'P22_dual': ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63/global-P22.json',
             'P22_review': ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW/report.json',
             'support': ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_Q3_SUPPORT_PREFLIGHT_20261104_v63/report.json',
             'source': Path(__file__)}
    cp, ur, pc = read(paths['checkpoint']), read(paths['unit_report']), read(paths['pair_certificate'])
    assert cp['completed_pivots'] == 8 and ur['frozen']
    assert sha(paths['checkpoint']) == ur['evidence_sha256']['checkpoint.json']
    assert sha(paths['pair_certificate']) == read(paths['pair_report'])['certificate_sha256']
    assert sha(paths['input']) == cp['source_input_sha256']
    assert read(paths['unit_review'])['all_passed'] and read(paths['P22_review'])['all_passed']
    inp = read(paths['input'])
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    decodep = lambda ts: R.from_dict({tuple(m): s.QQ(c) for m, c in ts})
    decode = lambda q: K(decodep(q[0]))/K(decodep(q[1]))
    encodep = lambda p: [[list(m), str(c)] for m, c in p.items()]
    encode = lambda q: [encodep(q.numer), encodep(q.denom)]
    DL = [[decode(q) for q in row] for row in pc['Lframe_dual6x2']]
    pair = [[decode(q) for q in row] for row in pc['Lframe_pair2x6']]
    assert [[sum((pair[i][j]*DL[j][k] for j in range(6)), K.zero) for k in range(2)] for i in range(2)] == [[K.one, K.zero], [K.zero, K.one]]
    rows, cols = cp['remaining_rows'], cp['remaining_columns']
    group = list(range(4, 10))
    av = [decode(cp['matrix25x32'][4][j]) for j in group]
    vv = [r[1] for r in DL]
    assert vv[1] == K.zero
    row4_steps = [(step, c) for step in cp['all_completed_unit_steps'] for r, c in step['eliminations'] if r == 4]
    assert len(row4_steps) == 1 and row4_steps[0][0]['pivot_row'] == 0
    coefficient = decode(row4_steps[0][1])
    for k, j in enumerate(group):
        initial = K(decodep(inp['matrix25x32'][4][j]))
        assert av[k] == initial-coefficient*decode(cp['matrix25x32'][0][j])
        # Target x0*Sym2(y) carries the same unscaled source pair coefficients.
        assert initial == pair[1][k]
        assert decode(cp['matrix25x32'][0][j]) == pair[0][k]/K(25*a/4)
    dot = sum((av[i]*vv[i] for i in range(6)), K.zero)
    exact = {'row4_coefficients_on_x0_Sym2_y': [encode(q) for q in av],
             'transported_pair_second_column': [encode(q) for q in vv],
             'all_six_row4_source_binding_identities': True,
             'only_row4_elimination_old_source_row0_at_target0': encode(coefficient),
             'literal_zero_dual_slot': 1, 'six_term_dot': encode(dot)}
    assert dot == K.one
    universal = universal_completion()
    # Every new scalar atom uses only prior certified FULLU units. One dual
    # coefficient has an a denominator; record each denominator literally.
    known = [a, d, 500*a+90*e**2+3*e-24,
             (50*a+18*e**2-9*e-2)*d+60*e-40,
             100*a+18*e**2-9*e-2, 5*a*d-6*e-4]
    denominator_proofs = []
    for q in vv+av:
        den, powers = q.denom, []
        for f in known:
            n = 0
            while not den.rem(f):
                den = den.exquo(f); n += 1
            powers.append(n)
        assert den.is_ground and den
        denominator_proofs.append({'constant': str(den), 'a_d_F_B_C_J_powers': powers})
    exact['all_12_new_atom_denominator_proofs'] = denominator_proofs
    lazy = Lazy()
    cp_atom = lambda i, j: 0 if not cp['matrix25x32'][i][j][0] else lazy.put(['checkpoint', i, j])
    v_atom = lambda k: 0 if not vv[k] else lazy.put(['accepted_pair_second', k])
    # Put the zero dual component last. The permutation is conjugated back
    # inside the six original target coordinates, so no sign ambiguity remains.
    perm = [0, 2, 3, 4, 5, 1]
    group_order = [group[i] for i in perm]
    b = [cp_atom(4, j) for j in group_order[:5]]
    v = [v_atom(k) for k in perm[:5]]
    k = cp_atom(4, group_order[5])
    U0 = [[lazy.sub(int(i == j), lazy.mul(v[i], b[j])) for j in range(5)]+[v[i]] for i in range(5)]
    U0.append(b+[0])
    E, Ei = lazy.eye(6), lazy.eye(6)
    for j in range(5):
        E[5][j] = lazy.sub(0, lazy.mul(k, b[j]))
        Ei[5][j] = lazy.mul(k, b[j])
    U6, V6 = lazy.mm(U0, E), lazy.mm(Ei, U0)
    U24, V24 = lazy.eye(24), lazy.eye(24)
    group_positions = [cols.index(j) for j in group_order]
    for i, ri in enumerate(group_positions):
        for j, cj in enumerate(group_positions):
            U24[ri][cj], V24[ri][cj] = U6[i][j], V6[i][j]
    pivot_column = group_order[5]
    pivot_position = cols.index(pivot_column)
    assert pivot_column == 5
    outside = [j for j in cols if j not in group and cp['matrix25x32'][4][j][0]]
    assert outside == [1, 2, 3]
    tail, taili = lazy.eye(24), lazy.eye(24)
    for j in outside:
        tail[pivot_position][cols.index(j)] = lazy.sub(0, cp_atom(4, j))
        taili[pivot_position][cols.index(j)] = cp_atom(4, j)
    target, targeti = lazy.mm(U24, tail), lazy.mm(taili, V24)
    A = [[cp_atom(i, j) for j in cols] for i in rows]
    transformed = lazy.mm(A, target)
    pivot_row_position = rows.index(4)
    source = lazy.eye(17)
    for i in range(17):
        if i != pivot_row_position:
            source[i][pivot_row_position] = lazy.sub(0, transformed[i][pivot_position])
    final_rows = [i for i in rows if i != 4]
    final_cols = [j for j in cols if j != pivot_column]
    reduced = [[transformed[rows.index(i)][cols.index(j)] for j in final_cols] for i in final_rows]
    # Pivot row becomes e_p by the proved universal identities and actual dot.
    # Subtracting its column value changes only column p, which is then removed.
    exact['outside_row4_target_coefficients'] = [[j, cp['matrix25x32'][4][j]] for j in outside]
    exact['all_other_row4_coefficients_zero_outside_group_and_tail'] = True
    certificate = {'actual_six_dot_and_binding': exact, 'universal_integer_proof': universal,
                   'scalar_DAG': {'nodes': lazy.nodes,
                                  'checkpoint_reference': str(paths['checkpoint'].relative_to(ROOT)),
                                  'pair_reference': str(paths['pair_certificate'].relative_to(ROOT)),
                                  'accepted_pair_second_definition': 'Lframe_dual6x2[k][1]'},
                   'six_target_original_columns': group,
                   'permuted_six_target_original_columns': group_order,
                   'zero_dual_original_target_column': 5,
                   'U0_6x6_roots': U0, 'E_6x6_roots': E, 'E_inverse6x6_roots': Ei,
                   'U6_roots': U6, 'U6_inverse_roots': V6,
                   'target24_roots': target, 'target24_inverse_roots': targeti,
                   'tail_clear_original_columns': outside,
                   'source17_roots': source,
                   'pivot_row_transformed_before_simplification_roots': transformed[pivot_row_position],
                   'pivot_row_exact_simplified': [int(i == pivot_position) for i in range(24)],
                   'old17_source_original25_rows': rows, 'old24_target_original32_columns': cols,
                   'new16_source_original25_rows': final_rows,
                   'new16_source_labels': [inp['source25_labels'][i] for i in final_rows],
                   'new23_target_original32_columns': final_cols,
                   'remaining16x23_scalar_roots': reduced,
                   'source_pivot': {'original25_row': 4, 'original32_column': 5, 'value': '1'},
                   'all_5_remaining_heavy_and_psi_columns_unchanged': [27, 28, 29, 30, 31],
                   'source_change_determinant': '1', 'target_change_determinant': '-1',
                   'actual_new_pivots': 1, 'total_pivots_in_transformed_presentation': 9,
                   'original9minor_not_claimed': True,
                   'same_cokernel_proof': 'Actual six-dot equals1. Universal U0/E proof and three tail elementary operations turn full row4 into a unit vector. Other rows are cleared with the saved source17 matrix. Deleting that unit row and column leaves precisely the saved16x23 roots.'}
    (out/'certificate.json').write_text(json.dumps(certificate, separators=(',', ':'))+'\n')
    report = {'id': ID, 'all_six_dot_source_binding_and_universal_integer_completion_identities_passed': True,
              'literal_six_term_dot': '1', 'new_global_unit_pivots': 1,
              'total_known_global_pivots_in_projected_augmented_presentation': 9,
              'actual_new_shape_same_cokernel': [16, 23],
              'all_full_residual17_rows_and24_targets_transported_as_finite_DAG': True,
              'heavy5_and_psi_unchanged': True,
              'only_actual_row4_data_and_six_dot_expanded_other_rows_lazy': True,
              'no_new_inverses_beyond_accepted_FULLU_units': True,
              'new_core_denominators_only_nonzero_QQ_and_certified_a_unit': True,
              'no_Gaussian_minor_search_point_GB_or_automatic_next_pivot': True,
              'no_full_FULLU57_quartic239_or_uniform_kernel_count_change': True,
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths.values()},
              'certificate_sha256': sha(out/'certificate.json'),
              'scalar_DAG_nodes': len(lazy.nodes),
              'preserved_failed_metadata_attempt': {str(p.relative_to(out)): sha(p) for p in (out/'failed-attempt').glob('*')},
              'seconds': time.monotonic()-began, 'frozen': True, 'independent_review_pending': True}
    assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file()) < 32*1024*1024
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())
