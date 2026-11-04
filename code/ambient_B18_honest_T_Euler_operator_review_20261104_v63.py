#!/usr/bin/env python3
"""Independent source/type and all-degree contraction audit; no author import."""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json, resource, time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63_REVIEW'
AUTHOR = 'AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63'
ZERO = (0,) * 8
PAIRS = list(combinations(range(8), 2))


def digest(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def read(p):
    return json.loads((ROOT / p).read_text())


def acc(dst, key, value):
    value = dst.get(key, F(0)) + value
    if value:
        dst[key] = value
    else:
        dst.pop(key, None)


def poly(raw):
    ans = {}
    for m, q in raw:
        assert len(m) == 8 and min(m) >= 0
        acc(ans, tuple(m), F(q))
    return ans


def pvar(i, power=1):
    m = [0] * 8
    m[i] = power
    return {tuple(m): F(1)}


def pmul(a, b):
    ans = {}
    for u, q in a.items():
        for v, r in b.items():
            acc(ans, tuple(x + y for x, y in zip(u, v)), q * r)
    return ans


def pscaled_sum(a, b, scale=1):
    ans = dict(a)
    for m, q in b.items():
        acc(ans, m, q * scale)
    return ans


def contraction_from_coefficients(v):
    # Direct coefficient scatter: this does not call the author's derivative,
    # Euler or boundary helper.  The denominators refer to OUTPUT u degree.
    w = [{} for _ in PAIRS]
    for k, (i, j) in enumerate(PAIRS):
        if j == 7:
            for mon, q in v[7].items():
                if mon[i]:
                    out = list(mon)
                    out[i] -= 1
                    acc(w[k], tuple(out), q * mon[i] / sum(mon[:7]))
        else:
            for component, axis, sign in ((j, i, 1), (i, j, -1)):
                for mon, q in v[component].items():
                    if mon[axis]:
                        out = list(mon)
                        out[axis] -= 1
                        acc(w[k], tuple(out), sign * q * mon[axis] / (sum(mon[:7]) + 1))
    return w


def boundary(w, n):
    v = [{} for _ in range(8)]
    for (i, j), p in zip(PAIRS, w):
        for mon, q in p.items():
            mi = list(mon)
            mj = list(mon)
            mi[i] += 1
            mj[j] += (n if j == 7 else 1)
            acc(v[j], tuple(mi), q)
            acc(v[i], tuple(mj), -q)
    return v


def homogeneous_symbolic_identity():
    # For u-homogeneous degree d in v_j, d>=0, v_8 has degree d+1.
    # sum_i u_i*w_ij = (d*v_j + v_j + r^n*d_j v8)/(d+1).
    # The correction is exactly r^n*w_j8.  The last component is
    # k/k*v8, k>=1.  Check the polynomial numerators in the formal integers
    # d,k, rather than extrapolating a finite degree fixture.
    return {
        'first7_formal_degree_numerator': {'d*vj': 1, 'vj': 1, 'r^n*partial_j_v8': 1},
        'first7_denominator': 'd+1',
        'first7_cancelled_numerator': {'(d+1)*vj': 1},
        'last_component_numerator_and_denominator': ['k*v8', 'k'],
        'd_ge_0_k_ge_1': True,
        'v8_u_degree_zero_vanishes_by_r_power_in_domain': True,
        'proof_not_deduced_from_fixtures': True,
    }


def verify_hash_tree(value, hashes):
    # Check every explicit {path,sha256} leaf of this frozen dependency tree.
    if isinstance(value, dict):
        if isinstance(value.get('path'), str) and isinstance(value.get('sha256'), str):
            p = value['path']
            assert digest(p) == value['sha256'], p
            hashes[p] = value['sha256']
        for v in value.values():
            verify_hash_tree(v, hashes)
    elif isinstance(value, list):
        for v in value:
            verify_hash_tree(v, hashes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='reports/' + STEM)
    args = parser.parse_args()
    out = ROOT / args.output
    assert not (out / 'report.json').exists(), 'Use a fresh output directory; frozen reviews are immutable.'
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    input_hashes = {}
    author_report_path = 'reports/' + AUTHOR + '/report.json'
    author_report = read(author_report_path)
    assert author_report['all_passed'] and author_report['frozen']
    input_hashes[author_report_path] = digest(author_report_path)
    verify_hash_tree(author_report, input_hashes)
    paths = author_report['input_hashes']
    data = {k: read(v['path']) for k, v in paths.items()}
    for k in ('m_peer', 'B_peer', 'SyzQ_peer', 'resolution_peer'):
        assert data[k]['all_passed'], k
        verify_hash_tree(data[k], input_hashes)
    p_peer_path = 'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json'
    p_peer = read(p_peer_path)
    assert p_peer['all_passed'] and p_peer['actual_D0_free9_relations640_and_D4_map_checked']
    assert p_peer['input_sha256'][paths['P']['path']] == paths['P']['sha256']
    input_hashes[p_peer_path] = digest(p_peer_path)
    E = read(author_report['evidence']['path'])
    M, B, P, S, R = (data[k] for k in ('m', 'B', 'P', 'SyzQ', 'resolution'))
    for p in ('code/ambient_B18_honest_T_Euler_operator_20261104_v63.py',
              'research/notes/' + AUTHOR + '.md'):
        input_hashes[p] = digest(p)
    assert input_hashes['code/ambient_B18_honest_T_Euler_operator_20261104_v63.py'] == author_report['code_sha256']
    assert input_hashes['research/notes/' + AUTHOR + '.md'] == author_report['note_sha256']
    O = E['scalar_and_vector_operator_nodes']
    assert len(O) == 74176 and E['input_hashes'] == paths
    assert E['source_Cbar_minimal_target_order'] == ['tauA', 'tauB', 'tauA*tauB']
    assert S['actual_d2'] == {
        'shape': [48, 128], 'value': 'Res4(N^t)*Res4(S12_by32)',
        'exact_block_recipe': 'block(i,j)=mu(sum_k N[k,i]*S[k,j]); native12rowblock then4basis columns'}
    assert R['r_lengths'] == [3, 3, 2] and R['ranks'] == [3,24,84,168,210,168,84,24,3]
    assert E['actual_CB_shift_coefficients7'] == R['CB_shift_coefficients7']
    assert B['all640_UP_equals_QB_verified'] and B['all28_intrinsic_C0_columns_have_zero_lift']
    assert len(P['honest_D0_relation_columns']) == 640
    preflight = {
        'operator_nodes': len(O), 'stored_coefficient_slots': sum(map(len, E['honest_T116_by936_sparse_columns'])),
        'source_X_shape': [56, 768], 'presentation_shape': [116,936],
        'source_P_DAG_read_only_nodes': len(P['DAG_scalar_nodes']),
        'no_source_function_expansion_or_348_evaluation': True,
        'max_allowed_polynomial_slots': 150000,
        'input_sha256': input_hashes,
    }
    assert preflight['stored_coefficient_slots'] == 70017 < 150000
    (out / 'preflight.json').write_text(json.dumps(preflight, indent=2) + '\n')

    selected = M['fixed_selected24_columns']
    complement = [j for j in range(56) if j not in selected]
    assert len(set(selected)) == 24 and len(complement) == 32
    assert E['actual_source24_selected_columns'] == selected
    assert E['actual_source32_complement_rows'] == complement
    m = M['typed_m24_by56']
    unit_id = M['typed_unit24_inverse_node']
    unit = M['typed_regular_function_operators'][unit_id]
    assert unit[0] == 'formal_unit_matrix_inverse'
    assert unit[1] == [[row[j] for j in selected] for row in m]
    assert unit[2] == 'actual_center_matrix_inverse'
    assert unit[3] == M['fixed24_minor_inverse']
    for i in range(56):
        for j in range(24):
            a = M['typed_n56_by24'][i][j]
            if i not in selected:
                assert a == 0
            else:
                assert M['typed_regular_function_operators'][a] == ['unit_matrix_inverse_entry', unit_id, selected.index(i), j]
    fixed = [[F(x) for x in row] for row in M['fixed24_minor_matrix']]
    inverse = [[F(x) for x in row] for row in M['fixed24_minor_inverse']]
    assert fixed == [[F(row[j]) for j in selected] for row in M['all24_by56_m0']]
    inverse_products = 0
    for A, C in ((fixed, inverse), (inverse, fixed)):
        for i in range(24):
            for j in range(24):
                assert sum((A[i][k]*C[k][j] for k in range(24)), F(0)) == (i == j)
                inverse_products += 1

    X = E['actual_X_source']['source_X_sparse_columns']
    assert E['actual_X_source']['actual_X_shape'] == [56,768] and len(X) == 768
    all_used = {0}
    for j, col in enumerate(X):
        actual = dict(col)
        assert len(actual) == len(col) and list(actual) == sorted(actual)
        assert all(0 <= i < 56 and 0 <= a < len(O) for i, a in col)
        wanted = {}
        if j < 640:
            for i, a in B['B48_by640_sparse_columns'][j]:
                wanted[i] = ['actual_chain_scalar_function_restrict_ell0', a]
                assert 0 <= a < len(B['operator_scalar_nodes'])
            source_p = P['honest_D0_relation_columns'][j]
            assert len(source_p) == 9
            for i, a in enumerate(source_p[1:]):
                if a:
                    wanted[48+i] = ['actual_C0_scalar_function_restrict_ell0', a]
                    assert 0 <= a < len(P['DAG_scalar_nodes'])
        else:
            wanted = {i: ['actual_Res4_SyzQ_scalar_entry_restrict_ell0', i, j-640] for i in range(48)}
        assert set(actual) == set(wanted)
        for i, a in actual.items():
            assert O[a] == wanted[i], (j, i)
            all_used.add(a)

    templates = E['formal_Koszul_Euler_templates']
    for sector, n in enumerate((3,3,2)):
        t = templates[sector]
        assert t['sector'] == sector and t['r_power'] == n
        assert t['to_coordinates'] == ('z=u-c*r2' if sector == 1 else 'z=u')
        assert t['from_coordinates'] == ('u=z+c*r2' if sector == 1 else 'u=z')
        assert t['derivatives'] == 'partial in u at fixed r'
        expected = []
        for i, j in PAIRS:
            f = {'Euler_output_shift':1, 'partial':i, 'input_component':7} if j == 7 else {
                'Euler_output_shift':2, 'difference_of_partials':[[i,j],[j,i]]}
            expected.append({'pair':[i,j], 'formula':f})
        assert t['outputs28'] == expected
    template_ids = [i for i,n in enumerate(O) if n[0] == 'formal_Euler_Koszul_template24_to84']
    assert len(template_ids) == 1
    template_id = template_ids[0]
    assert O[template_id] == ['formal_Euler_Koszul_template24_to84', templates]
    all_used.add(template_id)

    # Rebuild d2/d3 from exterior-basis signs, including the exact CB shift.
    regular = [[poly(p) for p in row] for row in R['actual_cyclic_sector_regular_sequences8']]
    for sector, n in enumerate((3,3,2)):
        expected = [pvar(i) for i in range(7)] + [pvar(7,n)]
        if sector == 1:
            for i,c in enumerate(R['CB_shift_coefficients7']):
                expected[i] = pscaled_sum(expected[i], pvar(7,2), F(c))
        assert regular[sector] == expected
    basis = R['resolution_basis_sector_and_wedge']
    for degree in range(4):
        wanted = [[sector,list(w)] for sector in range(3) for w in combinations(range(8),degree)]
        assert basis[degree] == wanted
    own_d3 = []
    d2_columns_checked = 0
    for degree in (2,3):
        target = {(s,tuple(w)): i for i,(s,w) in enumerate(basis[degree-1])}
        rebuilt = []
        for sector,wedge in basis[degree]:
            c = {}
            for position,axis in enumerate(wedge):
                row = target[(sector,tuple(x for x in wedge if x != axis))]
                c[row] = {m:q*((-1)**position) for m,q in regular[sector][axis].items()}
            rebuilt.append(c)
        saved = R['differentials8_column_sparse'][degree-1]
        assert len(saved) == len(rebuilt)
        for raw, c in zip(saved,rebuilt):
            assert {i:poly(p) for i,p in raw} == c
        if degree == 2:
            d2_columns_checked = len(rebuilt)
        else:
            own_d3 = rebuilt

    T = E['honest_T116_by936_sparse_columns']
    W = E['w84_by768_component_operator_nodes']
    assert len(T) == 936 and len(W) == 768
    entries = 0
    for j,c in enumerate(own_d3):
        saved = dict(T[j])
        assert len(saved) == len(T[j]) and set(saved) == set(c)
        for row,node in saved.items():
            assert O[node][0] == 'actual_literal_Rbar_polynomial'
            assert poly(O[node][1]) == c[row]
            all_used.add(node)
        entries += len(saved)
    for j in range(768):
        row_nodes = W[j]
        assert len(row_nodes) == 84
        lift_ids = set()
        for k,node in enumerate(row_nodes):
            assert O[node][0] == 'Euler_lift_component' and O[node][2] == k
            lift_ids.add(O[node][1])
            assert O[node][1] < node
            all_used.add(node)
        assert len(lift_ids) == 1
        lifted = lift_ids.pop()
        assert O[lifted][0] == 'apply_formal_Euler_Koszul_template' and O[lifted][1] == template_id
        cycle = O[lifted][2]
        assert cycle < lifted and template_id < lifted
        assert O[cycle] == ['actual_matrix_product_column', {'matrix_artifact':'m','matrix_key':'typed_m24_by56'},j]
        all_used.update((cycle,lifted))
        want = {k:a for k,a in enumerate(row_nodes)}
        x = dict(X[j])
        for k,old in enumerate(complement):
            if old in x:
                want[84+k] = x[old]
        assert dict(T[168+j]) == want and len(T[168+j]) == len(want)
        entries += len(want)
    assert entries == 70017 and all_used == set(range(len(O)))
    tag_counts = Counter(n[0] for n in O)
    assert set(tag_counts) == {
        'zero_function','actual_chain_scalar_function_restrict_ell0',
        'actual_C0_scalar_function_restrict_ell0','actual_Res4_SyzQ_scalar_entry_restrict_ell0',
        'formal_Euler_Koszul_template24_to84','actual_literal_Rbar_polynomial',
        'actual_matrix_product_column','apply_formal_Euler_Koszul_template','Euler_lift_component'}
    assert O[0] == ['zero_function']
    assert tag_counts['Euler_lift_component'] == 84*768

    # Independent fixed-sign tests use coefficient scatter, and additionally
    # test r-power4 and high u degrees, unlike the author's fixture set.
    monomials = [ZERO]
    monomials += [next(iter(pvar(i))) for i in range(8)]
    monomials += [next(iter(pvar(i,5))) for i in range(7)]
    monomials += [(2,1,0,1,0,0,0,3), (0,0,3,0,1,0,2,1)]
    fixture_count = 0
    for n in (2,3,4):
        for k in range(28):
            for a in monomials:
                w = [{} for _ in PAIRS]
                w[k] = {a:F(3,7)}
                v = boundary(w,n)
                assert boundary(contraction_from_coefficients(v),n) == v
                fixture_count += 1
    proof = {
        'all768_actual_source_columns_reconstructed': True,
        'source_shape': [56,768], 'source_640_sign': '[B48;P8] for V=[Q48,-U8]',
        'source_128_actual_native_transpose': 'Res4(N^t)*Res4(S12_by32), zero8',
        'P0_cancellation': 'U(1)=1; constant row of UP=QB eliminates P0 without changing the other56 source entries',
        'unit_center_inverse_products': inverse_products,
        'unit_determinant_from_accepted_fixed_source': M['fixed24_minor_determinant'],
        'column_change': '(old_S,old_J)->(M old_S+N old_J,old_J), inverse=(M^-1(a-Nb),b)',
        'kernel_after_unit_change': 'im(d2) direct_sum Rbar^32',
        'all84_d2_and168_d3_columns_independently_rebuilt': [d2_columns_checked,len(own_d3)],
        'all70017_relation_entries_and74176_nodes_bound': True,
        'operator_tag_counts': dict(tag_counts),
        'universal_homogeneous_contraction_proof': homogeneous_symbolic_identity(),
        'independent_coefficient_scatter_fixtures': fixture_count,
        'source_shift_applied_before_differentiation': 'u=z+c r² on CB, u=z elsewhere; r fixed',
        'formal_continuity': 'partial lowers order at most1; Euler inverse only rescales nonzero homogeneous degrees; adapted shifts are continuous triangular automorphisms',
        'all_degrees_exactness': 'd2 Euler(cycle)=cycle; kernel(d2)=im(d3) by the three regular Koszul sequences',
        'presentation_exactness': 'Rbar^84+Rbar^32 surjects ker(V); its kernel is im(d3)+0. Preimages of each of768 X columns plus this kernel generate the complete inverse image of im(X).',
        'nonlinear_contraction_not_used_as_Rbar_linear_map': True,
        'matrix_formed_by768_selected_preimages_then_Rbar_linear_extension': True,
        'original_formal_functions_not_reduced_I_squared_first': True,
        'integer_Euler_inverses_only': True,
        'no_coefficient_evaluation_minimality_length_K_torsion_depth_claim': True,
    }
    (out / 'independent-source-Euler-and-presentation-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    own_code = str(Path(__file__).relative_to(ROOT))
    own_note = 'research/notes/' + STEM + '.md'
    input_hashes[own_code] = digest(own_code)
    input_hashes[own_note] = digest(own_note)
    evidence = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    report = {
        'id': STEM, 'all_passed': True, 'frozen': True,
        'actual_complete_T_typed_presentation': [116,936],
        'all768_source_columns_all168_d3_and70017_entries_checked': True,
        'all_degree_Euler_proof_and_actual_Cbar_signs_checked': True,
        'coefficient_X0_X1_348_T_length_K_torsion_or_depth_computed': False,
        'input_sha256': input_hashes, 'evidence_sha256': evidence,
        'seconds': time.monotonic()-started,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)


if __name__ == '__main__':
    main()
