#!/usr/bin/env python3
"""Independent fixed-source ell identities and finite DVR module proof audit.

No producer imports, source equation solves, source jets, or CAS calls.
The two fixed 9x9 matrices are independently inverted over Fraction; both
inverse products and all coordinate readbacks are checked.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
AUTHOR = 'AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63'
STEM = AUTHOR + '_REVIEW'


def digest(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def read(p):
    return json.loads((ROOT / p).read_text())


def mm(a, b):
    return [[sum((x*y for x, y in zip(row, col)), F(0))
             for col in zip(*b)] for row in a]


def inverse(a):
    """Own augmented inversion, followed by BOTH exact product checks."""
    n = len(a)
    identity = [[F(i == j) for j in range(n)] for i in range(n)]
    aug = [a[i][:] + identity[i][:] for i in range(n)]
    det = F(1)
    pivots = []
    for j in range(n):
        r = next(i for i in range(j, n) if aug[i][j])
        pivots.append(r)
        if r != j:
            aug[j], aug[r] = aug[r], aug[j]
            det = -det
        value = aug[j][j]
        det *= value
        aug[j] = [x / value for x in aug[j]]
        for i in range(n):
            if i != j:
                value = aug[i][j]
                aug[i] = [x - value*y for x, y in zip(aug[i], aug[j])]
    answer = [row[n:] for row in aug]
    assert mm(a, answer) == identity
    assert mm(answer, a) == identity
    return answer, det, pivots


def encode(a):
    if isinstance(a, list):
        return [encode(x) for x in a]
    return str(a)


def on_line(poly, velocity):
    result = {}
    for exp, *coefficient in poly:
        c = F(coefficient[0]) if len(coefficient) == 1 else F(*coefficient)
        degree = sum(exp)
        for e, v in zip(exp, velocity):
            c *= v**e
        result[degree] = result.get(degree, F(0)) + c
    return {e: c for e, c in result.items() if c}


def main():
    started = time.monotonic()
    out = ROOT / 'reports' / STEM
    assert not (out / 'report.json').exists()
    out.mkdir(parents=True, exist_ok=True)
    hashes = {}

    def get(p, expected=None):
        hashes[p] = digest(p)
        if expected is not None:
            assert hashes[p] == expected, p
        return read(p)

    primary = get('reports/' + AUTHOR + '/report.json',
                  '6178d9b7e8990a01a9640f91e066d713157b954de52697f639c80ae3b5435a45')
    assert primary['all_passed']
    for p, h in primary['input_sha256'].items():
        assert digest(p) == h, p
        hashes[p] = h
        if p.endswith('/report.json'):
            assert read(p)['all_passed'], p
    for name, h in primary['evidence_sha256'].items():
        p = 'reports/' + AUTHOR + '/' + name
        assert digest(p) == h
        hashes[p] = h
    claimed = read('reports/' + AUTHOR + '/actual-ell-line-values-and-nil-DVR-proof.json')

    z = read('reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63/finite-cokernel-Z-and-actions.json')
    interface = read('reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json')
    assert z['physical_R0_rows9'] == [27,108,109,112,113,117,118,119,61]
    assert interface['retained_source9_order'] == ['a','d','v2','v3','v4','v5','v6','v7','v8']
    w = [[F(c) for c in row] for row in z['actual_physical_R0_firstjets_on_CA_F9_9by9']]
    k = [[F(c) for c in row] for row in interface['actual_source9_by_physical9_Jacobian']]
    r = [[F(c) for c in row] for row in interface['physical_R0_first9_by9']]
    wi, wd, wp = inverse(w)
    ki, kd, kp = inverse(k)
    assert kd == F(5,48)
    ell_f9 = [F(0)]*7 + [F(-5), F(5)]
    ell_physical = mm([ell_f9], wi)[0]
    ell_c = mm([ell_physical], r)[0]
    ell_source = mm([ell_c], ki)[0]
    assert mm([ell_physical], w) == [ell_f9]
    assert mm([ell_source], k) == [ell_c]
    assert ell_physical == [0,0,0,F(-15,2),0,F(10,3),0,0,0]
    assert ell_source == [F(-5,6),F(5,3),F(15,2),0,0,0,F(-10,3),0,0]
    for field, actual in [('ell_F9_first',ell_f9), ('ell_physical_first',ell_physical),
                          ('ell_regular_C_first',ell_c), ('ell_proper_source9_first',ell_source)]:
        assert [F(x) for x in claimed[field]] == actual

    # Bind the exact support lines to their accepted source-polynomial primes.
    five_report = read('reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json')
    assert five_report['E4_Ass_six_and_embedded_N_P_accepted']
    for p,h in five_report['input_sha256'].items():
        if p.endswith('/five-line-and-generic-models.json'):
            five = get(p,h)
    assert five['normalized_E4']['Ass'] == ['(d)','N','C','R','P','L']
    assert five['actual_core']['geometric_generic_lengths'] == [3,4,2,2,1]
    lines = {'N':[0,0,0,1], 'C':[0,1,1,0], 'R':[-2,1,-1,-2],
             'P':[2,0,-2,1], 'L':[1,1,-1,F(-1,2)]}
    ell4 = [ell_source[i] for i in (0,1,2,6)]
    values = {}
    for name, velocity in lines.items():
        for prime_generator in five['support']['primes'][name]:
            assert on_line(prime_generator, velocity) == {}
        for p in ('source_Q0','source_Q1','source_f'):
            assert on_line(five[p], velocity) == {}
        values[name] = sum((a*b for a,b in zip(ell4,velocity)), F(0))
        assert values[name] != 0
    assert values == {'N':F(-10,3),'C':F(55,6),'R':F(5,2),'P':F(-20),'L':F(-5)}
    assert {n:F(v) for n,v in claimed['all5_proper_core_line_first_values'].items()} == values
    assert any(ell_source[i] for i in (0,2,6))  # ell not in main prime (d).

    # The accepted rational source pivots differ from u5 by order >=2.
    core_peer = get('reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW/report.json')
    assert core_peer['all_passed'] and core_peer['five_roots_all7_residuals_source_verified']
    root_path = 'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json'
    assert root_path in core_peer['input_sha256']
    roots = get(root_path, core_peer['input_sha256'][root_path])
    assert roots['actual_core4_variable_order'] == ['a','d','v2','v6']
    pivot_orders = {}
    for name, data in roots['five_actual_rational_roots'].items():
        order = min(sum(e) for e,n,d in data['numerator'] if n)
        assert order >= 2 and data['D_power'] >= 0 and data['P_power'] >= 0
        pivot_orders[name] = order
    assert set(pivot_orders) == {'v3','v4','v5','v7','v8'}
    assert all(ell_source[i] == 0 for i in (3,4,5,7,8))

    sq = read('reports/AMBIENT_B18_NORMALIZER_NIL_SQUARE_THREE_SUPPORTS_20261104_v63_REVIEW/report.json')
    assert sq['Nil_E4_square_Ass_exact'] == ['N','P','C']
    assert sq['generic_lengths'] == [1,1,1] and sq['dimension_depth_multiplicity'] == [1,1,3]
    rank_k = read('reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json')
    rank_c = read('reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json')
    rank_cp = read('reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json')
    hf = read('reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json')
    assert rank_k['rank_D_K'] == 4
    assert rank_c['actual_comparison_C_generic_DVR_rank'] == 4
    assert rank_cp['generic_DVR_rank_Cprime'] == 2
    interval = hf['actual_T_length_interval_accepted']
    assert interval == [11,17]
    generic_core = dict(zip(lines, five['actual_core']['geometric_generic_lengths']))
    nil_lengths = {name: length-(name not in ('N','P')) for name,length in generic_core.items()}
    assert nil_lengths == {'N':3,'C':3,'R':1,'P':2,'L':0}
    free_e = sum(nil_lengths.values())
    free_image = free_e - 4 + 2
    assert (free_e,free_image,free_image+4) == (9,7,11)

    proof = {
        'method':'Own Fraction augmented inverses, both product checks, source-prime substitutions; accepted full source/Ass gates, no producer import.',
        'matrix_W_determinant':str(wd), 'matrix_W_inverse':encode(wi), 'matrix_W_pivot_rows':wp,
        'matrix_K_determinant':str(kd), 'matrix_K_inverse':encode(ki), 'matrix_K_pivot_rows':kp,
        'four_9by9_inverse_products_checked':True,
        'ell_F9':encode(ell_f9), 'ell_same_physical_R0':encode(ell_physical),
        'ell_regular_C9':encode(ell_c), 'ell_same_proper_source9':encode(ell_source),
        'source_pivot_root_numerator_orders':pivot_orders,
        'all5_exact_prime_and_Q0_Q1_f_line_substitutions_zero':True,
        'five_actual_line_ell_first_values':{n:str(v) for n,v in values.items()},
        'allorders_regularity_proof':'Each nonzero line coefficient makes ell a parameter by formal inversion, regardless of higher terms. Non-d first terms exclude main(d). All SIX accepted Ass, including embedded N/P, are missed; thus ell is regular on E4 and Ered.',
        'finite_module_argument':'NilE4 is finite over the complete source ring and supported on these curves. Its ell quotient has finite length; complete Nakayama/formal coefficient lifting proves finiteness over k[[ell]]. E4-injectivity of ell makes it torsionfree, hence free. Each reduced line has field degree one over k((ell)).',
        'generic_nil_lengths':nil_lengths, 'normalized_nil_D_free_rank':9,
        'embedded_lengths_explanation':'At N/P NilE4=d*(core) as a MODULE, so its lengths are the full core lengths3/2, not core maximal-ideal lengths. At C/R/L d is a unit, subtract one. Main(d) contributes0.',
        'normalized_nil_square_D_free_rank':sum(sq['generic_lengths']),
        'nil_pullback_exactness':'K is nil. If an element maps to nil, a power lies in K and a further power is zero. Therefore preimage(NilE4)=NilD0 and 0->K->NilD0->NilE4->C->Cprime->0 is exact.',
        'image_free_rank':free_image,
        'split_proof':'The image in freeD9 is finite torsionfree of rank9-(4-2)=7. It is free over the DVR, so 0->K->NilD0->D7->0 splits NONCANONICALLY as D-modules.',
        'NilD0_D_free_plus_K':'K direct_sum D^7, D-module ONLY',
        'rank_D_NilD0':11, 'mu_D_NilD0_formula':'L+3',
        'mu_D_NilD0_interval':[x+3 for x in interval],
        'torsion_submodule_equality':'tors_D NilD0=tors_D K=H0(D0) under the actual inclusion',
        'torsion_mu_formula':'L-8', 'torsion_mu_interval':[x-8 for x in interval],
        'future_Tnil_length_L_plus5':'Conditional on constructing the actual honest source comparison; NOT constructed here.',
        'not_claimed':['R0/D0/algebra split','explicit seven nil lifts','IAIB annihilates all nil','entire full graded nil is D-finite','torsion length/exponents','cube answer or exact nil index','104x948 presentation']}
    ep = out / 'independent-source-ell-and-module-proof.json'
    ep.write_text(json.dumps(proof,indent=2)+'\n')
    for p in [str(Path(__file__).relative_to(ROOT)), 'research/notes/'+STEM+'.md']:
        hashes[p] = digest(p)
    report = {
        'id':STEM, 'all_passed':True,
        'same_physical_ell_two_fixed9_source_maps_independently_checked':True,
        'whole_E4_and_Ered_ell_regularity_all6_Ass_accepted':True,
        'normalized_nil_DVR_finite_free_rank':9,
        'normalized_nil_square_DVR_finite_free_rank':3,
        'actual_NilD0_DVR_finite_rank':11,
        'actual_NilD0_noncanonical_DVR_split':'K plus D^7, D-module ONLY',
        'actual_NilD0_mu_interval':[14,20], 'actual_NilD0_torsion_mu_interval':[3,9],
        'no_new_source_solve_jet_CAS_point_or_producer_rerun':True,
        'new_nil_lifts_Tnil_presentation_cube_or_torsion_length_claimed':False,
        'input_sha256':hashes,
        'evidence_sha256':{ep.name:hashlib.sha256(ep.read_bytes()).hexdigest()},
        'seconds':time.monotonic()-started,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}))


if __name__ == '__main__':
    main()
