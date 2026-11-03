#!/usr/bin/env python3
"""Independent exact common-base odd conversion and tangent-cone support review."""
import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ACTUAL_D_TANGENT_CONE_SUPPORT_20261103_v63_REVIEW'


def matrix(A):
    return [[F(q) for q in row] for row in A]


def inverse(A):
    n = len(A)
    B = [list(row)+[F(i == j) for j in range(n)] for i, row in enumerate(A)]
    determinant = F(1)
    for i in range(n):
        p = next(j for j in range(i, n) if B[j][i])
        if p != i:
            B[i], B[p] = B[p], B[i]
            determinant = -determinant
        pivot = B[i][i]
        determinant *= pivot
        B[i] = [a/pivot for a in B[i]]
        for j in range(n):
            if j != i and B[j][i]:
                a = B[j][i]
                B[j] = [b-a*c for b, c in zip(B[j], B[i])]
    return [row[n:] for row in B], determinant


def main(output):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    assert not output.exists()
    authorpath = ROOT/'reports/AMBIENT_B18_ACTUAL_D_TANGENT_CONE_SUPPORT_20261103_v63/report.json'
    reviewpaths = [ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json']
    author, *reviews = [own.accepted(p) for p in [authorpath, *reviewpaths]]
    cpath = authorpath.parent/'certificate.json.gz'
    redpath = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    firstpath = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    tcpath = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz'
    c, red, first, tc = [own.read(p) for p in [cpath, redpath, firstpath, tcpath]]
    R = red['common28_linear_target_rows']+[red['common29th_target_row']]
    odd = tc['common10_odd_target_rows']
    assert R == c['common_R29_physical_coordinate_rows'] and odd == c['common_odd10_physical_coordinate_rows']
    even, blinear, nu = [], [], []
    # A fresh augmented inverse is different from the producer's direct use of the
    # saved 'actual' b-linear arrays and reconstructs the old two-coordinate conversion.
    for name in ['CA', 'CB']:
        data = first['CA_CB_full_firstjet_data'][name]
        J = matrix(red['families'][name]['linear_J'])
        W, Z = matrix(data['W']), matrix(data['rawZ'])
        augmented = [row+w for row, w in zip(J, W)]
        inv, det = inverse([augmented[i] for i in R])
        E = own.matrix_product(augmented, inv)
        raw = own.matrix_product([row+z for row, z in zip(matrix(data['mixed']), Z)], inv)
        correction = own.matrix_product(E, [raw[i] for i in R])
        b = [[a-z for a, z in zip(row, corr)] for row, corr in zip(raw, correction)]
        assert E == matrix(tc['families'][name]['full150_by29_even_cotangent_embedding'])
        assert b == matrix(data['actual']) and inv[-1] == matrix(data['discriminant'])[0]
        assert E == matrix(c['common_even_embedding150_by29'])
        even.append(E); blinear.append([b[i] for i in odd]); nu.append(inv[-1])
    assert even[0] == even[1]
    assert nu == matrix(first['CA_CB_discriminant_differentials']) == matrix(c['actual_discriminant_linear_forms2_by29'])
    assert blinear[0] == matrix(c['actual_odd_linear_forms_CA10_by29'])
    assert blinear[1] == matrix(c['actual_odd_linear_forms_CB10_by29'])
    assert [even[0][i] for i in R] == [[F(i == j) for j in range(29)] for i in range(29)]
    O = [row[:19] for row in matrix(red['families']['CA']['linear_J'])]
    OR = [O[i] for i in R]
    assert OR == matrix(c['orbit19_in_R29']) and own.matrix_product(even[0], OR) == O
    complement = c['orbit_complement_R29_indices']
    assert complement == [6, 13, 14, 16, 17, 18, 19, 20, 27, 28]
    Q = [[F(i == j) for j in complement] for i in range(29)]
    full = [a+b for a, b in zip(OR, Q)]
    fullinv, orbitdet = inverse(full)
    projection = fullinv[19:]
    assert projection == matrix(c['orbit_quotient_projection10_by29'])
    assert own.matrix_product(projection, OR) == [[0]*19 for _ in range(10)]
    assert own.matrix_product(projection, Q) == [[F(i == j) for j in range(10)] for i in range(10)]
    assert own.matrix_product(nu, OR) == [[0]*19 for _ in range(2)]
    assert all(own.matrix_product(b, OR) == [[0]*19 for _ in range(10)] for b in blinear)
    assert [R[i] for i in complement] == [27, 108, 109, 112, 113, 117, 118, 119, 140, 61]
    variables = own.symbols([f'r{i}' for i in range(10)])
    values = list(variables.values())
    linear = lambda row: sum((a*x for a, x in zip(row, values)), own.Poly(0))
    ns = [linear(row) for row in own.matrix_product(nu, Q)]
    bs = [[linear(row) for row in own.matrix_product(b, Q)] for b in blinear]
    assert all(a == own.evaluate(q, variables) for a, q in zip(ns, c['discriminant_linear_forms']))
    for a, key in zip(bs, ['odd_linear_forms_CA', 'odd_linear_forms_CB']):
        assert all(x == own.evaluate(q, variables) for x, q in zip(a, c[key]))
    L = 2*variables['r0']-variables['r9']
    active = [0, 1, 2, 4, 6, 7, 9]
    qm = [values[i]*values[j] for i, j in itertools.combinations_with_replacement(active, 2)]
    exponents = [next(iter(m.terms)) for m in qm]
    pairs = list(itertools.combinations_with_replacement(range(10), 2))
    K = [ns[1]*bs[1][i]*bs[1][j]-ns[0]*bs[0][i]*bs[0][j] for i, j in pairs]
    assert all(a == own.evaluate(q, variables) for a, q in zip(K, c['all55_pure_even_cubic_projections']))
    quotients, nonzero = [], []
    for pair, cubic in zip(pairs, K):
        if not cubic.terms:
            continue
        remainder, quotient = cubic, own.Poly(0)
        while remainder.terms:
            e = max(remainder.terms)
            assert e[0] >= 1
            f = (e[0]-1,)+e[1:]
            term = own.Poly({f: remainder.terms[e]/2})
            quotient += term
            remainder -= term*L
        assert all(e in exponents for e in quotient.terms)
        quotients.append(quotient); nonzero.append(list(pair))
    assert len(quotients) == 28 and nonzero == c['nonzero_product_pairs']
    assert all(a == own.evaluate(q, variables) for a, q in zip(quotients, c['all28_quotient_quadrics']))
    M = [[q.terms.get(e, F(0)) for e in exponents] for q in quotients]
    assert M == matrix(c['full28_by28_quotient_coefficient_matrix'])
    Mt = list(map(list, zip(*M)))
    invMt, determinant = inverse(Mt)
    assert determinant == F(c['quotient_coefficient_determinant'])
    target = [(L*L).terms.get(e, F(0)) for e in exponents]
    combo = own.matvec(invMt, target)
    assert combo == list(map(F, c['exact_linear_combination_for_L_cube']))
    assert sum((a*q*L for a, q in zip(combo, quotients)), own.Poly(0)) == L**3
    LR = [2*a-b for a, b in zip(projection[0], projection[9])]
    assert [LR] == matrix(c['full_R29_linear_support_equation'])
    # Independent actual N pinch projection in the same physical R29 order.
    JN = [[q for j, q in enumerate(row) if j != 20] for row in matrix(red['families']['CA']['linear_J'])]
    WN = matrix(first['N_full_firstjet_data']['W'])
    augN = [row+w for row, w in zip(JN, WN)]
    inverseN, detN = inverse([augN[i] for i in R[:28]])
    EN = own.matrix_product(augN, inverseN)
    Nlinear = [-q for q in EN[61]]+[F(1)]
    assert Nlinear == [-q for q in LR] == matrix(c['actual_N_pinch_projection_linear_form'])[0]
    plane = matrix(c['full150_by28_tangent_cone_support_plane_basis'])
    planeR = [plane[i] for i in R]
    assert own.matrix_product(even[0], planeR) == plane
    assert own.matvec(list(zip(*planeR)), LR) == [0]*28
    # A nonzero LR and 28 independent columns form the complete support hyperplane.
    pivot = next(i for i, q in enumerate(LR) if q)
    _, planedet = inverse([row for i, row in enumerate(planeR) if i != pivot])
    paths = [Path(__file__), Path(own.__file__), ROOT/f'research/notes/{STEM}.md', authorpath,
             cpath, redpath, firstpath, tcpath, *reviewpaths]
    result = dict(id=STEM, all_passed=True,
        independent_Fraction_augmented29_inverse_reconstructs_full_CA_CB_even_and_odd_embeddings=True,
        actual_discriminants_trace_compatible_blinear_and_common_even_E_rebuilt=True,
        independent_orbit19_quotient_and_all_forms_annihilate_orbit=True,
        independent_all55_cubic_projections_and_all28_exact_linear_divisions=True,
        nonzero_cubic_projection_count=28, exact_quotient_rank=28,
        exact_quotient_determinant=str(determinant), exact_Lcube_combo_reconstructed=True,
        actual_N_pinch_and_full150_support_plane_binding_checked=True,
        true_product_cancellation_lifts_include_unknown_Yquadratic_terms=True,
        radical_of_full_true_initial_ideal='(Y1,...,Y10,L)',
        tangent_cone_reduced_support_dimension=28, orbit_transverse_reduced_support_dimension=9,
        tangent_cone_nilradical_power6_zero=True,
        actual_D_nilradical_index_full_initial_ideal_or_exact_depth_claimed=False,
        no_Groebner_or_new_source_jets=True,
        verification_kind='Independent pure Fraction inverses/polynomials/determinant and conceptual tangent-cone proof; producer not rerun.',
        input_sha256={str(p.relative_to(ROOT)): own.sha(p) for p in paths}, seconds=time.monotonic()-start)
    output.write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)
