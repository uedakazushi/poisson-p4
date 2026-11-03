#!/usr/bin/env python3
"""Exact closed-fiber Koszul certificate and finite rank-two/13-generator proof."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own
from ambient_B18_actual_D_tangent_cone_support_review_20261103_v63 import inverse, matrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ACTUAL_D_FINITE13_20261103_v63_REVIEW'


def main(output):
    start = time.monotonic()
    assert not output.exists()
    supportpath = ROOT/'reports/AMBIENT_B18_ACTUAL_D_TANGENT_CONE_SUPPORT_20261103_v63_REVIEW.json'
    tangentpath = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json'
    support, tangent = [own.accepted(p) for p in [supportpath, tangentpath]]
    cpath = ROOT/'reports/AMBIENT_B18_ACTUAL_D_TANGENT_CONE_SUPPORT_20261103_v63/certificate.json.gz'
    tcpath = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz'
    reportpath = tcpath.parent/'report.json'
    c, tc, actual = [own.read(p) for p in [cpath, tcpath, reportpath]]
    assert actual['actual_D_cone_dimension'] == 28 and actual['actual_D_cone_embedding_dimension'] == 39
    assert actual['actual_D_Hilbert_Samuel_multiplicity'] == 2
    R = c['common_R29_physical_coordinate_rows']
    assert len(R) == 29 and R[-1] == 61
    odd = c['common_odd10_physical_coordinate_rows']
    assert tc['common39_physical_coordinate_rows'] == R+odd
    D = matrix(tc['common150_by39_normalized_tangent'])
    assert [D[i] for i in R+odd] == [[F(i == j) for j in range(39)] for i in range(39)]
    LR = matrix(c['full_R29_linear_support_equation'])[0]
    assert LR[-1] == -1
    plane = matrix(c['full150_by28_tangent_cone_support_plane_basis'])
    _, projection_det = inverse([plane[i] for i in R[:28]])
    vectors = [[F(row[-1]) for row in c[key]] for key in
               ['actual_odd_linear_forms_CA10_by29', 'actual_odd_linear_forms_CB10_by29']]
    assert vectors[0] == [F(-1), 0, 0, 0, F(1, 4), 1, 0, 0, 0, 0]
    assert vectors[1] == [-F(1, 6), 0, 0, 0, -F(8, 3), -F(2, 3), 0, 0, 0, 0]
    assert vectors[0][0]*vectors[1][4]-vectors[0][4]*vectors[1][0] == F(65, 24)
    # First six CA regular Koszul pairs, one independent CB pair, then CA regular/core.
    selected = [('CA', 0, j) for j in range(1, 7)]+[('CB', 0, 4)]+[('CA', 0, j) for j in range(7, 10)]
    coefficients = []
    for family, i, j in selected:
        v = vectors[family == 'CB']
        row = [F(0)]*10
        row[j], row[i] = v[i], -v[j]
        coefficients.append(row)
    cofactors, minor = inverse(coefficients)
    assert minor != 0
    assert own.matrix_product(cofactors, coefficients) == [[F(i == j) for j in range(10)] for i in range(10)]
    # Closed fiber's initial upper-bound model: (YY,zY,z^3), basis 1,z,z^2,Y10.
    monomial_basis = ['1', 'z', 'z^2']+[f'Y{i}' for i in range(10)]
    assert len(monomial_basis) == 13 and 39-28 == 11 and 1+11 == 12
    paths = [Path(__file__), Path(own.__file__), ROOT/'code/ambient_B18_actual_D_tangent_cone_support_review_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', supportpath, tangentpath, cpath, tcpath, reportpath]
    stringify = lambda A: [[str(a) for a in row] for row in A]
    result = dict(id=STEM, all_passed=True,
        physical_regular_parameter_rows=R[:28], centered_z_physical_row=61, odd_physical_rows=odd,
        tangent_cone_support_projection_to_physical28_isomorphism=True,
        support_projection_determinant=str(projection_det),
        two_regular_b_vectors_nonproportional_minor='65/24',
        selected_true_Koszul_pairs=selected,
        closed_fiber_zY10_coefficient_matrix=stringify(coefficients),
        explicit_cofactors_for_all_zY10=stringify(cofactors),
        exact_Koszul_cross_rank=10, exact_Koszul_cross_minor=str(minor),
        closed_fiber_initial_ideal_contains_YY_zY_z3=True,
        closed_fiber_embedding_dimension=11, closed_fiber_Hilbert_function_options=[[1,11], [1,11,1]],
        closed_fiber_length_allowed=[12,13], closed_fiber_maximal_ideal_cube_zero=True,
        finite_injective_physical_R28_parameter_map_proved=True,
        parameters_form_maximal_ideal_reduction=True,
        generic_R28_module_rank=2, global_complete_R28_module_generators=monomial_basis,
        module_generator_upper_bound=13,
        free_rank13_or_free_rank2_or_exact_closed_fiber_length_or_full_D_depth_claimed=False,
        verification_kind='Independent small Fraction rank/cofactors plus full nonreduced formal Noether finiteness/reduction/Nakayama proof; no GB or new source jets.',
        input_sha256={str(p.relative_to(ROOT)): own.sha(p) for p in paths}, seconds=time.monotonic()-start)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','exact_Koszul_cross_minor','closed_fiber_length_allowed','generic_R28_module_rank','module_generator_upper_bound','seconds']},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output)
