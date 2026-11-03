#!/usr/bin/env python3
"""Exact pure-even cubic projections; reduced actual tangent cone and nil bound."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ACTUAL_D_TANGENT_CONE_SUPPORT_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def mat(A):
    return [[str(q) for q in row] for row in A.tolist()]


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    reduction_path = ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    first_path = ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    tangent_path = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/report.json'
    tangent_cert = tangent_path.parent/'certificate.json.gz'
    reviews = [
        ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json',
    ]
    assert all(read(p)['all_passed'] for p in reviews)
    reduction, first, tangent, tc = map(read, [reduction_path, first_path, tangent_path, tangent_cert])
    assert tangent['actual_D_cone_dimension'] == 28 and tangent['actual_D_cone_embedding_dimension'] == 39
    E = [s.Matrix(tc['families'][name]['full150_by29_even_cotangent_embedding']) for name in ['CA', 'CB']]
    assert E[0] == E[1]
    R = reduction['common28_linear_target_rows']+[reduction['common29th_target_row']]
    odd = tc['common10_odd_target_rows']
    assert len(R) == 29 and len(odd) == 10
    J = s.Matrix(reduction['families']['CA']['linear_J'])
    O = J[:, :19]
    OR = O[R, :]
    assert E[0]*OR == O
    pivot = list(OR.T.rref()[1])
    complement = [i for i in range(29) if i not in pivot]
    assert complement == [6, 13, 14, 16, 17, 18, 19, 20, 27, 28]
    Q = s.zeros(29, 10)
    for j, i in enumerate(complement):
        Q[i, j] = 1
    selector = s.zeros(19, 29)
    for i, j in enumerate(pivot):
        selector[i, j] = 1
    projection = (s.eye(29)-OR*OR[pivot, :].inv()*selector)[complement, :]
    assert projection*OR == s.zeros(10, 19) and projection*Q == s.eye(10)
    nu = s.Matrix(first['CA_CB_discriminant_differentials'])
    bb = [s.Matrix(first['CA_CB_full_firstjet_data'][name]['actual'])[odd, :] for name in ['CA', 'CB']]
    assert nu*OR == s.zeros(2, 19) and all(b*OR == s.zeros(10, 19) for b in bb)
    variables = s.symbols('r0:10')
    linear = lambda row: sum(a*v for a, v in zip(row, variables))
    ns = [linear(row) for row in (nu*Q).tolist()]
    bs = [[linear(row) for row in (b*Q).tolist()] for b in bb]
    L = 2*variables[0]-variables[9]
    active = [0, 1, 2, 4, 6, 7, 9]
    qmons = [variables[i]*variables[j] for i, j in itertools.combinations_with_replacement(active, 2)]
    pairs = list(itertools.combinations_with_replacement(range(10), 2))
    cubics, quotients, nonzero_pairs = [], [], []
    for i, j in pairs:
        k = s.expand(ns[1]*bs[1][i]*bs[1][j]-ns[0]*bs[0][i]*bs[0][j])
        cubics.append(k)
        if not k:
            continue
        quotient, rem = s.div(k, L, *variables)
        assert rem == 0
        assert all(s.diff(quotient, variables[i]) == 0 for i in [3, 5, 8])
        quotients.append(s.expand(quotient))
        nonzero_pairs.append([i, j])
    assert len(quotients) == 28
    M = s.Matrix([[s.Poly(q, *variables).coeff_monomial(m) for m in qmons] for q in quotients])
    det = s.QQ.to_sympy(DomainMatrix.from_Matrix(M).convert_to(s.QQ).det())
    assert det != 0
    target = s.Matrix([s.Poly(L*L, *variables).coeff_monomial(m) for m in qmons])
    combo = M.T.inv()*target
    assert s.expand(sum(a*q*L for a, q in zip(combo, quotients))-L**3) == 0
    L10 = s.Matrix([[2]+[0]*8+[-1]])
    LR = L10*projection
    # Match the same hyperplane to the already actual common-R29 N projection.
    Jn = J[:, [i for i in range(28) if i != 20]]
    Wn = s.Matrix(first['N_full_firstjet_data']['W'])
    Y = reduction['common28_linear_target_rows']
    EN = Jn.row_join(Wn)
    EN = EN*EN[Y, :].inv()
    Npinch_linear = s.zeros(1, 29)
    Npinch_linear[0, 28] = 1
    for j, a in enumerate(EN[61, :]):
        Npinch_linear[0, j] -= a
    assert Npinch_linear == -LR
    plane_basis = s.Matrix.hstack(*LR.nullspace())
    assert plane_basis.shape == (29, 28)
    physical_plane = E[0]*plane_basis
    assert physical_plane.rank() == 28
    evidence = dict(
        common_R29_physical_coordinate_rows=R, common_odd10_physical_coordinate_rows=odd,
        common_even_embedding150_by29=mat(E[0]),
        common_odd10_linear_even_offsets=mat(E[0][odd, :]),
        orbit19_in_R29=mat(OR), orbit_complement_R29_indices=complement,
        orbit_complement_physical_rows=[R[i] for i in complement],
        orbit_quotient_projection10_by29=mat(projection),
        actual_discriminant_linear_forms2_by29=mat(nu),
        actual_odd_linear_forms_CA10_by29=mat(bb[0]), actual_odd_linear_forms_CB10_by29=mat(bb[1]),
        quotient10_variable_order=list(map(str, variables)), discriminant_linear_forms=list(map(str, ns)),
        odd_linear_forms_CA=list(map(str, bs[0])), odd_linear_forms_CB=list(map(str, bs[1])),
        all55_pure_even_cubic_projections=list(map(str, cubics)), nonzero_product_pairs=nonzero_pairs,
        common_linear_factor=str(L), active_linear_space_indices=active,
        all28_quotient_quadrics=list(map(str, quotients)), quotient_quadratic_monomials=list(map(str, qmons)),
        full28_by28_quotient_coefficient_matrix=mat(M), quotient_coefficient_determinant=str(det),
        exact_linear_combination_for_L_cube=list(map(str, combo)),
        full_R29_linear_support_equation=mat(LR), actual_N_pinch_projection_linear_form=mat(Npinch_linear),
        full150_by28_tangent_cone_support_plane_basis=mat(physical_plane),
    )
    with (out/'certificate.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', reduction_path, first_path,
             tangent_path, tangent_cert, *reviews]
    result = dict(
        id=STEM, all_passed=True,
        exact_common_even_linear_embedding_CA_equals_CB=True,
        full_actual_discriminant_and_odd_linear_forms_bound_in_common_R29=True,
        pure_even_projections_of_true_cubic_product_cancellations=55,
        nonzero_cubic_projection_count=28, cubic_projection_QQ_rank=28,
        exact_cubic_projection_ideal='L*(r0,r1,r2,r4,r6,r7,r9)^2', L='2*r0-r9',
        common_factor_equals_negative_actual_N_pinch_linear_projection=True,
        full_unknown_cubic_Y_terms_not_serialized_or_assumed_zero=True,
        actual_tangent_cone_reduced_support_is_single_linear28_plane=True,
        orbit_transverse_tangent_cone_reduced_support_is_single_linear9_plane=True,
        exact_radical_of_true_initial_ideal='(Y1,...,Y10,L)',
        tangent_cone_nilradical_sixth_power_zero=True,
        tangent_cone_nil_bound_from_shared_Yproducts_and_existence_of_Lcube_plus_Yquadratics=True,
        actual_nilradical_index_full_initial_ideal_exact_depth_or_vertex_absence_claimed=False,
        no_Groebner_basis_or_additional_source_jet_computed=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'certificate.json.gz': sha(out/'certificate.json.gz')},
        seconds=time.monotonic()-start,
    )
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'nonzero_cubic_projection_count',
          'exact_cubic_projection_ideal', 'tangent_cone_nilradical_sixth_power_zero', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())
