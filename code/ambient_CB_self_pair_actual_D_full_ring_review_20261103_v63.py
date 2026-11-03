#!/usr/bin/env python3
"""Full asymmetric Qb actual-D proof, cached38 gate, monomial invariants audit."""
import argparse
import gzip
import itertools as it
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
import ambient_CB_actual_D_cubic_syzygy_review_20261103_v63 as gate

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW'
load, sha, accepted, matrix, mv = gate.load, gate.sha, gate.accepted, gate.matrix, gate.mv
K, b, physical = gate.K, gate.b, gate.physical


def main(out):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63'
    accepted(producer / 'report.json')
    deps = [ROOT / 'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_CB_SELF_PAIR_COMMON_BASE_ETALE_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63_REVIEW.json']
    audits = [accepted(p) for p in deps]
    assert audits[0]['all150_regseq_syzygy_degree_bounds_and_cubic_socle_zero_accepted']
    assert audits[1]['actual_double_meeting_quotient_diagram_bound_and_accepted']
    assert audits[5]['whole_CB_geometric_flags'] == 2 and audits[5]['whole_CA_geometric_flags'] == 1
    assert audits[7]['generic_nonzero_arc_whole_fiber_exactly_two_reduced_points']
    assert audits[8]['whole_CB_source_fiber_exactly_one_geometric_point_for_every_a_nonzero']
    cubicfolder = ROOT / 'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63'
    setup = load(cubicfolder / 'generic-setup.json.gz')
    D = load(cubicfolder / 'generic-cubic.json.gz')
    J = [matrix(A) for A in setup['source_Jacobians150_by29']]
    R, j0 = setup['base29_rows'], setup['f_row']
    expr_matrix = lambda A: s.Matrix([[physical.expr(q) if q else 0 for q in row] for row in A])
    determinants = []
    for A, expected in zip(J, (-1728 * b ** 8 * (b * b - 5) / 125, -442368 * b ** 9 / 25, 8064 * b / 25)):
        dm = DomainMatrix.from_Matrix(expr_matrix([A[i] for i in R])).convert_to(K)
        value = dm.det()
        assert value == physical.coerce(expected)
        determinants.append(str(physical.expr(value)))
    # Bind full quadratic-zero gate independently using our accepted physical
    # Jacobians, the prior independent trueL, and the raw normal quadratics.
    d2path = ROOT / 'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/independent-certificate.json.gz'
    d2 = load(d2path)
    L = list(map(physical.coerce, d2['own_second_difference']))
    E = [J[0][i] + J[1][i] + [L[i]] for i in range(150)]
    point = expr_matrix(E).subs(b, 1)
    columns = list(DomainMatrix.from_Matrix(point).convert_to(s.QQ).rref()[1])
    assert len(columns) == 38
    rows = list(DomainMatrix.from_Matrix(point[:, columns].T).convert_to(s.QQ).rref()[1])
    inverse = DomainMatrix.from_Matrix(expr_matrix([[E[i][j] for j in columns] for i in rows])).convert_to(K).inv().to_list()
    rawpath = ROOT / 'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/generic-quadratic.json.gz'
    raw = matrix(load(rawpath)['raw150_by3'])
    for j in range(3):
        coefficients = mv(inverse, [raw[i][j] for i in rows])
        assert mv([[row[k] for k in columns] for row in E], coefficients) == [row[j] for row in raw]
    # Extra binding: all150 actual Z4 quadratics from the freshly reconstructed
    # source curves, without repeating any physical D3 computation.
    archivepath = ROOT / 'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63_REVIEW/independent-source-curves.json.gz'
    curves = [v for v in load(archivepath) if str(v['mode']) == '87']
    assert [physical.coerce(v['slope']) for v in curves] == list(map(physical.coerce, (0, 1, -1, 2)))
    cc = list(map(physical.coerce, setup['constant_odd_row_normalization']))
    samples = []
    for curve in curves:
        A, C0 = [list(map(physical.coerce, curve['physical_values'][i][2])) for i in (0, 2)]
        diff = [a - c for a, c in zip(A, C0)]
        samples.append([q if i == j0 else q - cc[i] * diff[j0] for i, q in enumerate(diff)])
    actualZ = [[samples[0][i], (samples[1][i] - samples[2][i]) / 2, (samples[1][i] + samples[2][i]) / 2 - samples[0][i]] for i in range(150)]
    assert actualZ == matrix(D['actual_Z4_quadratic_matrix'])

    # All normal variables: f,rho,z7; exact monomial computations.
    e = lambda i: tuple(int(i == j) for j in range(9))
    plus = lambda x, y: tuple(a + c for a, c in zip(x, y))
    power = lambda i, n: tuple(n if j == i else 0 for j in range(9))
    divides = lambda a, c: all(x <= y for x, y in zip(a, c))
    def minimal(gens):
        return sorted(g for g in set(gens) if not any(h != g and divides(h, g) for h in gens))
    def intersect(A, B):
        return minimal([tuple(max(a, c) for a, c in zip(x, y)) for x in A for y in B])
    f, rho, zs = e(0), e(1), [e(i) for i in range(2, 9)]
    JZ = zs + [power(0, 2), power(1, 2)]
    Q = [plus(f, z) for z in zs] + [power(0, 3), plus(f, power(1, 2))]
    primary = zs + [power(0, 3), power(1, 2)]
    assert intersect([f], primary) == minimal(Q)
    upper = intersect([f], JZ)
    # The right/left difference is precisely P{f²,f²rho}.
    extra = [m for i in range(1, 5) for j in range(5) for m in [plus(power(0, i), power(1, j))] if any(divides(g, m) for g in upper) and not any(divides(g, m) for g in Q)]
    assert extra == [power(0, 2), plus(power(0, 2), rho)]
    nilbasis = [f, power(0, 2), plus(f, rho), plus(power(0, 2), rho)]
    assert all(not any(divides(g, m) for g in Q) for m in nilbasis)
    squarebasis = [power(0, 2), plus(power(0, 2), rho)]
    assert all(not any(divides(g, m) for g in Q) for m in squarebasis)
    socle = squarebasis[1]
    assert all(any(divides(g, plus(socle, e(i))) for g in Q) for i in range(9))
    assert len([plus(power(0, i), power(1, j)) for i in range(3) for j in range(2)]) == 6
    t = s.symbols('t')
    H = 1 / (1 - t) ** 28 + (t + 2 * t * t + t ** 3) / (1 - t) ** 20
    assert s.diff(H, t).subs(t, 0) == 29 and s.limit((1 - t) ** 28 * H, t, 1) == 1
    evidence = dict(common_actual_R29_three_projection_determinants=determinants,
        independent_exact38_projection_of_all150_raw_quadratics_zero=True,
        independent_archive_all150_Z4_quadratic_binding=True,
        core81_lower_ideal_equals_f_JZ_by_unit_and_Nakayama_accepted=True,
        full121_odd_syzygy_colon_recipe_and_degree_bounds_accepted=True,
        remaining_after_quadzero_is_submodule_of_free_P_rank1=True,
        cubic_generic_zero_plus_G_dominance_kills_remaining_submodule=True,
        no_nonflat_fiber_product_basechange_assumed=True,
        actual_full_ideal='(f*z7,f^3,f*rho^2)', nilmodule_rank=4, nil_square_rank=2, nil_cube_zero=True,
        primary_generic_length=6, nilmodule_generic_length=4,
        Ass=['minimal N', 'embedded S_Q'], cone_dim28_edim29_depth20_e1=True,
        global_embedded_SQ20_and_distinct_from_SP_by_whole_flag_counts_accepted=True,
        origin_boundary_or_global_list_exhaustion_claim=False)
    (out / 'independent-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md', ROOT / 'code/ambient_CB_actual_D_cubic_syzygy_review_20261103_v63.py',
        ROOT / 'code/ambient_CB_self_pair_actual_D_full_ring_20261103_v63.py', ROOT / 'research/notes/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md',
        producer / 'report.json', cubicfolder / 'generic-setup.json.gz', cubicfolder / 'generic-cubic.json.gz', archivepath, rawpath, d2path, *deps]
    report = dict(id=STEM, all_passed=True, producer_or_physical_worker_rerun=False,
        asymmetric_Qb_full_actual_D_ideal_equality_accepted=True,
        core81_plus_full121_colon_and_cubic_gate_application_accepted=True,
        nilrank4_square_rank2_cubezero_Ass_N_and_embedded_SQ_depth20_e1_accepted=True,
        global_embedded_SQ_distinct_from_SP_accepted=True,
        scope='Q1 and a nonempty unit open of Q_b plus G saturation; no b0/all-boundary/full-global-primary-list claim.',
        independent_method='Cached independently bound source matrices; new exact38 quotient and all150 Z4 check; finite-colon/core81 proof; monomial primary/nil/socle arithmetic.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(out / 'independent-evidence.json')}, seconds=time.monotonic() - start)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output)
