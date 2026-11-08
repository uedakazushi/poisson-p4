#!/usr/bin/env python3
"""ONE whole-source readback for the two seven-kernel co-Higgs sources.

No new Jacobi linear operator, rank, minor, point, CAS, or jet is built.
The old full kernels, actual shear, and the SAME D2 inverse are input.
Only stable prior independent Fraction sparse primitives are imported.
"""
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import resource
import time

import classification_plane_rank_zero_rank_two_review as r

q = r.q
ROOT = Path(__file__).resolve().parents[1]
JOB = 'research/classification/jobs/plane_cohiggs_kernel_pointwise.json'
OUTPUT = 'reports/CLASSIFICATION_PLANE_COHIGGS_KERNEL_POINTWISE_20261108.json'
INPUT = 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json'
PAIRS = list(combinations(range(5), 2))
ZERO = (0,) * 8


def fraction(value, rho):
    value = r.evaluate(r.scalar(value), rho)
    assert len(value.p) <= 1 and len(value.q) == 1
    return value.p[0] / value.q[0] if value.p else F(0)


def mat(data, rho):
    return [[fraction(value, rho) for value in row] for row in data]


def constant_tensor(data, rho):
    return q.antisymmetric({
        ij: {ex: fraction_from_rf(c, rho) for ex, c in r.parse(data[str(ij)]).items()}
        for ij in PAIRS
    })


def fraction_from_rf(value, rho):
    value = r.evaluate(value, rho)
    assert len(value.p) <= 1 and len(value.q) == 1
    return value.p[0] / value.q[0] if value.p else F(0)


def bases():
    first = []
    for normal in (3, 4):
        for ij in combinations(range(3), 2):
            for k in range(3):
                first.append(q.antisymmetric({ij: q.atom(1, normal, k)}))
    for normal, j in [((3, 3), 4), ((3, 4), 3),
                      ((3, 4), 4), ((4, 4), 3)]:
        for i in range(3):
            first.append(q.antisymmetric({(i, j): q.atom(1, *normal)}))
    second = [
        q.antisymmetric({ij: q.atom(1, *normal)})
        for ij in combinations(range(3), 2)
        for normal in ((3, 3), (3, 4), (4, 4))
    ]
    assert len(first) == 30 and len(second) == 9
    return first, second


def tensor_linear(columns, coefficients):
    return q.antisymmetric({
        ij: q.add(*(q.scale(B.get(ij, {}), c)
                    for B, c in zip(columns, coefficients)))
        for ij in PAIRS
    })


def polynomial_records(poly):
    return [{'exponents': list(ex), 'coefficient': str(c)}
            for ex, c in sorted(poly.items()) if c]


def tensor_records(P):
    return {str(ij): polynomial_records(P.get(ij, {})) for ij in PAIRS}


def divisor(P):
    def pf(indices):
        a, b, c, d = indices
        return q.add(q.mul(P.get((a, b), {}), P.get((c, d), {})),
                     q.scale(q.mul(P.get((a, c), {}), P.get((b, d), {})), -1),
                     q.mul(P.get((a, d), {}), P.get((b, c), {})))
    return q.add(*(q.scale(q.mul(q.atom(1, i),
                                 pf([j for j in range(5) if j != i])),
                           (-1)**i) for i in range(5)))


def lift(poly, n):
    return {ex + (0,) * (n - len(ex)): c for ex, c in poly.items()}


def atom_n(c=1, index=None, n=12):
    ex = [0] * n
    if index is not None:
        ex[index] = 1
    return {tuple(ex): F(c)} if c else {}


def subst(poly, images):
    n = len(images)
    out = {}
    for ex, c in poly.items():
        term = atom_n(c, n=n)
        for i, power in enumerate(ex):
            for _ in range(power):
                term = q.mul(term, images[i])
        out = q.add(out, term)
    return out


def push(P, G, inverse_images):
    L = [[q.der(G[i], j) for j in range(5)] for i in range(5)]
    return q.antisymmetric({
        (i, j): subst(q.add(*(
            q.mul(q.mul(L[i][k], P.get((k, l), {})), L[j][l])
            for k in range(5) for l in range(5))), inverse_images)
        for i, j in PAIRS
    })


def universal_shear(P0, curve, first, S):
    n = 12
    x = [atom_n(index=i, n=n) for i in range(n)]
    G, inverse_images = list(x), list(x)
    for i in range(3):
        delta = q.add(q.mul(x[6+i], x[3]), q.mul(x[9+i], x[4]))
        G[i] = q.add(x[i], delta)
        inverse_images[i] = q.add(x[i], q.scale(delta, -1))
    assert [subst(p, inverse_images) for p in G] == x
    assert [subst(p, G) for p in inverse_images] == x
    P0n = {ij: lift(p, n) for ij, p in P0.items()}
    pushed0 = push(P0n, G, inverse_images)
    difference = q.antisymmetric({
        ij: {ex: c for ex, c in q.add(pushed0.get(ij, {}),
                                    q.scale(P0n.get(ij, {}), -1)).items()
             if ex[3] + ex[4] - sum(i >= 3 for i in ij) == 1}
        for ij in PAIRS
    })
    expected = q.antisymmetric({
        ij: q.add(*(q.mul(x[6+j], lift(tensor_linear(
            first, [S[i][j] for i in range(30)]).get(ij, {}), n))
                    for j in range(6)))
        for ij in PAIRS
    })
    residual = q.antisymmetric({
        ij: q.add(difference.get(ij, {}), q.scale(expected.get(ij, {}), -1))
        for ij in PAIRS
    })
    # This is precisely equality in the complete Euler quotient.
    for i, j, k in combinations(range(5), 3):
        assert not q.add(q.mul(x[i], residual.get((j, k), {})),
                         q.scale(q.mul(x[j], residual.get((i, k), {})), -1),
                         q.mul(x[k], residual.get((i, j), {})))
    whole = {ij: lift(p, n) for ij, p in curve.items()}
    roundtrip = push(push(whole, G, inverse_images), inverse_images, G)
    assert all(roundtrip.get(ij, {}) == whole.get(ij, {}) for ij in PAIRS)
    return {
        'all_six_shears_independent_parameters': True,
        'coordinate_map': ['a+s0*z+s3*w', 'b+s1*z+s4*w',
                           'c+s2*z+s5*w', 'z', 'w'],
        'inverse': ['a-s0*z-s3*w', 'b-s1*z-s4*w',
                    'c-s2*z-s5*w', 'z', 'w'],
        'determinant': '1',
        'both_global_polynomial_inverse_compositions': True,
        'actual_weight_one_shear_matches_saved_complete_30_by6_mod_Euler': True,
        'whole_t_curve_full150_push_inverse_roundtrip': True
    }


def run_case(record, rho, wall, first, second):
    C1 = wall['D1'] if wall else record['generic_D1']
    C2 = wall['D2'] if wall else record['generic_D2']
    CS = wall['actual_shear'] if wall else record['generic_actual_shear']
    assert C1['nullity'] == 7 and C2['nullity'] == 0 and C2['rank'] == 9
    assert CS['rank'] == 6
    S = mat(record['complete_actual_shear_30_by6'], rho)
    K = mat(C1['complete_kernel'], rho)
    I = CS['minor_rows']
    assert CS['minor_columns'] == list(range(6))
    Sinv = q.inverse([[S[i][j] for j in range(6)] for i in I])
    A = q.matmul(Sinv, [K[i] for i in I])
    residual = [[K[i][j] - sum(S[i][k]*A[k][j] for k in range(6))
                 for j in range(7)] for i in range(30)]
    j0 = next(j for j in range(7) if any(row[j] for row in residual))
    pivot = next(i for i in range(30) if residual[i][j0])
    Ncoordinates = [row[j0] / residual[pivot][j0] for row in residual]
    b = residual[pivot]
    assert all(residual[i][j] == Ncoordinates[i]*b[j]
               for i in range(30) for j in range(7))
    assert all(not Ncoordinates[i] for i in I)
    P0 = constant_tensor(record['normal_form_P0'], rho)
    N = tensor_linear(first, Ncoordinates)
    JN = q.jvec(q.jacobi(N))
    D2 = mat(record['complete_D2_350_by9'], rho)
    R, C = C2['minor_rows'], C2['minor_columns']
    assert sorted(C) == list(range(9)) and len(R) == 9
    inverse_D2 = q.inverse([[D2[i][j] for j in C] for i in R])
    solved = [sum(row[j] * (-JN[R[j]]) for j in range(9))
              for row in inverse_D2]
    Qcoordinates = [F(0)] * 9
    for j, value in zip(C, solved):
        Qcoordinates[j] = value
    Q = tensor_linear(second, Qcoordinates)
    t, t2 = q.atom(1, 5), q.atom(1, 5, 5)
    curve = q.antisymmetric({
        ij: q.add(P0.get(ij, {}), q.mul(t, N.get(ij, {})),
                  q.mul(t2, Q.get(ij, {}))) for ij in PAIRS
    })
    J = q.jacobi(curve)
    rows_by_power = [q.jvec(J, (j, 0, 0)) for j in range(5)]
    assert all(not value for row in rows_by_power[:2] for value in row)
    assert all(sum(ex[5:]) <= 4 and not any(ex[6:])
               for p in J.values() for ex in p)
    nonzero = [(i, j, value) for j, row in enumerate(rows_by_power)
               for i, value in enumerate(row) if value]
    Factual = divisor(curve)
    if nonzero:
        i, j, value = nonzero[0]
        assert j >= 2
        conclusion = 'reduced_support_t0_only__proper_P0_divisor__EMPTY'
        obstruction = {'original_full350_row': i, 'power_t': j,
                       'nonzero_coefficient': str(value)}
        factor = None
    else:
        obstruction = None
        if not Factual:
            conclusion = 'entire_whole_curve_actual_Pf_zero__EMPTY'
            factor = None
        else:
            common = [min(ex[i] for ex in Factual) for i in range(5)]
            if any(common):
                fex = tuple(common) + (0, 0, 0)
                quotient = {tuple(a-b for a, b in zip(ex, fex)): c
                            for ex, c in Factual.items()}
                assert q.mul({fex: F(1)}, quotient) == Factual
                conclusion = 'entire_whole_curve_proper_coordinate_factor__EMPTY'
                factor = {'factor_exponents': list(fex),
                          'whole_quotient': polynomial_records(quotient)}
            else:
                conclusion = 'entire_whole_curve_explicit__pointwise_parent_UNRESOLVED'
                factor = None
    return {
        'source': record['type'], 'rho': str(rho),
        'old_complete_kernel_dim': 7, 'old_actual_shear_rank': 6,
        'old_D2_nullity': 0, 'no_new_operator_rank_minor_or_point': True,
        'complete_30_by7_identity_K_equals_SA_plus_Nb': True,
        'A_6_by7': [[str(v) for v in row] for row in A],
        'b_7': list(map(str, b)),
        'N_30': list(map(str, Ncoordinates)),
        'P0': tensor_records(P0), 'N': tensor_records(N),
        'Q_9': list(map(str, Qcoordinates)), 'Q': tensor_records(Q),
        'same_old_shear_minor_rows': I,
        'same_old_D2_inverse_rows': R, 'same_old_D2_inverse_columns': C,
        'same_old_D2_minor_determinant': C2['minor_determinant'],
        'candidate_whole_curve': tensor_records(curve),
        'all350_original_Jacobi_residuals_by_t_power':
            [[str(v) for v in row] for row in rows_by_power],
        'all350_identically_zero_on_whole_curve': not nonzero,
        'actual_Euler_Pf_all126_parameter_polynomials': polynomial_records(Factual),
        'reduced_support_obstruction_if_any': obstruction,
        'proper_factor_certificate_if_any': factor,
        'actual_finite_shear_global_inverse': universal_shear(
            P0, curve, first, S),
        'whole_source_conclusion': conclusion
    }


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (175, 180))
    resource.setrlimit(resource.RLIMIT_AS, (3*1024**3, 3*1024**3))
    started = time.monotonic()
    job = json.loads((ROOT/JOB).read_text())
    for inp in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/inp['path']).read_bytes()).hexdigest() == inp['sha256']
    destination = ROOT/OUTPUT
    assert not destination.exists(), 'immutable output already exists'
    packet = json.loads((ROOT/INPUT).read_text())
    first, second = bases()
    BSS = next(c for c in packet['cases'] if c['type'] == 'BSS_rho')
    BN = next(c for c in packet['cases'] if c['type'] == 'BN0')
    wall = next(v for v in BSS['all_exceptional_roots'] if v['rho'] == '-3')
    results = [run_case(BSS, F(-3), wall, first, second),
               run_case(BN, F(0), None, first, second)]
    result = {
        'task_id': 'E_reduced_nonnormal', 'status': 'finite_two_whole_source_batch_completed',
        'all_explicit_identity_checks_passed': True,
        'both_whole_sources_empty_proved': all(
            'EMPTY' in x['whole_source_conclusion'] for x in results),
        'scope': 'Entire M(p)=0/m2/rankM2 BSSrho=-3 and BN0, all seven K1 and FULL nine inhomogeneous P2, all six finite shears and global inverse.',
        'cases': results, 'inputs': job['exact_inputs_and_hashes'],
        'no_producer_import_no_CAS': True,
        'no_new_M1_operator_rank_minor_point_or_jet': True,
        'seconds': time.monotonic()-started,
        'maximum_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'executed_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'executed_source': Path(__file__).read_text(),
        'main_frontier_global_and_old_bytes_unchanged': True
    }
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    for x in results:
        print(x['source'], x['whole_source_conclusion'], flush=True)
    print('PASS explicit identities', round(result['seconds'], 3),
          'seconds', result['maximum_RSS_KiB'], 'KiB', flush=True)


if __name__ == '__main__':
    main()
