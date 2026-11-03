#!/usr/bin/env python3
"""Full octic coefficient relations and same-point two-adjoint recovery."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import sys
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

from q_new_family_degree16_invariant_probe_20261103_v63 import character
from weight1_structure import J, sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63'
sys.set_int_max_str_digits(0)


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):
    return [[str(x) for x in row] for row in M.tolist()]


def piv(M):
    return list(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])


def witness(M, rank):
    columns = piv(M)
    assert len(columns) == rank
    rr = piv(M[:, columns].T)
    W = M.extract(rr, columns)
    determinant = W.det(method='domain-ge')
    assert determinant
    return dict(rank=rank, rows=rr, columns=columns, minor=rows(W), determinant=str(determinant))


def add(polys):
    out = defaultdict(int)
    for poly in polys:
        for m, c in poly.items():
            out[m] += c
    return {m: c for m, c in out.items() if c}


def multiply(A, B):
    out = defaultdict(int)
    for m, a in A.items():
        for n, b in B.items():
            out[m + n] += a * b
    return {m: c for m, c in out.items() if c}


def polynomial_mm(A, B):
    return [add([multiply(A[4 * i + k], B[4 * k + j]) for k in range(4)])
            for i in range(4) for j in range(4)]


def trace_pair(A, B):
    return add([multiply(A[4 * i + j], B[4 * j + i]) for i in range(4) for j in range(4)])


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    spinpath = ROOT / 'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    pointpath = ROOT / 'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63/certificate.json.gz'
    invpath = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    spin, point, info = map(read, [spinpath, pointpath, invpath])
    monomials = spin['quadratic_monomials']
    powers = [9 ** i for i in range(16)]
    def integer_entries(key):
        out = []
        for row in spin[key]:
            poly = {}
            for i, value in enumerate(row):
                if value != '0':
                    c = 4 * Fraction(value)
                    assert c.denominator == 1
                    poly[sum(powers[j] for j in monomials[i])] = c.numerator
            out.append(poly)
        return out
    # All polynomials below are scaled by 4^4=256, uniformly for all 12 columns.
    Cp, Dp = map(integer_entries, ['spin_adjoint_C_coefficients16x136', 'vector_adjoint_D_coefficients16x136'])
    CC, CD, DD = polynomial_mm(Cp, Cp), polynomial_mm(Cp, Dp), polynomial_mm(Dp, Dp)
    traces = [trace_pair(A, B) for A, B in [(CC, CC), (CC, CD), (CC, DD), (CD, CD), (CD, DD), (DD, DD)]]
    quadrics = [add([A[5 * i] for i in range(4)]) for A in [CC, CD, DD]]
    products = [multiply(quadrics[i], quadrics[j]) for i in range(3) for j in range(i, 3)]
    polynomials = traces + products
    keys = sorted(set().union(*(P.keys() for P in polynomials)))
    relation_rows = []
    for product_index, coefficient in [(0, -4), (1, 2), (2, -2), (3, -1), (4, 1)]:
        r = [0] * 12
        r[6 + product_index], r[11] = 1, coefficient
        relation_rows.append(r)
    relation_rows += [[1, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, -3],
                      [-1, 0, 6, 0, 12, 8, 0, 0, 0, 0, 0, 0]]
    assert len(piv(s.Matrix(relation_rows))) == 7
    assert all(sum(c * P.get(m, 0) for c, P in zip(relation, polynomials)) == 0
               for relation in relation_rows for m in keys)
    base, selected = {}, []
    for m in keys:
        row = [P.get(m, 0) % 32003 for P in polynomials]
        for pivot, vector in sorted(base.items()):
            z = row[pivot]
            row = [(a - z * b) % 32003 for a, b in zip(row, vector)]
        pivot = next((i for i, c in enumerate(row) if c), None)
        if pivot is None:
            continue
        z = pow(row[pivot], -1, 32003)
        base[pivot] = [c * z % 32003 for c in row]
        selected.append(m)
        if len(selected) == 5:
            break
    coefficient_witness = witness(s.Matrix([[P.get(m, 0) for P in polynomials] for m in selected]), 5)
    # Reuse exactly the single QQ point already used in the mixed-discriminant test.
    assert point['sample_index'] == 0
    psi = s.Matrix(point['psi5x4'])
    beta = s.Matrix(point['accepted_source_point']['beta16'])
    Q = s.Matrix(spin['quadratic_form5'])
    C = psi.T * Q * psi * J.inv()
    V = psi * J.inv() * psi.T * Q
    gs = sp_basis()
    vector_reps = [s.Matrix(A) for A in spin['action5_matrices']]
    R5 = s.Matrix.hstack(*(A.reshape(25, 1) for A in vector_reps))
    r5rows = spin['lie_lift_selected_rows']
    Dcoordinates = R5[r5rows, :].inv() * V.reshape(25, 1)[r5rows, :]
    D = sum((c * A.T for c, A in zip(Dcoordinates, gs)), s.zeros(4))
    Cprime = J.inv() * C * J
    L = s.Matrix(20, 20, lambda r, c:
                 V[r // 4, c // 4] * int(r % 4 == c % 4)
                 - Cprime[c % 4, r % 4] * int(r // 4 == c // 4))
    gammas = [s.Matrix(A) for A in spin['gamma_matrices']]
    L = L.col_join(s.Matrix(4, 20, lambda r, c: gammas[c // 4][r, c % 4]))
    assert L * psi.reshape(20, 1) == s.zeros(24, 1)
    recovery_witness = witness(L, 19)
    free_column = next(i for i in range(20) if i not in recovery_witness['columns'])
    psi_free = psi.reshape(20, 1)[free_column]
    assert psi_free
    normalized_kernel = psi.reshape(20, 1) / psi_free
    assert L * normalized_kernel == s.zeros(24, 1) and normalized_kernel[free_column] == 1
    G = s.Matrix.hstack(*(A.T.reshape(16, 1) for A in gs))
    grows = piv(G.T)
    Gleft = G[grows, :].inv()
    CCcoeff = s.Matrix(spin['spin_adjoint_C_coefficients16x136'])
    DDcoeff = s.Matrix(spin['vector_adjoint_D_coefficients16x136'])
    paircoeff = (Gleft * CCcoeff[grows, :]).col_join(Gleft * DDcoeff[grows, :])
    Jac = s.zeros(20, 16)
    for k, (i, j) in enumerate(monomials):
        Jac[:, i] += paircoeff[:, k] * beta[j]
        Jac[:, j] += paircoeff[:, k] * beta[i]
    jacobian_witness = witness(Jac, 16)
    def trword(word, variation=None):
        M = s.eye(4)
        for i, c in enumerate(word):
            M = M * (variation[1] if variation and variation[0] == i else {'C': C, 'D': D}[c])
        return s.trace(M)
    def derivative(word, letter, E):
        return sum(trword(word, (i, E)) for i, c in enumerate(word) if c == letter)
    qD = trword('DD')
    gradient = s.zeros(4, 20)
    for i, letter in enumerate(['C'] * 10 + ['D'] * 10):
        E = gs[i % 10].T
        dCC, dCD, dDD, dCCCC, dDDDD = [derivative(w, letter, E) for w in ['CC', 'CD', 'DD', 'CCCC', 'DDDD']]
        gradient[:, i] = s.Matrix([dCC - 2 * dDD, dCD + dDD,
                                  dCCCC + 4 * dDDDD - 6 * qD * dDD,
                                  -dCCCC + 6 * derivative('CCDD', letter, E)
                                  + 12 * derivative('CDDD', letter, E) + 8 * dDDDD])
    gradient_witness = witness(gradient, 4)
    assert gradient * Jac == s.zeros(4, 16)
    nonzero_entry_index = next(i for i, c in enumerate(list(C) + list(D)) if c)
    nonzero_entry = (list(C) + list(D))[nonzero_entry_index]
    characters = character(list(map(tuple, info['Cartan_weights'])), 24)
    degree22, degree24 = [characters[d]['invariant_dimension'] for d in [22, 24]]
    cert = dict(trace_words=['CCCC', 'CCCD', 'CCDD', 'CDCD', 'CDDD', 'DDDD'],
                product_word_pairs=[['CC', 'CC'], ['CC', 'CD'], ['CC', 'DD'], ['CD', 'CD'], ['CD', 'DD'], ['DD', 'DD']],
                polynomial_scale=256, exponent_encoding='sum exponent_i*9^i; total degree8',
                trace_polynomials=[[[m, c] for m, c in sorted(P.items())] for P in traces],
                product_polynomials=[[[m, c] for m, c in sorted(P.items())] for P in products],
                seven_complete_relation_rows=relation_rows, coefficient_selected_monomials=selected,
                coefficient_rank5_witness=coefficient_witness,
                same_point_index=0, actual_C4=rows(C), actual_D4=rows(D), actual_V5=rows(V),
                actual_psi5x4=rows(psi), linear_recovery24x20=rows(L), recovery_rank19_witness=recovery_witness,
                recovery_free_column=free_column, recovery_normalized_kernel_at_same_point=list(map(str, normalized_kernel)),
                pair_coordinate_basis4=list(map(rows, [A.T for A in gs])),
                full_pair_coordinate_coefficients20x136=rows(paircoeff), jacobian20x16=rows(Jac),
                jacobian_rank16_witness=jacobian_witness, four_constraints_gradient4x20=rows(gradient),
                constraint_gradient_rank4_witness=gradient_witness,
                scale_fixing_entry_in_C_then_D_row_major=nonzero_entry_index,
                scale_fixing_nonzero_entry=str(nonzero_entry), character_records_through24=characters)
    with gzip.open(out / 'certificate.json.gz', 'wt') as f:
        json.dump(cert, f, indent=2)
        f.write('\n')
    result = dict(id=STEM, all_passed=True, exact_field='QQ',
                  all_12_octic_polynomials_rank=5, full_exact_relation_space_dimension=7,
                  primitive_octic_relations_after_quartic_products=2,
                  same_existing_point_only=True, recovery24x20_rank=19, recovery_kernel_dimension=1,
                  pair_map20x16_jacobian_rank=16, four_constraint_gradient_rank=4,
                  gradient_times_pair_jacobian_zero=True,
                  generic_pair_map_degree=2, generic_fiber_signs_only=True,
                  rational_image_field_equals_sign_invariant_field=True,
                  rational_Sp4_invariant_fields_agree=True,
                  image_is_one_generically_smooth_codimension4_component=True,
                  global_four_constraint_locus_exhaustion_or_normality_claimed=False,
                  all_regular_invariant_generation_claimed=False,
                  character_only_Inv22_dimension=degree22, character_only_Inv24_dimension=degree24,
                  no_new_beta_points_or_degree22plus_evaluations=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                    [Path(__file__), ROOT / 'code/weight1_structure.py', ROOT / 'code/q_new_family_degree16_invariant_probe_20261103_v63.py', spinpath, pointpath, invpath]},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'recovery24x20_rank', 'pair_map20x16_jacobian_rank',
          'four_constraint_gradient_rank', 'character_only_Inv22_dimension', 'character_only_Inv24_dimension', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())
