#!/usr/bin/env python3
"""One existing accepted NEW25 QQ point: exact quadratic-field biform test."""
import argparse
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

from weight1_structure import COORDS, J, constraints

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63'
sys.set_int_max_str_digits(0)


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class K:
    q = None

    def __init__(self, a=0, b=0):
        if isinstance(a, K):
            self.a, self.b = a.a, a.b
        else:
            if isinstance(a, bool):
                a = int(a)
            self.a, self.b = Fraction(str(a)), Fraction(str(b))

    def __add__(self, other):
        other = K(other)
        return K(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-K(other))

    def __rsub__(self, other):
        return K(other) - self

    def __mul__(self, other):
        other = K(other)
        return K(self.a * other.a + self.q * self.b * other.b,
                 self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.a * self.a - self.q * self.b * self.b
        assert norm
        return K(self.a / norm, -self.b / norm)

    def __truediv__(self, other):
        return self * K(other).inverse()

    def conjugate(self):
        return K(self.a, -self.b)

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, other):
        other = K(other)
        return self.a == other.a and self.b == other.b

    def save(self):
        return [str(self.a), str(self.b)]


def transpose(A):
    return list(map(list, zip(*A)))


def mm(A, B):
    return [[sum((x * y for x, y in zip(row, col)), K(0))
             for col in zip(*B)] for row in A]


def eye(n):
    return [[K(i == j) for j in range(n)] for i in range(n)]


def inverse(A):
    n = len(A)
    M = [list(row) + e for row, e in zip(A, eye(n))]
    for i in range(n):
        pivot = next(j for j in range(i, n) if M[j][i])
        M[i], M[pivot] = M[pivot], M[i]
        scale = M[i][i].inverse()
        M[i] = [x * scale for x in M[i]]
        for j in range(n):
            if j != i:
                scale = M[j][i]
                M[j] = [x - scale * y for x, y in zip(M[j], M[i])]
    assert [row[:n] for row in M] == eye(n)
    return [row[n:] for row in M]


def determinant4(A):
    answer = K(0)
    for perm in itertools.permutations(range(4)):
        term = K((-1) ** sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4)))
        for i in range(4):
            term *= A[i][perm[i]]
        answer += term
    return answer


def save_matrix(A):
    return [[x.save() for x in row] for row in A]


def mat(A):
    return [[K(x) for x in row] for row in A.tolist()]


def circuit_evaluate(payload, inputs):
    values = []
    for op, arg in payload['nodes']:
        if op == 'q':
            x = K(arg)
        elif op == 'v':
            x = inputs[arg]
        elif op == '+':
            x = sum((values[i] for i in arg), K(0))
        elif op == '*':
            x = K(1)
            for i in arg:
                x *= values[i]
        elif op == 'inv':
            x = values[arg].inverse()
        else:
            raise ValueError(op)
        values.append(x)
    return values


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    ip = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    samplespath = ROOT / 'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    acceptedpath = ROOT / 'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/report.json'
    spinpath = ROOT / 'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    discpath = ROOT / 'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63/mixed-discriminant-circuit.json.gz'
    reviewpath = ROOT / 'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW/report.json'
    info, sampledata, oldaccepted, spin, payload, reviewed = map(read, [ip, samplespath, acceptedpath, spinpath, discpath, reviewpath])
    assert oldaccepted['all_passed'] and reviewed['all_passed']
    # Deterministic first accepted point. No search over the 36-point packet.
    sample = sampledata['actual_QQ_source_points'][0]
    assert s.Rational(sample['value_pivot']) and s.Rational(sample['value_pfaffian'])
    B = s.Matrix(info['primitive24_by16_basis'])
    beta = B * s.Matrix([s.Rational(x) for x in sample['beta16']])
    ps = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    psi = s.Matrix(5, 4, lambda i, j: beta[COORDS.index((ps[i], j))])
    chi = s.Matrix([(-1) ** i * psi[[j for j in range(5) if j != i], :].det() for i in range(5)])
    Q = s.Matrix(spin['quadratic_form5'])
    normal = Q.inv() * chi
    q = (normal.T * Q * normal)[0]
    assert q < 0  # In particular nonzero and nonsquare over QQ for this fixed point.
    K.q = Fraction(str(q))
    p0 = (psi.T * Q * psi * J.inv()).det()
    assert q == 16 * p0
    gammas = [s.Matrix(a) for a in spin['gamma_matrices']]
    N = sum((normal[i] * gammas[i] for i in range(5)), s.zeros(4))
    assert N * N == q * s.eye(4)
    nmat, I, root = mat(N), eye(4), K(0, 1)
    projectors = [[[((I[i][j] + sign * nmat[i][j] / root) / K(2))
                    for j in range(4)] for i in range(4)] for sign in (1, -1)]
    omega = mat(J.inv())
    def pairing(a, b):
        return sum((a[i] * omega[i][j] * b[j] for i in range(4) for j in range(4)), K(0))
    eigenvectors, pair_indices, pairings = [], [], []
    for P in projectors:
        cols = transpose(P)
        i = next(i for i, col in enumerate(cols) if any(col))
        v = cols[i]
        j = next(j for j, col in enumerate(cols) if pairing(v, col))
        h = pairing(v, cols[j])
        w = [-x / h for x in cols[j]]  # omega_01 = omega_23 = -1.
        eigenvectors += [v, w]
        pair_indices.append([i, j])
        pairings.append(h)
    frame = transpose(eigenvectors)
    assert mm(mm(transpose(frame), omega), frame) == omega
    frame_inverse = inverse(frame)
    g = transpose(frame_inverse)
    assert mm(mm(transpose(g), mat(J)), g) == mat(J)
    assert determinant4(g) == K(1)
    diagonal = [[root * K((1 if i < 2 else -1) if i == j else 0) for j in range(4)] for i in range(4)]
    assert mm(mm(frame_inverse, nmat), frame) == diagonal
    beta_forms = []
    for k in range(4):
        F = s.zeros(4)
        for (i, j), kk in COORDS:
            if kk == k:
                F[i, j], F[j, i] = beta[COORDS.index(((i, j), k))], -beta[COORDS.index(((i, j), k))]
        beta_forms.append(mat(F))
    new_forms = []
    for k in range(4):
        h = [[sum((g[l][k] * beta_forms[l][i][j] for l in range(4)), K(0)) for j in range(4)] for i in range(4)]
        new_forms.append(mm(mm(transpose(g), h), g))
    new_beta = [new_forms[k][i][j] for (i, j), k in COORDS]
    new_psi = [[new_forms[k][i][j] for k in range(4)] for i, j in ps]
    assert new_psi[0] == [K(0)] * 4
    constraint = s.Matrix.vstack(*constraints())
    assert all(sum((K(constraint[i, j]) * new_beta[j] for j in range(24)), K(0)) == K(0)
               for i in range(8))
    assert determinant4(new_psi[1:]) == root
    for j in range(2):
        assert new_psi[3 + j][0] == new_psi[1 + j][1]
    for i in range(2):
        assert new_psi[2 + 2 * i][2] == new_psi[1 + 2 * i][3]
    uv = []
    for i in range(3):
        for j in range(2):
            uv.append(new_psi[1 + j][0] if i == 0 else
                      2 * new_psi[1 + j][1] if i == 1 else new_psi[3 + j][1])
    for i in range(2):
        uv += [new_psi[1 + 2 * i][2], 2 * new_psi[1 + 2 * i][3], new_psi[2 + 2 * i][3]]
    values = circuit_evaluate(payload, uv)
    conjugate_values = circuit_evaluate(payload, [x.conjugate() for x in uv])
    assert conjugate_values == [x.conjugate() for x in values]
    outputs = payload['outputs']
    R = values[outputs['resultant']]
    disc = values[outputs['quintic_discriminant']]
    D = values[outputs['mixed_discriminant']]
    assert R and D.b == 0 and D == conjugate_values[outputs['mixed_discriminant']]
    assert disc == D * R * R
    assert D.a.denominator % 32003
    cert = dict(sample_index=0, accepted_source_point=sample,
                beta24=list(map(str, beta)), psi5x4=[[str(x) for x in row] for row in psi.tolist()],
                normal_cofactor_chi=list(map(str, chi)), normal_vector5=list(map(str, normal)),
                normal_square_q=str(q), p0_detC=str(p0),
                field='QQ[s]/(s^2-q)', field_element_encoding='[rational_constant,rational_s_coefficient]',
                Clifford_normal4x4=[[str(x) for x in row] for row in N.tolist()],
                projectors=list(map(save_matrix, projectors)), eigenvector_projector_column_pairs=pair_indices,
                symplectic_pairings_before_normalization=[x.save() for x in pairings],
                eigenframe4x4=save_matrix(frame), physical_pullback_g4x4=save_matrix(g),
                transformed_beta24=[x.save() for x in new_beta], transformed_psi5x4=save_matrix(new_psi),
                biform12=[x.save() for x in uv],
                binary_quintic_coefficients=[values[i].save() for i in outputs['binary_quintic_coefficients']],
                V_quadratic_resultant=R.save(), binary_quintic_discriminant=disc.save(),
                mixed_discriminant=D.save(), all170_circuit_nodes=[x.save() for x in values],
                all170_field_conjugate_nodes_match=True)
    with gzip.open(out / 'certificate.json.gz', 'wt') as f:
        json.dump(cert, f, indent=2)
        f.write('\n')
    result = dict(id=STEM, all_passed=True, sample_index=0, number_existing_NEW25_points_tested=1,
                  field_quadratic_radicand_negative=True, nonisotropic_normal=True,
                  exact_symplectic_eigenframe_and_full_physical_beta_transform=True,
                  V_basepoint_resultant_nonzero=bool(R), mixed_discriminant_nonzero=bool(D),
                  mixed_discriminant_descends_to_QQ=D.b == 0,
                  mixed_discriminant_QQ_numerator_bits=D.a.numerator.bit_length(),
                  mixed_discriminant_QQ_denominator_bits=D.a.denominator.bit_length(),
                  mixed_discriminant_QQ_mod32003=D.a.numerator * pow(D.a.denominator, -1, 32003) % 32003,
                  field_involution_checked_on_all_circuit_nodes=True,
                  candidate_rejected_for_NEW25_hypersurface=bool(D),
                  vanishing_at_one_point_would_not_prove_equation=True,
                  no_degree22_or_higher_evaluation_matrix=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                      [Path(__file__), ROOT / 'code/weight1_structure.py', ip, samplespath, acceptedpath, spinpath, discpath, reviewpath]},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'mixed_discriminant_nonzero',
          'mixed_discriminant_descends_to_QQ', 'candidate_rejected_for_NEW25_hypersurface', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output.resolve())
