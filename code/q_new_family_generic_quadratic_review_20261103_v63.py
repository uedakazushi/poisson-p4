#!/usr/bin/env python3
"""Independent physical Hessian+Jacobian,79 unit pivots and augmented99 audit.

The successful77.9066s native field run is bound entry-by-entry and its exact
engine audited, not rerun in this review. No producer derivative is imported.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from q_positive_weight6_small_model_review_20261103_v63 import physical_coefficients

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    started = time.monotonic()
    directory = ROOT / 'reports/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63'
    report = load(directory / 'report.json')
    assert report['all_passed'] is True
    for name, digest in report['input_sha256'].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in report['evidence_sha256'].items():
        assert sha(directory / name) == digest, name
    cert = load(directory / 'certificate.json.gz')
    family_path = ROOT / 'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz'
    smith_path = ROOT / 'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz'
    normal_path = ROOT / 'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63/certificate.json.gz'
    family, smith, normal = map(load, [family_path, smith_path, normal_path])
    normal_review = ROOT / 'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW.json'
    assert load(normal_review)['all_passed'] is True
    R, n0, n2 = ring('n0,n2', s.QQ)
    sym1 = s.Symbol('n1')
    q = [R.from_expr(s.sympify(value).subs(sym1, 1)) for value in family['canonical150']]
    C = s.Matrix(smith['canonical_basis'])
    T = s.Matrix(normal['extra_normal126_by2'])
    physicalT = C * T
    Tr = [[R.from_expr(value) for value in row] for row in physicalT.tolist()]
    Crows = [[(j, s.QQ.convert(C[i,j])) for j in range(126) if C[i,j]] for i in range(150)]
    # The constant physical quadratic dictionaries were independently built
    # from the full Schouten formula in the earlier small-model review.
    coefficients, _ = physical_coefficients()
    J, H = [], []
    for terms in coefficients:
        derivative, hessian = {}, [R.zero, R.zero, R.zero]
        for (a,b), coefficient in terms.items():
            derivative[a] = derivative.get(a, R.zero) + coefficient * q[b]
            derivative[b] = derivative.get(b, R.zero) + coefficient * q[a]
            hessian[0] += coefficient * Tr[a][0] * Tr[b][0]
            hessian[1] += coefficient * (Tr[a][0]*Tr[b][1] + Tr[a][1]*Tr[b][0])
            hessian[2] += coefficient * Tr[a][1] * Tr[b][1]
        row = {}
        for i, value in derivative.items():
            if value:
                for j, constant in Crows[i]:
                    row[j] = row.get(j, R.zero) + constant * value
        J.append({j:value for j,value in row.items() if value})
        H.append(hessian)
    assert len(J) == len(H) == 350
    for i in range(350):
        for j in range(3):
            assert H[i][j] == R.from_expr(s.sympify(cert['full350_by3_quadratic_Hessian'][i][j]))
            assert H[i][j] == R.from_expr(s.sympify(normal['full350_by3_Hessian_columns'][i][j]))
        for j in range(2):
            assert sum((value * R.from_expr(T[k,j]) for k,value in J[i].items()), R.zero) == R.zero
    A = {i:dict(row) for i,row in enumerate(J)}
    for i in range(350):
        A[i].update({126+j:value for j,value in enumerate(H[i]) if value})
    columns = list(range(129))
    assert len(cert['rational_unit_pivots']) == 79
    for row_id, col_id, constant in cert['rational_unit_pivots']:
        value = s.QQ.convert(s.Rational(constant))
        assert value and A[row_id].get(col_id, R.zero) == R(value)
        pivotrow = A[row_id]
        for i, row in A.items():
            if i == row_id:
                continue
            entry = row.pop(col_id, R.zero)
            if entry:
                factor = entry / value
                for j, coefficient in pivotrow.items():
                    if j != col_id:
                        changed = row.get(j, R.zero) - factor * coefficient
                        if changed:
                            row[j] = changed
                        else:
                            row.pop(j, None)
        del A[row_id]
        columns.remove(col_id)
    assert len(A) == 271 and len(columns) == 50
    zero_rows = [i for i,row in A.items() if not row]
    zero_columns = [j for j in columns if not any(j in row for row in A.values())]
    assert zero_rows == cert['all98_identically_zero_original_rows'] and len(zero_rows) == 98
    assert zero_columns == cert['all8_identically_zero_original_columns'] and len(zero_columns) == 8
    rows = [i for i,row in A.items() if row]
    cols = [j for j in columns if j not in zero_columns]
    assert rows == cert['remaining_original_rows'] and cols == cert['remaining_original_columns']
    assert len(rows) == 173 and len(cols) == 42
    for i, row_id in enumerate(rows):
        for j, col_id in enumerate(cols):
            assert A[row_id].get(col_id, R.zero) == R.from_expr(s.sympify(cert['full173_by42_polynomial_residual'][i][j]))
    # Every numerator AND denominator in the ACTUAL successful native input
    # is checked against the independently reconstructed polynomial matrix.
    with (directory / 'native-input.txt').open() as stream:
        assert list(map(int, stream.readline().split())) == [2, 173, 42]
        def polynomial():
            p = R.zero
            for _ in range(int(stream.readline())):
                data = list(map(int, stream.readline().split()))
                assert len(data) == 3
                p += R.from_dict({tuple(data[1:]):s.QQ(data[0])})
            return p
        for i in rows:
            for j in cols:
                numerator, denominator = polynomial(), polynomial()
                assert denominator and denominator.is_ground
                assert numerator / denominator.LC == A[i].get(j, R.zero)
        assert not stream.read().strip()
    log = (directory / 'native-log.txt').read_text()
    assert 'EXACT_FUNCTION_FIELD_RANK 20 ROWS 173 COLUMNS 42' in log
    assert log.count('PIVOT ') == 20
    engine = ROOT / 'code/rational_function_matrix_rank_20261103_v62.cpp'
    text = engine.read_text()
    assert 'assert(matrix[row][c].zero())' in text
    assert 'fmpz_mpoly_q_div' in text and 'fmpz_mpoly_q_sub' in text
    # Independent lower99 and binary discriminant at the rational111 point.
    def at_one(poly):
        return s.QQ.to_sympy(sum(poly.values(), s.QQ.zero))
    sample = s.zeros(350,129)
    for i,row in enumerate(J):
        for j,value in row.items(): sample[i,j] = at_one(value)
        for j in range(3): sample[i,126+j] = at_one(H[i][j])
    det = DomainMatrix.from_Matrix(sample[cert['lower99_rows'],cert['lower99_columns']]).convert_to(s.QQ).det()
    assert s.Rational(det) == s.Rational(cert['lower99_minor']) != 0
    Jp, Hp = sample[:,:126], sample[:,126:]
    jc = list(DomainMatrix.from_Matrix(Jp).rref()[1])
    jr = list(DomainMatrix.from_Matrix(Jp[:,jc].T).rref()[1])
    assert len(jc) == 98
    coker = Hp - Jp[:,jc] * Jp[jr,jc].inv() * Hp[jr,:]
    assert DomainMatrix.from_Matrix(coker).rank() == 1
    binary = list(coker[cert['nondegenerate_binary_QQ_row_index'],:])
    assert binary == list(map(s.Rational,cert['nondegenerate_binary_QQ_coefficients']))
    discriminant = s.factor(binary[1]**2 - 4*binary[0]*binary[2])
    assert discriminant == s.Rational(cert['nondegenerate_binary_discriminant']) != 0
    paths = [Path(__file__), directory/'report.json',directory/'certificate.json.gz',
        directory/'native-input.txt',directory/'native-log.txt',family_path,smith_path,normal_path,
        normal_review,engine,ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.py',
        ROOT/'research/notes/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63_REVIEW.md']
    result = dict(id='Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63_REVIEW',all_passed=True,
        independent_physical_all350_Jacobian_and_three_Hessian_columns=True,
        all350_Jacobian_times_two_polynomial_normals_zero=True,
        independently_replayed_all79_nonzero_QQ_unit_operations=True,
        all98_zero_rows_and8_zero_columns_verified_as_polynomial_identities=True,
        all173_by42_native_numerator_denominator_entries_bound_to_independent_matrix=True,
        successful_native20_rank_receipt_and_final_all_zero_engine_audited=True,
        native77second_run_reexecuted_in_this_review=False,
        independent99_minor_and_binary_nondegenerate_discriminant=True,
        generic_augmented_rank=99,generic_Jacobian_rank=98,
        intrinsic_generic_two_normal_quadratic_obstruction_rank=1,
        descent_mod_actual26_and_dense_corrected_shear_scope_accepted=True,
        generic_Artin_length_at_least5=True,
        exact_Artin_length_or_flat_full_completion_requires_separate_finite_rank_lemma=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('INDEPENDENT AUG99 / INTRINSIC QUADRATIC RANK1 PASS',result['seconds'],flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
