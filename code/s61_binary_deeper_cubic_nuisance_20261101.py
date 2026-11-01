#!/usr/bin/env python3
"""Tracked cubic constraints on the deeper seven-variable binary chart.

The JSON sidecar supplies candidate rational coefficients, which are checked
against the reconstructed quadratic rows. FLINT performs exact QQ polynomial
and QQ(e,d) field arithmetic. No ignored work file is mathematical input.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement

import numpy as np
import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_deeper_uniform_quadratic_20261101 import (
    build_quadratics, elimination_rank, PIVOT_ROWS, PIVOT_COLS,
)

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / 'code/certificates/s61_binary_deeper_cubic_20261101.json'
LINEAR = ROOT / 'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
QUADRATIC = ROOT / 'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
SUPPORT = [0, 1, 4, 7, 8, 9, 12, 18, 21, 25, 26, 27, 30, 32, 34, 35, 40]
MS2 = list(combinations_with_replacement(range(9), 2))
MS3 = list(combinations_with_replacement(range(9), 3))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normal_exponent(indices):
    result = [0] * 9
    for index in indices:
        result[index] += 1
    return result


def write_poly(stream, terms):
    stream.write(str(len(terms)) + '\n')
    for monomial, coefficient in terms:
        assert len(monomial) == 16
        stream.write(str(coefficient) + ' ' + ' '.join(map(str, monomial)) + '\n')


def write_core_poly(stream, poly):
    stream.write(str(len(poly)) + ' ' + ' '.join(
        str(coefficient) + ' ' + ' '.join(map(str, monomial))
        for monomial, coefficient in poly.items()) + '\n')


def prepare(data, certificate, native_input, module_input):
    started = time.time()
    R, quad = build_quadratics(data)
    K = s.QQ.frac_field('e', 'd')
    R2 = K.field.ring
    R5, *_ = ring('p,q,r,h,k', K)
    zero = R.zero
    core = lambda raw: R2.from_dict({tuple(m): s.QQ.from_sympy(s.Rational(c)) for m, c in raw})
    fraction = lambda raw: K.field.raw_new(core(raw['numerator']), core(raw['denominator']))
    outer = lambda raw: R5.from_dict({tuple(m): fraction(c) for m, c in raw})
    rr = certificate['quadratic_rref']
    assert rr['columns'] == PIVOT_COLS
    basis = [[K.from_sympy(s.sympify(x)) for x in row] for row in rr['rref']]
    assert len(basis) == 24
    for i, row in enumerate(basis):
        assert [row[j] for j in PIVOT_COLS] == [K.one if i == k else K.zero for k in range(24)]

    def convert(poly):
        grouped = {}
        for monomial, coefficient in poly.items():
            grouped.setdefault(monomial[2:], {})[monomial[:2]] = coefficient
        return R5.from_dict({m: K.field.raw_new(R2.from_dict(c), R2.one)
                             for m, c in grouped.items()})

    def clear_outer(poly, denominator):
        result = {}
        for nuisance, coefficient in poly.items():
            numerator = coefficient.numer * denominator.exquo(coefficient.denom)
            for monomial, value in numerator.items():
                result[monomial + nuisance] = value
        return R.from_dict(result)

    # Change the quadratic implicit lift only by multiples of quadratic
    # initial equations. Cubic classes change by their linear multiples.
    reduced_lift = []
    for original in quad['W']:
        row = [convert(value) for value in original]
        previous = row[:]
        for k, pivot in enumerate(PIVOT_COLS):
            if previous[pivot]:
                row = [a - previous[pivot].mul_ground(b) for a, b in zip(row, basis[k])]
        assert all(not row[j] for j in PIVOT_COLS)
        reduced_lift.append(row)
    assert all(not any(row) for i, row in enumerate(reduced_lift) if i not in SUPPORT)
    assert not any(reduced_lift[44])
    assert all(value == 0 for value in quad['Tn'][44])
    lift_denominator = R2.one
    for i in SUPPORT:
        for value in reduced_lift[i]:
            for coefficient in value.values():
                lift_denominator = lift_denominator.lcm(coefficient.denom)
    lift = [[clear_outer(value, lift_denominator) for value in reduced_lift[i]] for i in SUPPORT]
    print('quadratic lift reduced with zero t component', time.time() - started, flush=True)

    # Reconstruct the bilinear Hessian operators from the full original rows.
    model = S61()
    beta = {int(k[1:]) - 1: s.sympify(v) for k, v in data['beta'].items()}
    _, equations = specialize(model, beta, True)
    assert len(equations) == 135
    pol = lambda value: R.from_expr(s.sympify(value))
    J = [[pol(value) for value in row] for row in data['J']]
    T, U, Du = quad['Tn'], quad['Un'], quad['Uden']
    selected = data['pivot_rows']
    operators = []
    for equation in equations:
        matrix = [[zero] * len(SUPPORT) for _ in range(9)]
        for monomial, coefficient in equation.items():
            if len(monomial) != 2 or data['fixed'] in monomial:
                continue
            i, j = [data['cols'].index(v) for v in monomial]
            coefficient = pol(coefficient)
            for a in range(9):
                if i in SUPPORT and T[j][a]:
                    matrix[a][SUPPORT.index(i)] += coefficient * T[j][a]
                if j in SUPPORT and T[i][a]:
                    matrix[a][SUPPORT.index(j)] += coefficient * T[i][a]
        operators.append(matrix)
    projection = [[sum((c * U[j][k] for j, c in enumerate(row) if c), zero)
                   for k in range(33)] for row in J]
    raw_operators = {}
    for i in PIVOT_ROWS + [92, 134]:
        raw_operators[i] = [[Du * operators[i][a][j] - sum(
            (c * operators[selected[k]][a][j] for k, c in enumerate(projection[i]) if c), zero)
            for j in range(len(SUPPORT))] for a in range(9)]

    relation_data = certificate['two_quadratic_row_relations']
    assert relation_data['rows'] == PIVOT_ROWS
    relations = {int(i): [outer(value) for value in values]
                 for i, values in relation_data['relations'].items()}
    assert set(relations) == {92, 134}
    # Coefficient inputs are checked below as exact polynomial identities.
    syzygy_multipliers = {}
    for i, values in relations.items():
        coefficients = {i: R5.one}
        coefficients.update({j: -value for j, value in zip(PIVOT_ROWS, values) if value})
        denominator = R2.one
        for value in coefficients.values():
            for coefficient in value.values():
                denominator = denominator.lcm(coefficient.denom)
        syzygy_multipliers[i] = {j: clear_outer(value, denominator)
                                for j, value in coefficients.items()}

    # Remove operator contents before clearing the single core denominator.
    primitive_operators, ratios = {}, {}
    for i, matrix in raw_operators.items():
        content = zero
        for row in matrix:
            for value in row:
                content = content.gcd(value)
        assert content and all(not any(m[2:]) for m in content)
        primitive_operators[i] = [[value.exquo(content) for value in row] for row in matrix]
        ratios[i] = K.from_sympy(content.as_expr()) / K.from_sympy(quad['contents'][i].as_expr())
    annihilators = {}
    for i, values in relations.items():
        coefficients = {i: R5(ratios[i])}
        for j, value in zip(PIVOT_ROWS, values):
            coefficient = -value.mul_ground(ratios[j])
            if coefficient:
                coefficients[j] = coefficient
        denominator = R2.one
        for value in coefficients.values():
            for coefficient in value.values():
                denominator = denominator.lcm(coefficient.denom)
        multipliers = {j: clear_outer(value, denominator) for j, value in coefficients.items()}
        numerator = [[sum((c * primitive_operators[j][a][b]
                           for j, c in multipliers.items()), zero)
                      for b in range(len(SUPPORT))] for a in range(9)]
        annihilators[i] = (numerator, denominator)
    print('two tracked cubic operators reconstructed', time.time() - started, flush=True)

    rr_denominator = R2.one
    for row in basis:
        for coefficient in row:
            rr_denominator = rr_denominator.lcm(coefficient.denom)
    normal_terms = lambda row: [(list(m) + normal_exponent(normal), c)
                               for normal, poly in zip(MS2, row) for m, c in poly.items()]
    parameter_terms = lambda poly: [(list(m) + [0] * 9, c) for m, c in poly.items()]
    core_terms = lambda poly: [(list(m) + [0] * 14, c) for m, c in poly.items()]
    with native_input.open('w') as stream:
        stream.write('24 53 2\n')
        write_poly(stream, core_terms(rr_denominator))
        for row in basis:
            terms = []
            for normal, value in zip(MS2, row):
                numerator = value.numer * rr_denominator.exquo(value.denom)
                terms += [(list(m) + [0] * 5 + normal_exponent(normal), c)
                          for m, c in numerator.items()]
            write_poly(stream, terms)
        assert len(quad['nonzero']) == 53
        for i in quad['nonzero']:
            stream.write(str(i) + '\n')
            write_poly(stream, normal_terms(quad['Ag'][i]))
            for pivot in PIVOT_COLS:
                write_poly(stream, parameter_terms(quad['Ag'][i][pivot]))
        for i in [92, 134]:
            multipliers = syzygy_multipliers[i]
            stream.write(str(len(multipliers)) + '\n')
            for j, value in multipliers.items():
                stream.write(str(j) + '\n')
                write_poly(stream, parameter_terms(value))
        stream.write(str(len(SUPPORT)) + '\n')
        for row in lift:
            write_poly(stream, normal_terms(row))
        stream.write('2\n')
        for i, (matrix, denominator) in annihilators.items():
            stream.write(str(i) + '\n')
            total_denominator = denominator * lift_denominator * R2.from_expr(Du.as_expr()) ** 2
            write_poly(stream, core_terms(total_denominator))
            for j in range(len(SUPPORT)):
                terms = []
                for a, row in enumerate(matrix):
                    terms += [(list(m) + normal_exponent((a,)), c) for m, c in row[j].items()]
                write_poly(stream, terms)
    index3 = {monomial: i for i, monomial in enumerate(MS3)}
    with module_input.open('w') as stream:
        stream.write('216 165\n')
        for row in basis:
            for a in range(9):
                multiplied = {index3[tuple(sorted((a,) + monomial))]: value
                              for monomial, value in zip(MS2, row) if value}
                stream.write(str(len(multiplied)) + '\n')
                for j, value in multiplied.items():
                    stream.write(str(j) + '\n')
                    write_core_poly(stream, value.numer)
                    write_core_poly(stream, value.denom)
    return quad, relations


def point_witness(data, quad, relations):
    prime = 32003
    values = dict(zip(s.symbols('e d p q r h k'), [2, 3, 5, 7, 11, 13, 17]))
    def mod(value):
        value = s.Rational(value)
        return int(value.p) % prime * pow(int(value.q), -1, prime) % prime
    convert = lambda value: mod(s.sympify(value).subs(values))
    T = np.array([[convert(v) for j, v in enumerate(row) if j in data['normal']]
                  for row in data['T']], dtype=np.int64)
    U = np.array([[convert(v) for v in row] for row in data['U']], dtype=np.int64)
    J = np.array([[convert(v) for v in row] for row in data['J']], dtype=np.int64)
    projection = J @ U % prime
    model = S61()
    beta = {int(k[1:]) - 1: s.sympify(v).subs(values) for k, v in data['beta'].items()}
    _, equations = specialize(model, beta, True)
    H = np.zeros((135, 45, 45), dtype=np.int64)
    for i, equation in enumerate(equations):
        for monomial, coefficient in equation.items():
            if data['fixed'] in monomial:
                continue
            a, b = [data['cols'].index(k) for k in monomial] + [44] * (2 - len(monomial))
            coefficient = mod(coefficient)
            if a == b:
                H[i, a, b] += coefficient
            else:
                H[i, a, b] += coefficient * pow(2, -1, prime) % prime
                H[i, b, a] += coefficient * pow(2, -1, prime) % prime
    H %= prime
    Q = np.zeros((135, 45), dtype=np.int64)
    for i in range(135):
        matrix = T.T @ (H[i] @ T % prime) % prime
        for k, (a, b) in enumerate(MS2):
            Q[i, k] = matrix[a, b] * (1 if a == b else 2) % prime
    rows = data['pivot_rows']
    W2 = -U @ Q[rows, :] % prime
    A2 = (Q - projection @ Q[rows, :]) % prime
    cubic = np.zeros((135, 165), dtype=np.int64)
    index3 = {monomial: i for i, monomial in enumerate(MS3)}
    for i in range(135):
        matrix = 2 * T.T @ (H[i] @ W2 % prime) % prime
        for a in range(9):
            for k, monomial in enumerate(MS2):
                cubic[i, index3[tuple(sorted((a,) + monomial))]] += matrix[a, k]
    cubic %= prime
    A3 = (cubic - projection @ cubic[rows, :]) % prime
    selected = A2[np.ix_(PIVOT_ROWS, PIVOT_COLS)]
    inverse = np.array([[mod(v) for v in row] for row in
                        s.Matrix(selected.tolist()).inv_mod(prime).tolist()], dtype=np.int64)
    coefficients = A2[:, PIVOT_COLS] @ inverse % prime
    assert not np.any((A2 - coefficients @ A2[PIVOT_ROWS, :]) % prime)
    true_cubic = (A3 - coefficients @ A3[PIVOT_ROWS, :]) % prime
    multiples = []
    for row in A2[PIVOT_ROWS, :]:
        for a in range(9):
            result = np.zeros(165, dtype=np.int64)
            for monomial, value in zip(MS2, row):
                result[index3[tuple(sorted((a,) + monomial))]] = value
            multiples.append(result)
    assert elimination_rank(np.array(multiples), prime)[0] == 129
    assert elimination_rank(np.vstack([multiples, true_cubic[92], true_cubic[134]]), prime)[0] == 131
    # The stored primitive-row syzygies specialize to the unique coefficients
    # independently obtained above from the original 135 equations.
    for i, relation in relations.items():
        gi = convert(quad['contents'][i].as_expr())
        assert gi
        for j, value in zip(PIVOT_ROWS, relation):
            gj = convert(quad['contents'][j].as_expr())
            assert gj
            evaluated = convert(value.as_expr())
            assert gi * pow(gj, -1, prime) * evaluated % prime == int(coefficients[i, PIVOT_ROWS.index(j)])
    return dict(prime=prime, parameters={str(k): str(v) for k, v in values.items()},
                quadratic_linear_rank=129, rank_with_two_tracked_cubics=131,
                original_135_equations_reconstructed=True,
                candidate_syzygy_specializations_checked=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--flint-include', type=Path, default=Path(os.environ.get(
        'FLINT_INCLUDE_DIR', ROOT / 'work/flint-headers/flint/3.5.0/include')))
    parser.add_argument('--mpfr-include', type=Path, default=Path(os.environ.get(
        'MPFR_INCLUDE_DIR', ROOT / 'work/flint-headers/mpfr/4.2.2/include')))
    parser.add_argument('--flint-library', type=Path, default=Path(os.environ.get(
        'FLINT_LIBRARY_FILE', ROOT / 'work/m2-runtime/lib/libflint.24.0.dylib')))
    args = parser.parse_args()
    assert not args.output.exists(), 'output must be new'
    started = time.time()
    certificate = json.loads(CERTIFICATE.read_text())
    for path, digest in certificate['input_sha256'].items():
        assert sha(ROOT / path) == digest, path
    data = json.loads(LINEAR.read_text())['linear']
    assert (args.flint_include / 'flint/fmpq_mpoly.h').is_file()
    assert args.flint_library.is_file()
    compiler = shutil.which('clang++') or shutil.which('g++')
    assert compiler, 'a C++17 compiler and FLINT headers/library are required'
    environment = os.environ.copy()
    environment['DYLD_LIBRARY_PATH'] = ':'.join([
        str(args.flint_library.parent), '/opt/homebrew/lib',
        environment.get('DYLD_LIBRARY_PATH', '')])
    sources = ['code/s61_binary_deeper_cubic_products_20261101.cpp',
               'code/s61_binary_deeper_cubic_field_20261101.cpp']
    with tempfile.TemporaryDirectory(prefix='poisson-deeper-cubic-', dir='/tmp') as directory:
        directory = Path(directory)
        binaries = []
        for source in sources:
            binary = directory / Path(source).stem
            command = [compiler, '-O3', '-std=c++17', '-I' + str(args.flint_include),
                       '-I' + str(args.mpfr_include), '-I/opt/homebrew/include',
                       str(ROOT / source), str(args.flint_library), '-L/opt/homebrew/lib',
                       '-lgmp', '-o', str(binary)]
            subprocess.run(command, check=True, timeout=60)
            binaries.append(binary)
        native_input, module_input = directory / 'polynomials.txt', directory / 'module.txt'
        quad, relations = prepare(data, certificate, native_input, module_input)
        sample = point_witness(data, quad, relations)
        cubic_output, field_report = directory / 'cubics.txt', directory / 'field.json'
        subprocess.run([str(binaries[0]), str(native_input), str(cubic_output)],
                       check=True, env=environment, timeout=600)
        subprocess.run([str(binaries[1]), str(module_input), str(cubic_output), str(field_report)],
                       check=True, env=environment, timeout=1200)
        result = json.loads(field_report.read_text())
        assert result['all_passed']
    dependencies = [str(Path(__file__).relative_to(ROOT)), *sources,
                    'code/s61_binary_deeper_cubic_field_20261101.hpp',
                    str(CERTIFICATE.relative_to(ROOT)), str(LINEAR.relative_to(ROOT)),
                    str(QUADRATIC.relative_to(ROOT)), 'code/s61_binary_deeper_uniform_quadratic_20261101.py',
                    'code/s61_model.py', 'code/s61_flag_transverse_obstruction_20261101.py']
    report = dict(id='S61_BINARY_DEEPER_CUBIC_NUISANCE_20261101', all_passed=True,
                  coefficient_field='QQ(e,d)(p,q,r,h,k)', normal_variable_count=9,
                  field_checks=result, independent_original_equation_point_check=sample,
                  claim='Two genuine cubic initial classes descend to QQ(e,d); their generic quotient classes are independent.',
                  limitations=['No full cubic initial-ideal exhaustion is asserted.',
                               'No uniform four-jet bound or whole deeper-stratum removal is asserted.',
                               'The identities are over the stated rational coefficient field.'],
                  proof_note='research/notes/S61_BINARY_DEEPER_CUBIC_NUISANCE_20261101.md',
                  input_sha256={path: sha(ROOT / path) for path in dependencies},
                  flint_library_sha256=sha(args.flint_library), seconds=time.time() - started)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
    print('PASS', report['seconds'], args.output, flush=True)


if __name__ == '__main__':
    main()
