"""Reconstruct the actual six-normal-variable ZERO12 jets from all raw135 equations."""
import json
import os
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as s

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT = Path(__file__).resolve().parents[1]


def normal_jets():
    start = time.time()
    low = json.loads((ROOT / 'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json').read_text())
    old = json.loads((ROOT / 'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text())
    model = S61()
    b = s.symbols('b0:15')
    parameters = list(map(s.Symbol, low['parameters']))
    normal = [2, 4, 5, 6, 7, 8]
    np, n, maximum_degree = len(parameters), len(normal), 6
    sub = {s.Symbol(k): s.sympify(v) for k, v in low['substitutions'].items()}
    H = s.Matrix(old['flag']['beta_kernel']).applyfunc(s.Rational)
    beta = (H * s.Matrix(b)).subs(sub)
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    q0 = s.Matrix(old['flag']['point61']).applyfunc(s.Rational)[two, :]
    fixed = two.index(37)
    columns = [i for i in range(45) if i != fixed]
    _, equations = specialize(model, dict(zip(one, beta)), True)
    J = s.zeros(135, 45)
    for row, f in enumerate(equations):
        for mon, coefficient in f.items():
            if len(mon) == 2:
                i, j = mon
                if i != fixed:
                    J[row, columns.index(i)] += coefficient * q0[j]
                if j != fixed:
                    J[row, columns.index(j)] += coefficient * q0[i]
            elif len(mon) == 1:
                J[row, 44] += coefficient * q0[mon[0]]
    rows, pivots = low['pivot_rows'], low['pivot_columns']
    block = J[rows, pivots].to_DM().to_field()
    assert s.cancel(block.domain.to_sympy(block.det())) == -s.Rational(1400832, 25)
    inverse = block.inv().to_Matrix()
    U = s.zeros(45, 36)
    for i, column in enumerate(pivots):
        U[column, :] = inverse[i, :]
    T = s.Matrix(low['tangent']).applyfunc(s.sympify)[:, normal]
    assert (J * T).applyfunc(s.expand) == s.zeros(135, n)
    assert (J[rows, :] * U).applyfunc(s.expand) == s.eye(36)

    def write_poly(stream, expression):
        terms = [(e, c) for e, c in s.Poly(expression, *parameters).terms() if c]
        stream.write(str(len(terms)) + '\n')
        for exponent, coefficient in terms:
            stream.write(str(coefficient) + ' ' + ' '.join(map(str, tuple(exponent) + (0,) * n)) + '\n')

    env = os.environ.copy()
    env['DYLD_LIBRARY_PATH'] = str(ROOT / 'work/m2-runtime/lib') + ':/opt/homebrew/lib'
    with tempfile.TemporaryDirectory(prefix='zero12-normal-six-') as directory:
        directory = Path(directory)
        input_path, output_path, binary = [directory / name for name in ('input', 'output', 'engine')]
        compiled = subprocess.run([
            'c++', '-O3', '-std=c++17',
            '-I' + str(ROOT / 'work/flint-headers/flint/3.5.0/include'),
            '-I/opt/homebrew/include', str(ROOT / 'code/parametric_implicit_jets_20261101.cpp'),
            '-L' + str(ROOT / 'work/m2-runtime/lib'), '-L/opt/homebrew/lib',
            '-lflint', '-lgmp', '-o', str(binary)], capture_output=True, text=True, timeout=60)
        assert compiled.returncode == 0, compiled.stderr
        with input_path.open('w') as stream:
            stream.write(f'135 {np} {n} {maximum_degree} 36\n')
            stream.write(' '.join(map(str, rows)) + '\n')
            write_poly(stream, 1)
            for matrix in (J, T, U):
                for row in matrix.tolist():
                    for value in row:
                        write_poly(stream, value)
            for f in equations:
                terms = []
                for mon, coefficient in f.items():
                    if fixed in mon:
                        continue
                    i, j = [columns.index(k) for k in mon] + [44] * (2 - len(mon))
                    terms.append((i, j, coefficient))
                stream.write(str(len(terms)) + '\n')
                for i, j, coefficient in terms:
                    stream.write(f'{i} {j}\n')
                    write_poly(stream, coefficient)
        native = subprocess.run([str(binary), str(input_path), str(output_path)],
                                env=env, capture_output=True, text=True, timeout=300)
        assert native.returncode == 0, native.stdout + native.stderr
        matrices = {}
        with output_path.open() as stream:
            assert list(map(int, next(stream).split())) == [135, np, n, maximum_degree]
            for degree in range(2, maximum_degree + 1):
                mons = list(combinations_with_replacement(range(n), degree))
                index = {mon: i for i, mon in enumerate(mons)}
                for key, count in [('W', 45), ('A', 135)]:
                    matrix = []
                    for _ in range(count):
                        row = [[] for _ in mons]
                        for _ in range(int(next(stream))):
                            fields = next(stream).split()
                            coefficient, exponent = fields[0], list(map(int, fields[1:]))
                            assert len(exponent) == np + n
                            mon = tuple(j for j, e in enumerate(exponent[np:]) for _ in range(e))
                            assert len(mon) == degree
                            row[index[mon]].append([exponent[:np], coefficient])
                        matrix.append(row)
                    matrices[key + str(degree)] = matrix
            assert not stream.read().strip()

    # Check the independent previously stored full-nine-variable QQ jets by restriction.
    ring = s.QQ.poly_ring(*parameters)
    PR = ring.ring
    for degree in (2, 3):
        mons = list(combinations_with_replacement(range(n), degree))
        for key, saved in [('W', low['lift']), ('A', low['residual'])]:
            for i, row in enumerate(saved[str(degree)]):
                expected = {tuple(normal.index(j) for j in mon): PR.from_expr(s.sympify(c))
                            for mon, c in row if all(j in normal for j in mon)}
                actual = {mon: PR.from_dict({tuple(e): s.QQ(c) for e, c in polynomial})
                          for mon, polynomial in zip(mons, matrices[key + str(degree)][i]) if polynomial}
                expected = {mon: c for mon, c in expected.items() if c}
                assert actual == expected, (key, degree, i)
    for degree in range(2, maximum_degree + 1):
        assert all(not any(matrices['A' + str(degree)][i]) for i in rows)
        assert all(not any(matrices['W' + str(degree)][i]) for i in low['free_coordinates'])
    matrices.update(parameters=list(map(str, parameters)), normal_columns=normal,
                    parameter_columns=[0, 1, 3], inverse_minor='-1400832/25',
                    native_log=native.stdout + native.stderr,
                    independent_low_jet_restriction_verified=True,
                    seconds=time.time() - start)
    return matrices
