"""Native exact polynomial operations used by the sixth-order certificate."""
import os
from itertools import combinations_with_replacement
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]

class Backend:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.environment = os.environ.copy()
        self.environment['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        self.binaries = {}
        self.timings = []
        for kind, source in [('contract', 'multivariate_jet_contraction_20261101.cpp'),
                             ('matrix', 'multivariate_polynomial_matrix_20261101.cpp')]:
            binary = self.directory/kind
            command = ['c++', '-O3', '-std=c++17',
                       '-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
                       '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),
                       '-I/opt/homebrew/include', str(ROOT/'code'/source),
                       str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),
                       '-L/opt/homebrew/lib', '-lgmp', '-o', str(binary)]
            subprocess.run(command, check=True, capture_output=True, timeout=60)
            self.binaries[kind] = binary

    @staticmethod
    def polynomial(stream, value):
        stream.write(str(len(value))+'\n')
        for exponent, coefficient in value.items():
            stream.write(str(coefficient)+' '+' '.join(map(str, exponent))+'\n')

    def matrix(self, stream, matrix):
        for row in matrix:
            entries = [(j, value) for j, value in enumerate(row) if value]
            stream.write(str(len(entries))+'\n')
            for j, value in entries:
                stream.write(str(j)+'\n')
                self.polynomial(stream, value)

    @staticmethod
    def read_matrix(path, ring):
        with Path(path).open() as stream:
            nr, nc, n = map(int, stream.readline().split())
            assert n == 5
            out = []
            for _ in range(nr):
                row = [ring.zero]*nc
                for _ in range(int(stream.readline())):
                    column = int(stream.readline())
                    terms = {}
                    for _ in range(int(stream.readline())):
                        line = stream.readline().split()
                        terms[tuple(map(int, line[1:]))] = s.QQ.from_sympy(s.Rational(line[0]))
                    row[column] = ring.from_dict(terms)
                out.append(row)
            assert not stream.read().strip()
        return out

    def execute(self, kind, source, output):
        started = time.time()
        result = subprocess.run([str(self.binaries[kind]), str(source), str(output)],
                                env=self.environment, text=True, capture_output=True, timeout=300)
        if result.returncode:
            raise RuntimeError(f'{kind} failed: {result.stderr[-4000:]}')
        self.timings.append(dict(kind=kind, seconds=time.time()-started, final_log=result.stderr.splitlines()[-1]))
        assert output.is_file() and output.stat().st_size

    def linear(self, ring, left, coefficients, right, scale=None):
        if scale is None:
            scale = ring.one
        nr, nc, nk = len(left), len(left[0]), len(right)
        assert len(coefficients) == nr and all(len(row) == nk for row in coefficients)
        assert all(len(row) == nc for row in left+right)
        with tempfile.TemporaryDirectory(prefix='matrix-', dir=self.directory) as temp:
            source, output = [Path(temp)/name for name in ['input.txt', 'output.txt']]
            with source.open('w') as stream:
                stream.write(f'5 {nr} {nk} {nc}\n')
                for matrix in [left, coefficients, right]:
                    self.matrix(stream, matrix)
                self.polynomial(stream, scale)
            self.execute('matrix', source, output)
            answer = self.read_matrix(output, ring)
        assert len(answer) == nr and len(answer[0]) == nc
        return answer

    def contract(self, ring, matrices, monomials, degree, normal_forms, output_columns, pairs, coefficients):
        """Contract sum H(W_i,W_j); the supplied coefficient matrix acts on equation rows."""
        coordinates = 45
        eq = len(coefficients)
        assert len(matrices) == degree-1 and all(len(matrix) == coordinates for matrix in matrices)
        full = list(combinations_with_replacement(range(11), degree))
        index = {mon: j for j, mon in enumerate(full)}
        assert len(full) == len(normal_forms)
        active = [j for j in range(len(pairs)) if any(row[j] for row in coefficients)]
        with tempfile.TemporaryDirectory(prefix='contract-', dir=self.directory) as temp:
            source, output = [Path(temp)/name for name in ['input.txt', 'output.txt']]
            with source.open('w') as stream:
                stream.write(f'5 {coordinates} {eq} {len(full)} {output_columns} {degree}\n')
                for matrix in matrices:
                    stream.write(str(len(matrix[0]))+'\n')
                    for row in matrix:
                        for value in row:
                            self.polynomial(stream, value)
                for da in range(degree//2):
                    db = degree-2-da
                    for left in monomials[da]:
                        stream.write(' '.join(str(index[tuple(sorted(left+right))]) for right in monomials[db])+'\n')
                for row in normal_forms:
                    stream.write(str(len(row))+'\n')
                    for column, value in row:
                        stream.write(str(column)+'\n')
                        self.polynomial(stream, value)
                stream.write(str(len(active))+'\n')
                for column in active:
                    i, j = pairs[column]
                    rows = [k for k in range(eq) if coefficients[k][column]]
                    stream.write(f'{i} {j} {len(rows)}\n')
                    for row in rows:
                        stream.write(str(row)+'\n')
                        self.polynomial(stream, coefficients[row][column])
            self.execute('contract', source, output)
            answer = self.read_matrix(output, ring)
        assert len(answer) == eq and len(answer[0]) == output_columns
        return answer
