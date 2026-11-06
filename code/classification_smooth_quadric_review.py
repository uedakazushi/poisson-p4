#!/usr/bin/env python3
"""Independent Fraction/polynomial review of the regular smooth-quadric source."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
ZERO = (0,)*10
PAIRS = list(itertools.combinations(range(5), 2))
PARAMETERS = ['p12', 'p14', 'p15', 'p16', 'p18']


def scalar(c):
    return {ZERO: F(c)} if c else {}


def variable(j):
    return {tuple(int(i == j) for i in range(10)): F(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for mon, c in poly.items():
            out[mon] = out.get(mon, F(0))+c
    return {m: c for m, c in out.items() if c}


def scale(poly, c):
    return {m: c*v for m, v in poly.items() if c*v}


def multiply(a, b):
    out = {}
    for m, c in a.items():
        for n, d in b.items():
            word = tuple(x+y for x, y in zip(m, n))
            out[word] = out.get(word, F(0))+c*d
    return {m: c for m, c in out.items() if c}


def power(poly, n):
    out = scalar(1)
    for _ in range(n):
        out = multiply(out, poly)
    return out


def derivative(poly, j):
    return {tuple(v-int(i == j) for i, v in enumerate(m)): c*m[j]
            for m, c in poly.items() if m[j]}


def parse(text, extra=None):
    env = {f'x{i}': variable(i) for i in range(5)}
    env.update({name: variable(i+5) for i, name in enumerate(PARAMETERS)})
    env.update(extra or {})
    def read(node):
        if isinstance(node, ast.Constant):
            return scalar(node.value)
        if isinstance(node, ast.Name):
            return env[node.id]
        if isinstance(node, ast.UnaryOp):
            assert isinstance(node.op, (ast.USub, ast.UAdd))
            return scale(read(node.operand), -1 if isinstance(node.op, ast.USub) else 1)
        assert isinstance(node, ast.BinOp)
        if isinstance(node.op, ast.Pow):
            return power(read(node.left), ast.literal_eval(node.right))
        a, b = read(node.left), read(node.right)
        if isinstance(node.op, ast.Add): return add(a, b)
        if isinstance(node.op, ast.Sub): return add(a, scale(b, -1))
        if isinstance(node.op, ast.Mult): return multiply(a, b)
        assert isinstance(node.op, ast.Div) and set(b) == {ZERO}
        return scale(a, 1/b[ZERO])
    return read(ast.parse(str(text), mode='eval').body)


def entry(P, i, j):
    return P.get((i, j), {}) if i < j else scale(P.get((j, i), {}), -1) if i > j else {}


def jacobi(P):
    out = {}
    for i, j, k in itertools.combinations(range(5), 3):
        terms = []
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            terms += [multiply(entry(P, a, r), derivative(entry(P, b, c), r)) for r in range(5)]
        out[(i, j, k)] = add(*terms)
    return out


def rank(rows):
    rows = [list(map(F, row)) for row in rows]
    if not rows: return 0
    at = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(at, len(rows)) if rows[i][column]), None)
        if pivot is None: continue
        rows[at], rows[pivot] = rows[pivot], rows[at]
        rows[at] = [v/rows[at][column] for v in rows[at]]
        for i in range(at+1, len(rows)):
            if rows[i][column]:
                c = rows[i][column]
                rows[i] = [a-c*b for a, b in zip(rows[i], rows[at])]
        at += 1
        if at == len(rows): break
    return at


def substitute(poly, values):
    out = {}
    for mon, c in poly.items():
        term = scalar(c)
        for i, exponent in enumerate(mon):
            term = multiply(term, power(values.get(i, variable(i)), exponent))
        out = add(out, term)
    return out


def pfaffians(P):
    out = []
    for missing in range(5):
        a, b, c, d = [i for i in range(5) if i != missing]
        out.append(scale(add(multiply(entry(P,a,b), entry(P,c,d)),
                            scale(multiply(entry(P,a,c),entry(P,b,d)), -1),
                            multiply(entry(P,a,d),entry(P,b,c))), (-1)**missing))
    return out


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (30, 31))
    resource.setrlimit(resource.RLIMIT_AS, (1024*1024**2, 1024*1024**2))
    begun = time.monotonic()
    path = ROOT/'reports/CLASSIFICATION_SMOOTH_QUADRIC_REGULAR_COVERAGE_20261107.json'
    source = json.loads(path.read_text())
    assert hashlib.sha256(source['producer_source'].encode()).hexdigest() == source['producer_sha256']
    weights = (0, 1, -1, 2, -2)
    q = parse('x0**2+x1*x2+x3*x4')
    xmons = {n: sorted((e for e in itertools.product(range(n+1), repeat=5) if sum(e) == n), reverse=True)
             for n in (1, 2, 3)}
    invariant = [(ij, mon) for ij in PAIRS for mon in xmons[2]
                 if sum(weights[k]*mon[k] for k in range(5)) == weights[ij[0]]+weights[ij[1]]]
    assert len(invariant) == 20
    def linear_map(P):
        div = [add(*(derivative(entry(P, i, j), i) for i in range(5))) for j in range(5)]
        contraction = [add(*(multiply(derivative(q, i), entry(P, i, j)) for i in range(5))) for j in range(5)]
        return div+contraction
    rows = [(j, mon+(0,)*5) for j in range(10) for mon in xmons[1 if j < 5 else 3]]
    columns = []
    for ij, mon in invariant:
        mapped = linear_map({ij: {mon+(0,)*5: F(1)}})
        columns.append([mapped[j].get(e, F(0)) for j, e in rows])
    matrix = [list(row) for row in zip(*columns)]
    rhs_polys = [{}]*5+[scale(multiply(q, variable(j)), weights[j]) for j in range(5)]
    rhs = [rhs_polys[j].get(e, F(0)) for j, e in rows]
    assert rank(matrix) == rank([row+[value] for row, value in zip(matrix, rhs)]) == 15
    P = {(int(key[0]), int(key[1])): parse(value) for key, value in source['brackets'].items()}
    assert all((ij, mon[:5]) in invariant for ij, poly in P.items() for mon in poly)
    mapped = linear_map(P)
    assert all(not add(a, scale(b, -1)) for a, b in zip(mapped, rhs_polys))
    directions = []
    for parameter in range(5, 10):
        direction = {ij: substitute(derivative(poly, parameter), {i: {} for i in range(5, 10)}) for ij, poly in P.items()}
        directions.append([direction[ij].get(mon+(0,)*5, F(0)) for ij in PAIRS for mon in xmons[2]])
    assert rank(directions) == 5
    raw = jacobi(P)
    coefficients = []
    for poly in raw.values():
        for xmon in xmons[3]:
            coeff = {mon[5:]: c for mon, c in poly.items() if mon[:5] == xmon}
            coefficients.append(coeff)
    conditions = [parse(text) for text in ['p12*p14', 'p12*(p15+1)',
        'p14*(p16+p18+3)', 'p18*(p15+p16-1)', '8*p18*(p15+p16)+30*(p15+p16)+7*p18+20']]
    conditions5 = [{mon[5:]: c for mon, c in poly.items()} for poly in conditions]
    parameter_mons = sorted(set().union(*(set(poly) for poly in coefficients+conditions5)))
    def coefficient_matrix(polys):
        return [[poly.get(mon, F(0)) for mon in parameter_mons] for poly in polys]
    assert rank(coefficient_matrix(coefficients)) == rank(coefficient_matrix(conditions5)) == rank(coefficient_matrix(coefficients+conditions5)) == 5
    branch_checks = []
    for branch in source['six_complete_sources']:
        values = {PARAMETERS.index(key)+5: parse(value, {'z': variable(5)}) for key, value in branch['substitution'].items()}
        assert all(not substitute(poly, values) for poly in conditions)
        branchP = {ij: substitute(poly, values) for ij, poly in P.items()}
        assert all(not poly for poly in jacobi(branchP).values())
        pf = pfaffians(branchP)
        D = scale(add(*(multiply(variable(i), pf[i]) for i in range(5))), F(1, 5))
        assert D == parse(branch['actual_divisor'], {'z': variable(5)})
        assert all(derivative(D, i) == pf[i] for i in range(5))
        degeneration = branch['diagonal_degeneration_weights']
        limit = {}
        powers = set()
        for ij, poly in branchP.items():
            value = {}
            for mon, c in poly.items():
                degree = degeneration[ij[0]]+degeneration[ij[1]]-sum(mon[k]*degeneration[k] for k in range(5))
                assert degree >= 0
                powers.add(degree)
                if degree == 0:
                    assert mon[:5] == tuple(int(k in ij) for k in range(5))
                    value[mon] = c
            limit[ij] = value
        limit_pf = pfaffians(limit)
        limit_D = scale(add(*(multiply(variable(i), limit_pf[i]) for i in range(5))), F(1, 5))
        assert limit_D == parse('10*x0*x1*x2*x3*x4/3')
        branch_checks.append(dict(name=branch['name'], all_parameter_Jacobi=True,
            actual_divisor_and_all_five_gradients=True, nonnegative_powers=sorted(powers),
            full_parameter_SNC_limit=True))
    report = dict(task_id='E_reduced_nonnormal', all_passed=True,
        independent_Fraction_polynomials=True, producer_and_SymPy_not_imported=True,
        invariant_quadratic_count=20, full_linear_coefficient_matrix_shape=[len(matrix), 20],
        full_linear_rank=15, augmented_rank=15, displayed_source_direction_rank=5,
        displayed_five_parameter_source_exhausts_linear_solutions=True,
        all_ten_strict_Jacobi_polynomials_reconstructed=True,
        Jacobi_coefficient_slots=len(coefficients), exact_five_condition_rowspan_rank=5,
        six_line_support_hand_proof='t(s-1)=0 and (8t+30)s+7t+20=0 yield (t,s)=(0,-2/3) or (-10/3,1); ab,a(u+1),b(v+t+3) split each into precisely the three recorded lines, including all intersections.',
        branch_checks=branch_checks,
        global_MPS_support_implication='Every source point has a regular PGL orbit degeneration to a nonzero SNC divisor; every component through it is PGL-invariant and therefore among the MPS40.',
        unhandled_Jordan_and_singular_quadric_cases_retained=True,
        input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest(),
                      'code/classification_smooth_quadric_review.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        producer_source=Path(__file__).read_text(), seconds=time.monotonic()-begun)
    with output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps({key: report[key] for key in ('all_passed', 'full_linear_rank', 'exact_five_condition_rowspan_rank', 'seconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)
