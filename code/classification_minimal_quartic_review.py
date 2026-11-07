#!/usr/bin/env python3
"""Independent sparse verification of the complete minimal-nilpotent source.

No producer or CAS imports. The reusable cofactor reviewer's sparse arithmetic
is extended only to Laurent monomials in the explicitly inverted parameter c.
"""
import argparse
import ast
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

from classification_quadratic_cofactor_review import Ring

ROOT = Path(__file__).resolve().parents[1]
PAIRS = tuple(itertools.combinations(range(4), 2))


class LaurentRing(Ring):
    def __init__(self, parameters, coordinates, invert=()):
        self.names = tuple(parameters)+tuple(coordinates)
        self.n = len(parameters)
        self.zero = (0,)*len(self.names)
        self.invert = set(invert)

    def inverse_monomial(self, polynomial):
        assert len(polynomial) == 1
        powers, coefficient = next(iter(polynomial.items()))
        assert all(not a or name in self.invert for name, a in zip(self.names, powers))
        return {tuple(-a for a in powers):1/coefficient}

    def parse(self, expression):
        def read(node):
            if isinstance(node, ast.Constant) and type(node.value) is int:
                return self.constant(node.value)
            if isinstance(node, ast.Name):
                return self.variable(node.id)
            if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
                return self.scale(read(node.operand), -1)
            if isinstance(node, ast.BinOp):
                a, b = read(node.left), read(node.right)
                if isinstance(node.op, ast.Add):
                    return self.add(a, b)
                if isinstance(node.op, ast.Sub):
                    return self.add(a, self.scale(b, -1))
                if isinstance(node.op, ast.Mult):
                    return self.mul(a, b)
                if isinstance(node.op, ast.Div):
                    return self.mul(a, self.inverse_monomial(b))
                if isinstance(node.op, ast.Pow):
                    assert set(b) <= {self.zero}
                    exponent = b.get(self.zero, Q(0))
                    assert exponent.denominator == 1
                    if exponent < 0:
                        a = self.inverse_monomial(a)
                    result = self.constant(1)
                    for _ in range(abs(int(exponent))):
                        result = self.mul(result, a)
                    return result
            raise ValueError(ast.dump(node))
        return read(ast.parse(expression, mode='eval').body)


def determinant(ring, matrix):
    """Sparse Laplace expansion, choosing the least populated remaining row."""
    n = len(matrix)

    @lru_cache(None)
    def minor(rows, columns):
        if not rows:
            return ring.constant(1)
        row_position = min(range(len(rows)), key=lambda p:sum(bool(matrix[rows[p]][j]) for j in columns))
        i = rows[row_position]
        reduced_rows = rows[:row_position]+rows[row_position+1:]
        result = {}
        for position, j in enumerate(columns):
            if matrix[i][j]:
                sign = (-1)**(row_position+position)
                term = ring.mul(matrix[i][j], minor(reduced_rows, columns[:position]+columns[position+1:]))
                result = ring.add(result, ring.scale(term, sign))
        return result

    return minor(tuple(range(n)), tuple(range(n)))


def check_pde(report):
    ring = LaurentRing(('c','u','v'), ('y','t','unused1','unused2'))
    monomials = [(i,j) for i in range(5) for j in range(5-i)]
    # Direct action of (c+uy+vy²+3yt)d_t-y²d_y on y^i t^j.
    columns = []
    for i, j in monomials:
        column = {}
        if j:
            column[(i,j-1)] = ring.parse(f'{j}*c')
            column[(i+1,j-1)] = ring.parse(f'{j}*u')
            column[(i+2,j-1)] = ring.parse(f'{j}*v')
        if 3*j-i:
            column[(i+1,j)] = ring.constant(3*j-i)
        columns.append(column)
    labels = sorted({powers for column in columns for powers in column}, reverse=True)
    matrix = [[column.get(powers,{}) for column in columns] for powers in labels]
    assert matrix == [[ring.parse(value) for value in row] for row in report['full_matrix']]
    actual_monomials = [ring.parse(f'y**{i}*t**{j}') for i,j in monomials]
    assert actual_monomials == [ring.parse(value) for value in report['all_degree_at_most_four_monomials']]
    I = ring.parse('y**3*t+c*y**2/2+u*y**3/3+v*y**4/4')
    kernel = [ring.constant(1), I]
    for polynomial in kernel:
        vector = [{e[:-4]+(0,)*4:value for e,value in polynomial.items() if e[-4:-2]==(i,j)}
                  for i,j in monomials]
        assert all(not ring.add(*(ring.mul(a,b) for a,b in zip(row,vector))) for row in matrix)
    assert [ring.parse(value) for value in report['complete_kernel_basis']] == kernel
    rows, cols = report['pivot_rows'], report['pivot_columns']
    assert len(rows) == len(cols) == 13 and len(set(rows)) == len(set(cols)) == 13
    value = determinant(ring, [[matrix[i][j] for j in cols] for i in rows])
    assert value == ring.constant(14929920) == ring.parse(report['pivot_minor'])
    # The two kernel vectors are independent for every c,u,v: the y³t
    # coefficient of I is1, and the constant coefficient of 1 is1.
    assert I.get((0,0,0,3,1,0,0)) == 1
    return dict(full_matrix_shape=[len(matrix),15], every_matrix_coefficient_checked=True,
                independent_kernel=['1','y³t+cy²/2+uy³/3+vy⁴/4'], rank=13,
                determinant=str(value[ring.zero]), determinant_is_nonzero_constant=True,
                all_parameter_values_covered=True, output_monomial_labels=[list(p) for p in labels])


def evaluate_parameters(ring, polynomial, point):
    result = {}
    for e, coefficient in polynomial.items():
        for name, power in zip(ring.names[:ring.n], e[:ring.n]):
            coefficient *= point[name]**power
        coordinate_powers = e[ring.n:]
        result[coordinate_powers] = result.get(coordinate_powers,Q(0))+coefficient
    return {e:c for e,c in result.items() if c}


def verify_source(report):
    parameters = ('a0','c','u','v','A')+tuple(f'b{i}' for i in range(6))
    ring = LaurentRing(parameters, ('x','y','z','w'), invert=('c',))
    f = ring.parse('c+u*y+v*y**2+3*y*w')
    I = ring.parse('y**3*w+c*y**2/2+u*y**3/3+v*y**4/4')
    F = ring.add(ring.parse('a0*c'),ring.scale(I,-1))
    a = ring.parse('a0+y*(-a0*u/c-3*a0*w/c)+A*y**2')
    numerator = ring.add(F,ring.scale(ring.mul(a,f),-1))
    d = {}
    for e, coefficient in numerator.items():
        assert e[ring.n+1] >= 2
        exponents = list(e)
        exponents[ring.n+1] -= 2
        d[tuple(exponents)] = coefficient
    B = ring.parse('b0+b1*y+b2*w+b3*y**2+b4*y*w+b5*w**2')
    b = ring.add(ring.parse('x*y'), B,
                 ring.mul(ring.variable('z'),ring.add(ring.parse('-y'),
                            ring.scale(ring.derivative(a,1),-1),ring.scale(ring.derivative(d,3),-1))))
    entries = [a,b,d,ring.parse('y**2'),{},f]
    assert entries == [ring.parse(value) for value in report['brackets']]
    assert F == ring.parse(report['actual_affine_pfaffian'])
    assert all(sum(e[-4:]) <= 2 for value in entries for e in value)
    P = [[{} for _ in range(4)] for _ in range(4)]
    for (i,j),value in zip(PAIRS,entries):
        P[i][j],P[j][i] = value,ring.scale(value,-1)
    divergence = [ring.add(*(ring.derivative(P[i][j],i) for i in range(4))) for j in range(4)]
    assert divergence == [ring.parse('y'),{},{},{}]
    jacobi = [ring.add(*(ring.mul(P[a][k],ring.derivative(P[b][d],k))
                        for a,b,d in [(i,j,l),(j,l,i),(l,i,j)] for k in range(4)))
              for i,j,l in itertools.combinations(range(4),3)]
    assert jacobi == [{},{},{},{}]
    pf = ring.add(ring.mul(a,f),ring.mul(d,ring.parse('y**2')))
    assert pf == F
    hamiltonian = [ring.add(*(ring.mul(P[i][j],ring.derivative(F,j)) for j in range(4))) for i in range(4)]
    assert hamiltonian == [ring.scale(ring.mul(ring.variable('y'),F),-1),{},{},{}]
    assert ring.mul(f,ring.derivative(I,3)) == ring.mul(ring.parse('y**2'),ring.derivative(I,1))
    assert ring.derivative(I,1) == ring.mul(ring.variable('y'),f)
    assert ring.derivative(I,3) == ring.parse('y**3')

    # Generate the old MINIMAL_M tensor from its independent legacy formula.
    # X=y,Y=6x,Z=w,W=-2z; new coordinate rows are (Y/6,X,-W/2,Z).
    old_ring = LaurentRing((),('X','Y','Z','W'))
    new_ring = LaurentRing((),('x','y','z','w'))
    old_entries = [old_ring.parse(value) for value in
                   ['1+3*X**2-6*X*Z','0','-2*X**2','3*X*Z-18*Z**2',
                    '6*X**2-2*X*Y-3*X*W+30*Z*W','1+6*X*Z']]
    substitutions = [new_ring.parse(value) for value in ['y','6*x','w','-2*z']]
    def substitute(polynomial):
        result = {}
        for powers, coefficient in polynomial.items():
            term = new_ring.constant(coefficient)
            for value,power in zip(substitutions,powers):
                for _ in range(power):
                    term = new_ring.mul(term,value)
            result = new_ring.add(result,term)
        return result
    old_P = [[{} for _ in range(4)] for _ in range(4)]
    for (i,j),value in zip(PAIRS,old_entries):
        old_P[i][j],old_P[j][i] = value,old_ring.scale(value,-1)
    order,scale = [1,0,3,2],[Q(1,6),Q(1),Q(-1,2),Q(1)]
    transformed = [new_ring.scale(substitute(old_P[order[i]][order[j]]),scale[i]*scale[j])
                   for i,j in PAIRS]
    point = dict(a0=Q(-1,6),c=Q(1,2),u=Q(0),v=Q(0),A=Q(-1,2),
                 b0=Q(0),b1=Q(0),b2=Q(0),b3=Q(-1,2),b4=Q(0),b5=Q(0))
    assert transformed == [evaluate_parameters(ring,value,point) for value in entries]
    assert point == {key:Q(value) for key,value in report['MINIMAL_M_source_point'].items()}
    return dict(full_parameter_count=11, every_source_bracket_reconstructed=True,
                polynomial_division_by_y_squared_exact=True, whole_divergence=["y","0","0","0"],
                all_four_strict_Jacobi_zero=True, all_source_Pfaffian_checked=True,
                complete_modular_identity_checked=True, I_derivative_identities_checked=True,
                all_coordinate_degrees_at_most_two=True,
                old_MINIMAL_M_from_legacy_formula=True, old_M_all_six_brackets_exactly_matched=True,
                old_M_parameter_point={name:str(value) for name,value in point.items()},
                only_display_denominator_parameter='c')


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    started = time.monotonic()
    geometry_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_RANK1_GEOMETRY_20261107.json'
    support_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_RANK1_SUPPORT_20261107.json'
    geometry = json.loads(geometry_path.read_text())
    support = json.loads(support_path.read_text())
    for report in [geometry,support]:
        assert hashlib.sha256(report['source_bytes'].encode()).hexdigest() == report['source_sha256']
    pde = check_pde(geometry['nil_rank_one_divisor_geometry'])
    source = verify_source(support['nil_rank_one_support'])
    dependencies = [Path(__file__),Path(__file__).with_name('classification_quadratic_cofactor_review.py'),
                    geometry_path,support_path,ROOT/'legacy/v07/minimal_thickness.tex']
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},
                  full_PDE=pde,whole_source=source,
                  independent_mathematical_review='IrreducibleF implies nonzero y=0 restriction. Modular z/w equations force a0*f0 constant nonzero, hence both are nonzero constants and l,m are multiples ofy. SL2 covers every nonzero(alpha,beta); zero case yields a binary reducible quartic. Complete degree4 Casimir solutions and modular x-equation force k=-1. Coprimality(f,y) and degreebounds recover all source coefficients. The eleven-parameter source is irreducible, contains the actual MINIMAL_M point, and the known unique ambient component at that point therefore contains the entire source.',
                  inverted_elements=['c in the displayed source','a0!=0 is forced by irreducibleF restriction','alpha or beta only for the SL2 cofactor normalization'],
                  complements=['F|y=0=0 contradicts irreducibility of quarticF','a0=0 or c=0 contradicts the nonzero restriction','alpha=beta=0 yields a polynomial in y alone and hence reducible projective quartic'],
                  external_dependency='Only the unique reduced ambient component consequence of legacy/v07/minimal_thickness.tex theorem minimal-full; nil exponent is not recomputed.',
                  scope='Whole nil211,translation0 irreducible-quartic support in MINIMAL_M; no other Jordan/translation types or global classification claimed.',
                  seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  reviewer_source_bytes=Path(__file__).read_text())
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key] for key in ('all_passed','seconds','maximum_RSS_KiB')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
