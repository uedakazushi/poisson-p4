#!/usr/bin/env python3
"""Full binary orbit transversal with the two positive first-jet directions retained."""
import itertools
import json
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
from s61_model import S61, X
from s61_flag_transverse_obstruction_20261101 import specialize
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import sp_basis

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003
CUTS = [5,19,20,23,27]
A,B,D,E = sp.symbols('a b d e')
PIVOT_STEPS = [(71, 46), (73, 8), (75, 43), (70, 3), (0, 9), (3, 1), (9, 4), (62, 25), (65, 12), (69, 18), (72, 21), (74, 7), (87, 26), (93, 32), (123, 35), (103, 30), (81, 42), (83, 40), (1, 10), (2, 11), (6, 2), (8, 17), (29, 28), (49, 15), (64, 29), (90, 38), (92, 39), (99, 6), (80, 41), (102, 44), (82, 31), (4, 13), (5, 14), (7, 16)]


def load_report(stem):
    data = json.loads((ROOT/'reports'/f'{stem}.json').read_text())
    assert data['all_passed']
    return data


def original_transversal():
    model = S61()
    chart = load_report('S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102')
    old = load_report('S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101')['linear']
    fixed, columns = old['fixed'], old['cols']
    weights = chart['local_coordinate_weights']+[1,1]
    grades = chart['equation_weights']
    pi = coords()['Pi']
    y = sp.Matrix(X[1:])
    tensor = sp.zeros(4)
    tensor[1,2] = y[0]*y[3]
    tensor -= tensor.T
    orbit = []
    for g in sp_basis():
        derivative = -g*tensor-tensor*g.T+sum(((g*y)[i]*tensor.diff(y[i])
                                                for i in range(4)),sp.zeros(4))
        physical = sp.Matrix([sp.Poly(derivative[i,j],*y).coeff_monomial(y[h]*y[k])
                              for i,j in pairs for h,k in m2])
        orbit.append((pi*physical)[columns,:])
    motion = sp.Matrix.hstack(*orbit)
    assert list(motion.T.to_DM().rref()[1]) == CUTS
    selected = list(motion[CUTS,:].to_DM().rref()[1])
    determinant = motion[CUTS,selected].det()
    assert determinant != 0
    normal = sp.symbols('n0 n1')
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    beta = dict.fromkeys(one,sp.S.Zero)
    beta.update({3:A,13:B,35:D,48:E,11:normal[0],12:normal[1]})
    _, equations = specialize(model,beta,True)
    rows = []
    for grade,equation in zip(grades,equations):
        local = {}
        for monomial,coefficient in equation.items():
            # The chosen mixed coefficient equals one in this projective chart.
            indices = tuple(columns.index(i) for i in monomial if i != fixed)+(44,)*(2-len(monomial))
            if any(i in CUTS for i in indices):
                continue
            for (n0,n1),coefficient in sp.Poly(coefficient,*normal).terms():
                mon = tuple(sorted(indices+(45,)*n0+(46,)*n1))
                local[mon] = local.get(mon,0)+coefficient
        rows.append({mon:coefficient for mon,coefficient in local.items() if coefficient})
        assert all(sum(weights[j] for j in mon) == grade for mon in rows[-1])
    assert [weights[i] for i in CUTS].count(0) == 2
    assert all(w > 0 for i,w in enumerate(weights) if i not in CUTS)
    return rows,weights,grades,selected,determinant


def reduce_test_line(rows,weights,grades):
    """Keep b=6,d=12,e=16 and a indeterminate. Invert only constants on this line."""
    variables = [i for i in range(47) if i not in CUTS]
    free = [i for i in variables if i not in {j for _,j in PIVOT_STEPS}]
    assert sorted(weights[i] for i in free) == [1]*7+[2]
    field = sp.GF(PRIME)
    polynomials,*x = ring('a,'+','.join(f'x{i}' for i in range(8)),field)
    coordinates = [polynomials.zero]*47
    for j,v in enumerate(free):
        coordinates[v] = x[j+1]
    constants = {B:6,D:12,E:16}
    coefficient_cache = {}
    def convert(f):
        if f not in coefficient_cache:
            v = polynomials.zero
            for (power,),scalar in sp.Poly(f.subs(constants),A).terms():
                coefficient = int(scalar.p)*pow(int(scalar.q),-1,PRIME) % PRIME
                v += polynomials(coefficient)*x[0]**power
            coefficient_cache[f] = v
        return coefficient_cache[f]
    def evaluate(row):
        value = polynomials.zero
        for mon,c in row.items():
            term = convert(c)
            for i in mon:
                term *= coordinates[i]
            value += term
        return value

    determinants = []
    for grade in range(1,5):
        steps = [(i,j) for i,j in PIVOT_STEPS if weights[j] == grade]
        inds,cols = zip(*steps)
        block = sp.Matrix([[rows[i].get((j,),sp.S.Zero).subs(constants) for j in cols] for i in inds])
        determinant = sp.factor(block.det(method='domain-ge'))
        assert determinant.is_Rational and int(determinant.p) % PRIME and int(determinant.q) % PRIME
        determinants.append(str(determinant))
        inverse = block.to_DM().to_field().inv().to_Matrix()
        values = [evaluate(rows[i]) for i in inds]
        for j,index in enumerate(cols):
            coordinates[index] = -sum((convert(inverse[j,h])*value for h,value in enumerate(values)
                                       if value),polynomials.zero)
    residual = [evaluate(row) if grade <= 4 else polynomials.zero for row,grade in zip(rows,grades)]
    assert all(not residual[i] for i,_ in PIVOT_STEPS)
    return residual, [weights[i] for i in free], determinants


def monomials(weights, degree):
    def extend(i, remaining):
        if i == len(weights):
            if not remaining:
                yield ()
            return
        for k in range(remaining//weights[i]+1):
            for tail in extend(i+1,remaining-k*weights[i]):
                yield (k,)+tail
    return list(extend(0,degree))
