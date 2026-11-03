#!/usr/bin/env python3
"""Rebuild old-curve augmented matrices and independently reduce row modules.

The Python reduction chooses the first leading position, independently of the
native producer's last-leading-position convention.  Its determinant support
is a finite-field all-parameter certificate, not an all-complex lift assertion.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM

from s61_marked_sigma_point_lift_20261103_v62 import load_model, coerce

ROOT = Path(__file__).resolve().parents[1]
P = 32003


def dm(rows, domain):
    return DomainMatrix.from_rep(DDM(rows, (len(rows), len(rows[0])), domain))


def row_determinant(rows, R):
    n = len(rows[0])
    basis = [None]*n
    def leading(row):
        choices = [(value.degree(), -j, j) for j, value in enumerate(row) if value]
        if not choices:
            return None
        degree, _, position = max(choices)
        return degree, position
    for raw in rows:
        row = list(raw)
        while (lead := leading(row)) is not None:
            degree, position = lead
            if basis[position] is None:
                basis[position] = row
                break
            old_degree = leading(basis[position])[0]
            if degree < old_degree:
                row, basis[position] = basis[position], row
                degree, old_degree = old_degree, degree
            factor = (row[position].LC/basis[position][position].LC)*R.gens[0]**(degree-old_degree)
            row = [a-factor*b for a, b in zip(row, basis[position])]
    assert all(row is not None for row in basis)
    positive = [i for i, row in enumerate(basis) if leading(row)[0] > 0]
    constant = [i for i in range(n) if i not in positive]
    residual = [list(basis[i]) for i in positive]
    for c in constant:  # first-leading convention: eliminate left to right
        pivot = basis[c][c]
        assert pivot.degree() == 0
        for row in residual:
            factor = row[c]/pivot.LC
            for j in range(n):
                row[j] -= factor*basis[c][j]
            assert not row[c]
    determinant = dm([[row[j] for j in positive] for row in residual], R.to_domain()).det()
    assert determinant.degree() == sum(leading(row)[0] for row in basis)
    return determinant.monic(), len(positive)


def verify(model, output):
    folder = ROOT/'reports/S61_MARKED_SIGMA_OLD_CURVE_COVER_20261103_v62'
    report_path = folder/'report.json'
    report = json.loads(report_path.read_text())
    ambient_path = ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_AMBIENT_20261103_v62/report.json'
    ambient = next(c for c in json.loads(ambient_path.read_text())['cases']
                   if c['case'] == 'distinct-u' and c['beta'] == 'actual')
    L = load_model(model)
    k = s.Symbol('k')
    ground = s.GF(P)
    F = ground.frac_field(k)
    R = F.field.ring
    results = []
    inputs = [Path(__file__), report_path, ambient_path]
    for case in report['results']:
        root = case['h_mod_P']
        normals = [F.zero]*13
        normals[9], normals[10], normals[12] = coerce(F, root), F.one, F.gens[0]
        point = L.lift(root, normals, domain=F)
        assert point['all135_zero'] and not point['t']
        frozen = {L.two[25], *[L.two[i] for i in ambient['q_cut_indices']],
                  *[L.one[i] for i in ambient['beta_cut_indices']]}
        variables = [i for i in range(62) if i not in frozen]
        assert len(variables) == 44
        direction = [F.from_sympy(s.sympify(v.replace(' mod 32003', ''))) for v in case['normal_direction_polynomials']]
        assert direction[10] and not direction[10].numer.degree()
        values = dict(zip(L.one, point['beta'])) | dict(zip(L.two, point['q'])) | {61:F.zero}
        delta = dict(zip(variables, direction))
        linear = [[F.zero]*44 for _ in range(135)]
        quadratic = []
        for row, equation in enumerate(L.original):
            qtotal = F.zero
            for mon, n, d in equation:
                ids = list(mon)+[61]*(2-sum(i in L.two for i in mon))
                coefficient = coerce(F, s.Rational(n, d))
                for variable, power in Counter(ids).items():
                    if variable not in delta:
                        continue
                    rest = list(ids);rest.remove(variable)
                    term = coefficient*power
                    for i in rest:
                        term *= values[i]
                    linear[row][variables.index(variable)] += term
                series = [coefficient, F.zero, F.zero]
                for i in ids:
                    v, e = values[i], delta.get(i, F.zero)
                    series = [series[0]*v, series[1]*v+series[0]*e, series[2]*v+series[1]*e]
                qtotal += series[2]
            quadratic.append(qtotal)
        assert all(not sum((a*b for a, b in zip(row, direction)), F.zero) for row in linear)
        assert len(case['records']) == 1 and case['records'][0]['omitted_column'] == 10
        path = folder/f'h{root}-omit10.in'
        tokens = iter(path.read_text().split())
        number = lambda: int(next(tokens))
        assert [number() for _ in range(6)] == [P,135,44,0,0,1]
        saved = []
        for _ in range(135):
            row = [R.zero]*44
            for _ in range(number()):
                column, count = number(), number()
                row[column] = R.from_dict({(i,):ground.convert(number()) for i in range(count)})
            saved.append(row)
        assert next(tokens, None) is None
        for index, row in enumerate(linear):
            expected = row[:10]+row[11:]+[quadratic[index]]
            for column, (old, new) in enumerate(zip(saved[index], expected)):
                # SymPy leaves nonmonic constant denominators in GF(p)(k),
                # so structural FracElement equality can be false for equal values.
                assert not (F.new(old)-new), (root, index, column, str(old), str(new))
        determinant, size = row_determinant(saved, R)
        discriminant = R.from_expr(9*k*k+2*root).monic()
        assert determinant == discriminant**2
        results.append(dict(h_mod_P=root, original135_augmented_matrix_matches=True,
                            omitted_direction_component_is_constant_unit=True,
                            independent_rowmodule_determinant=str(determinant),
                            independent_residual_size=size,
                            determinant_equals_discriminant_square=True,
                            all_finite_geometric_special_points_off_discriminant_upper2=True))
        inputs.append(path)
        print('PASS', root, 'DETERMINANT', determinant, flush=True)
    for path, digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    result = dict(id='S61_MARKED_SIGMA_OLD_CURVE_COVER_20261103_v62_REVIEW', all_passed=True,
                  results=results, all_complex_parameter_scope_accepted=False,
                  remaining_lift_requirements=['Special-fiber discriminant-zero ambient removal.',
                                               'Infinity ambient removal.',
                                               'Proper integral parameter family covering all weighted projective ratios, or a separate exact characteristic-zero all-parameter certificate.'],
                  input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs})
    with output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-folder', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.model_folder, args.output)
