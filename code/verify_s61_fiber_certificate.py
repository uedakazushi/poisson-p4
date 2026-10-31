#!/usr/bin/env python3
"""Independent centered-difference verification of saved full-fiber minors.

Weight-two variables occur with exponent at most two in each S61 equation.
Thus [f(p+e_j)-f(p-e_j)]/2 equals its partial derivative exactly. This verifier
does not call the sparse differentiation or family-normalization routines.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def centered_difference(equation, point, column):
    total = Fraction(0)
    for monomial, numerator, denominator in equation:
        if column not in monomial:
            continue  # The identical terms cancel in the two evaluations.
        assert monomial.count(column) <= 2
        plus = minus = Fraction(numerator, denominator)
        for index in monomial:
            plus *= point[index] + (1 if index == column else 0)
            minus *= point[index] - (1 if index == column else 0)
        total += (plus-minus)/2
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    source = ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    certificate = ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    data = json.loads(source.read_text())
    records = json.loads(certificate.read_text())
    for equation in data['equations']:
        for monomial, _, _ in equation:
            assert sum(data['parameter_weights'][index] == 2 for index in monomial) <= 2
    results = []
    for record in records['points']:
        point = [Fraction(x) for x in record['S61_coordinates_u1_through_u61']]
        witness = record['fiber_minor']
        rows, columns = witness['rows_zero_based'], witness['u_columns_zero_based']
        assert len(rows) == len(columns) == 45
        assert columns == [i for i,w in enumerate(data['parameter_weights']) if w == 2]
        matrix = [[centered_difference(data['equations'][row],point,column)
                   for column in columns] for row in rows]
        assert matrix == [[Fraction(x) for x in row] for row in witness['matrix']]
        determinant = s.Matrix(matrix).to_DM().det()
        assert Fraction(str(determinant)) == Fraction(witness['determinant']) != 0
        results.append({'component_id':record['component_id'],
                        'centered_difference_entries_matched':45*45,
                        'nonzero_determinant_matched':True})
    report = {'id':'S61_FIBER_INDEPENDENT_REVIEW_20261031','all_passed':True,
              'field':'QQ','method':'exact centered differences in each weight-two coordinate',
              'points':results,
              'limitations':['Verifies the reported minors from immutable equations and saved coordinates.',
                             'Family membership and component dimensions remain external inputs.',
                             'No other points or global exhaustion are certified.'],
              'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (Path(__file__),source,certificate)}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2)
        out.write('\n')
    print(json.dumps({'all_passed':True,'independently_matched_entries':6075,
                      'output':str(args.output)}))


if __name__ == '__main__':
    main()
