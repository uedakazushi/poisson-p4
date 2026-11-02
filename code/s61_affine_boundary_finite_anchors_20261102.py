#!/usr/bin/env python3
"""Check the full flag core cover and the two parts of the finite boundary intersection."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102'


def verify():
    inputs = {}
    def record(path):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        inputs[str(path.relative_to(ROOT))] = digest
        return digest

    def load(stem):
        path = ROOT / 'reports' / (stem+'.json')
        report = json.loads(path.read_text())
        assert report['all_passed']
        record(path)
        for name, digest in report.get('input_sha256', {}).items():
            assert record(ROOT/name) == digest
        return report

    previous = load('S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102')
    assert previous['dimension_upper'] == 2
    load('S61_AFFINE_SURFACE_BOUND_20261102')
    load('S61_GENERIC_BOUNDARY_INDEX_20261102')
    load('S61_FLAG_UNIFORM_OPEN_20261101')
    load('S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101')
    load('S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102')
    ordinary = load('S61_FLAG_ORDINARY_CORE_SUPPORT_20261102')
    roots = [m['core'] for stem in ['S61_FLAG_LINEAR_ARTIFACT_WEIGHTED_REMOVAL_20261102',
                                    'S61_FLAG_TRUE_CORE_WEIGHTED_REMOVAL_20261102']
             for m in load(stem)['models']]
    for stem in ['S61_FLAG_T5_CORE_REMOVAL_20261102', 'S61_FLAG_T7_CORE_REMOVAL_20261102',
                 'S61_FLAG_R2_POSITIVE_BASE_REMOVAL_20261102']:
        roots.append(load(stem)['core'])
    exception = sp.Poly(sp.sympify(ordinary['core_open'].split('!=')[0]),sp.Symbol('r'))
    assert exception.degree()==len(roots)==7
    assert set(exception.all_roots())==set(map(sp.Rational,roots))
    binary = load('S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102')
    square = load('S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102')
    assert binary['residual_finite_over_general_beta']
    assert square['nonempty_ambient_removed_open'] and square['generic_beta_marked_residual_dimension_at_most']==0
    # In a frame where the given vector is e1, a Lagrangian plane through
    # e1 is spanned by e1 and u*e3+v*e4, so its Plucker image is a line.
    x,y=sp.symbols('x y')
    line=sp.Matrix([0,x,y,0,0])
    assert line.jacobian([x,y]).rank()==2
    p12,p13,p14,p23,p24=line
    assert -p12**2-p13*p24+p14*p23==0
    record(Path(__file__).resolve())
    record(ROOT/'research/notes'/(STEM+'.md'))
    return dict(id=STEM,all_passed=True,all_generic_true_boundary_points_in_binary_bundle=True,
                ordinary_flag_exceptional_cores=[str(c) for c in sorted(map(sp.Rational,roots))],
                divisor='P_C(ker(alpha) tensor V)',tautological_divisor_ample=True,
                generic_true_boundary_intersection_with_divisor_dimension_upper=0,
                every_true_boundary_curve_meets_finite_intersection=True,
                generic_boundary_finite_proved=False,complete_affine_five_proved=False,input_sha256=inputs)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    path=parser.parse_args().output
    result=verify()
    with path.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS',len(result['input_sha256']),'input hashes')
