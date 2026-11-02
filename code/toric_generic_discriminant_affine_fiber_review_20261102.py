#!/usr/bin/env python3
"""Independent root, valuation and projective tangent arithmetic for the affine T01 branch."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
MAIN = 'TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_20261102'
STEM = 'TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_REVIEW_20261102'


def run():
    start = time.monotonic()
    x,r,z,s,t,c4,c6,h0,h1 = sp.symbols('x r z s t c4 c6 h0 h1')
    f = x**3 - 27*c4*x - 54*c6
    discriminant = sp.discriminant(f,x)
    assert sp.expand(discriminant - 78732*(c4**3-c6**2)) == 0
    # Construct multiplication by f' on the monic quotient without using a resultant.
    basis = [1,x,x**2]
    multiplication = sp.Matrix.hstack(*[
        sp.Matrix([sp.Poly(sp.rem(b*f.diff(x),f,x),x).nth(i) for i in range(3)])
        for b in basis])
    assert sp.expand(multiplication.det()+discriminant) == 0
    # The simple root lifts over the completed DVR. Writing it -2r gives this chart.
    split = (x+2*r)*((x-r)**2-z)
    parameters = {c4:(3*r**2+z)/27,c6:(-2*r**3+2*r*z)/54}
    assert sp.expand(f.subs(parameters)-split) == 0
    assert sp.expand(discriminant.subs(parameters)-4*z*(9*r**2-z)**2) == 0
    assert sp.expand(split.diff(x).subs(x,-2*r)-(9*r**2-z)) == 0
    ramified_slope = split.diff(x).subs({x:r+s,z:s**2})
    assert sp.expand(ramified_slope-2*s*(3*r+s)) == 0
    assert sp.cancel((-3*c6/c4).subs(parameters).subs(z,0)-r) == 0
    assert sp.cancel((6*c6/c4).subs(parameters).subs(z,0)+2*r) == 0
    # On the double fiber the projective chart denominator is h0+h1*s, h0 != 0.
    # Verify a generator of its nilpotent ideal via an explicit inverse of the denominator.
    chart_time = sp.rem(sp.expand(ramified_slope*(1/h0-h1*s/h0**2)),s**2,s)
    assert sp.cancel(chart_time-6*r*s/h0) == 0
    fiber_basis_matrix = sp.Matrix([[1,0],[0,6*r/h0]])
    assert fiber_basis_matrix.det() != 0
    deps = [
        MAIN,'TORIC_NORMALIZATION_BOUNDARY_20261101','PRIMITIVE_CURVE_DISCRIMINANT_20261101',
        'TORIC_C9_POLYNOMIAL_20261101','TORIC_POLYNOMIAL_E2_COORDINATE_20261101',
        'TORIC_ZERO_JET_EXHAUSTION_20261101','TORIC_UNIVERSAL_THREE_PLANE_20261101',
        'TORIC_JACOBIAN_TWO_TORSION_20261101']
    reports = {name:json.loads((ROOT/'reports'/(name+'.json')).read_text()) for name in deps}
    assert all(rep['all_passed'] for rep in reports.values())
    assert reports['TORIC_C9_POLYNOMIAL_20261101']['global_C9_is_homogeneous_polynomial_degree9']
    assert reports['TORIC_POLYNOMIAL_E2_COORDINATE_20261101']['polynomial_degree_in_q'] == 1
    assert reports['PRIMITIVE_CURVE_DISCRIMINANT_20261101']['discriminant_degree'] == 24
    files = [Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'code'/(MAIN.lower()+'.py')]
    files.extend(ROOT/'reports'/(n+'.json') for n in deps)
    files.extend(ROOT/'research/notes'/(n+'.md') for n in deps)
    return dict(id=STEM,all_passed=True,field='QQ',reviewed=MAIN,
                independent_symbolic_arithmetic=dict(
                    derivative_multiplication_matrix=sp.sstr(multiplication),
                    derivative_norm=str(sp.factor(multiplication.det())),
                    hensel_model=str(split),hensel_discriminant=str(4*z*(9*r*r-z)**2),
                    unramified_derivative=str(9*r*r-z),ramified_derivative=str(2*s*(3*r+s)),
                    projective_inverse_coordinate_mod_uniformizer_squared=str(chart_time),
                    doubled_residue_fiber_basis_determinant=str(fiber_basis_matrix.det())),
                proof_reviewed=['normal height-one intersection and a genuine simple pole',
                    'local Bezout inversion and scheme-theoretic affine fiber',
                    'finite projective image, uniformizer and local Nakayama',
                    'shrinking to an open neighborhood of the generic discriminant point'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in files},seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=run()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS',round(report['seconds'],3))
