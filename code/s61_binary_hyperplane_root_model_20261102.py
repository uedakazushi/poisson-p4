"""Shared immutable coordinates and original-equation reconstruction for root certificates."""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import os

import sympy as sp
from s61_model import S61, X, J as SYMPLECTIC
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / 'code'
REPORTS = ROOT / 'reports'
P = 32003
CHART = REPORTS / 'S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json'
OLD = REPORTS / 'S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
G6 = REPORTS / 'S61_UNIVERSAL_WEIGHT6_20261031.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes(paths):
    return {str(Path(path).resolve().relative_to(ROOT)): digest(path) for path in paths}


def source_paths():
    legacy = ROOT / 'legacy/v02/data'
    return [
        Path(__file__), CODE/'s61_model.py', CODE/'s61_flag_transverse_obstruction_20261101.py',
        CODE/'primitive_jet_elliptic_curve_20261101.py', CHART, OLD, G6,
        legacy/'symplectic_slice_certificate.json', legacy/'symplectic_slice_equations.json',
        legacy/'symplectic_slice_reduced_equations.json',
    ]


def output_argument(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--output', required=True, type=Path)
    path = parser.parse_args().output
    if path.exists():
        parser.error(f'output already exists: {path}')
    return path


def write_new(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def runtime_environment():
    return os.environ | {'DYLD_LIBRARY_PATH': str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}


def residue(value):
    q = Fraction(value)
    assert q.denominator % P
    return q.numerator * pow(q.denominator, -1, P) % P


def original_point(model, a, e=2, d=3):
    """Return all135 QQ local polynomials before the three cuts."""
    chart = json.loads(CHART.read_text())
    old = json.loads(OLD.read_text())['linear']
    substitution = dict(zip(sp.symbols('a e d'), map(sp.Rational, [a,e,d])))
    beta = {int(k[1:])-1: sp.sympify(v).subs(substitution) for k,v in chart['beta'].items()}
    two, source = specialize(model, beta, True)
    assert len(source) == 135
    fixed, cols = old['fixed'], old['cols']
    weights, row_weights = chart['local_coordinate_weights'], chart['equation_weights']
    equations = []
    for row, encoded in enumerate(source):
        q = {}
        for monomial, value in encoded.items():
            value = Fraction(str(value))
            if len(monomial) == 2 and monomial.count(fixed) == 1:
                local = (cols.index(next(j for j in monomial if j != fixed)),)
            elif monomial == (fixed,):
                local = (44,)
            elif fixed not in monomial:
                local = tuple(sorted([cols.index(j) for j in monomial] + [44]*(2-len(monomial))))
            else:
                assert not value
                continue
            q[local] = q.get(local, Fraction(0)) + value
        q = {m:c for m,c in q.items() if c}
        assert all(sum(weights[j] for j in m)==row_weights[row] for m in q)
        equations.append(q)
    return beta, two, old, weights, row_weights, equations


def old_cut_rows(old):
    inverse = sp.Matrix(old['chart_basis']).inv()
    raw = []
    for i in old['parameters']:
        row = inverse.row(old['free'][i])
        nonzero = [(j,v) for j,v in enumerate(row) if v]
        assert len(nonzero) == 1 and nonzero[0][1] == 1
        raw.append(nonzero[0][0])
    assert raw == [5,19,23]
    return raw


def actual_binary_tangent(model, beta, two, old):
    z, _, form = primitive_beta(model)
    form = form.subs({sp.Symbol(f'u{i+1}'):v for i,v in beta.items()})
    A,C,D,h,k = sp.symbols('planeA planeC planeD conicH conicK')
    origin = dict.fromkeys((A,C,D,h,k), 0)
    v1, v2 = sp.Matrix([A,1,0,C]), sp.Matrix([-C,0,1,D])
    section = sp.Matrix([(v1.T*form.subs(dict(zip(z,v)))*v2)[0] for v in (v1,v2)])
    assert section.subs(origin) == sp.zeros(2,1)
    jacobian = section.jacobian([A,C,D]).subs(origin)
    assert jacobian.rank() == 2
    direction = jacobian.row(0).cross(jacobian.row(1)).reshape(3,1)
    x,y,z,w = X[1:]
    l1,l2 = x-A*y+C*z, w-C*y-D*z
    family = (l1*l2+h*l1**2+k*l2**2)*(v1*v2.T-v2*v1.T)
    fixed, cols = old['fixed'], old['cols']
    center = sp.zeros(45,1)
    center[fixed] = 1
    derivatives = []
    for variable in (A,C,D,h,k):
        bivector = X[0]**2*SYMPLECTIC
        bivector[1:,1:] += family.diff(variable).subs(origin)
        parameters = model.normalize(bivector)['parameters']
        point = sp.Matrix([parameters[i] for i in two])
        point -= center*point[fixed]
        derivatives.append(point[cols,:].col_join(sp.zeros(1,1)))
    tangent = sp.Matrix.hstack(sp.Matrix.hstack(*derivatives[:3])*direction, *derivatives[3:])
    return jacobian, direction, tangent
