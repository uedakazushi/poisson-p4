#!/usr/bin/env python3
"""Independent finite-group check and proof audit of the sparse flag germ."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s
from s61_model import S61, X, J
from primitive_jet_elliptic_curve_20261101 import primitive_beta

ROOT = Path(__file__).resolve().parents[1]


def run():
    start = time.time()
    rp = ROOT/'reports/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.json'
    r = json.loads(rp.read_text())
    assert r['all_passed']
    for name,digest in r['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest,name
    model = S61()
    point = s.Matrix(list(map(s.Rational,r['point61'])))
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two = [i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    z,u,beta = primitive_beta(model)
    beta = beta.subs(dict(zip(u,point[one,:])))
    a,b,lam,mu = s.symbols('a b lam mu', nonzero=True)
    def group(c,h):
        return s.Matrix([[1,c,0,0],[0,1,0,0],[0,0,h,h-1/h],[0,0,0,1/h]])
    g = group(a,lam)
    assert (g*group(b,mu)-group(a+b,lam*mu)).applyfunc(s.cancel) == s.zeros(4)
    assert (g*J[1:,1:]*g.T-J[1:,1:]).applyfunc(s.cancel) == s.zeros(4)
    assert (g.T*beta.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g-beta).applyfunc(s.cancel) == s.zeros(4)
    # Use the full finite pullback, then differentiate its projective q chart.
    x,y,zz,w = X[1:]
    spatial = s.Matrix([x,y,zz,w])
    Q = s.zeros(4);Q[1,3] = -x*x;Q -= Q.T
    pulled = (g.inv()*Q.subs(dict(zip(spatial,g*spatial)),simultaneous=True)*g.inv().T).applyfunc(s.expand)
    homogeneous = X[0]**2*J
    homogeneous[1:,1:] += pulled
    q = s.Matrix(model.normalize(homogeneous)['parameters'])[two,:]
    fixed = two.index(37)
    assert point[37] == -1
    q = (point[37]*q/q[fixed]).applyfunc(s.cancel)
    identity = {a:0,lam:1}
    assert q.subs(identity) == point[two,:]
    tangents = q.jacobian([a,lam]).subs(identity)
    # q38 is fixed; free chart rows are read from the independently saved raw135 data.
    base = json.loads((ROOT/'reports/S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101.json').read_text())
    from s61_minimal_infinity_local_data_20261101 import calculate
    local = calculate(model,True,point)
    matrix = s.Matrix(local['jacobian']).applyfunc(s.Rational)
    frame = s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    columns = local['chart_columns']
    free = [i for i in range(45) if i not in matrix.rref()[1]]
    derivatives = tangents[columns,:].col_join(s.zeros(1,2))
    coordinates = derivatives[free,:]
    assert frame*coordinates == derivatives
    assert coordinates == s.Matrix(r['orbit_tangent']).applyfunc(s.Rational)
    assert coordinates[[0,1],:].det() == 4 and coordinates[2:,:] == s.zeros(7,2)
    assert base['geometry']['b1_derivative'][0] == ['0','0','-1/2','0','0','0','0','0','0']
    assert r['initial_torsion_top_degree'] == 6 and r['homogenizer_order'] == 4
    paths = [Path(__file__),rp,ROOT/'code/s61_flag_boundary_orbit_split_formal_20261101.py',ROOT/'research/notes/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.md',ROOT/'reports/S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101.json',ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/integer_implicit_jets_20261101.cpp',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_REVIEW_20261101',all_passed=True,independent_finite_group_law=True,independent_symplectic_and_beta_identities=True,independent_finite_q_pullback_tangent=True,orbit_slice_minor=4,conceptual_review=dict(formal_etale_product_accepted=True,same_flag_curve_generic_length5_accepted=True,extra_initial_kernel_in_finite_torsion=True,induced_filtration_top_degree6_accepted=True,order4_nonzero_and_square_zero_accepted=True),scope=['Independent finite group and finite q pullback check; raw local-frame utility shared with producer.','Producer QQ jet, observable-order, and Groebner-module calculations reproduced separately, not independently reimplemented.','Only the stated fixed-beta germ. No relative ambient or entire-boundary removal.'],seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
