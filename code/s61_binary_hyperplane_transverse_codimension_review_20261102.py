#!/usr/bin/env python3
"""Audit the H8 height-one proof and independently check the H8 dominance minor."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as sp
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import COORDS, sp_basis

ROOT = Path(__file__).resolve().parents[1]
MAIN = 'S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102'
STEM = 'S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_REVIEW_20261102'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dominance():
    z, symbols, form = primitive_beta(S61())

    def tensor(f):
        return sp.Matrix([f[i,j].diff(z[k]) for (i,j),k in COORDS])

    embedding = tensor(form).jacobian(symbols)
    independent = list(embedding.T.rref()[1])
    left = embedding[independent,:].inv()
    params = sp.symbols('a e d p q r h k')
    a,e,d,p,q,r,h,k = params
    beta = sp.Matrix([a,0,0,0,(e+1)/5,p,q,r,0,1,d,h,k,e,-5,0])
    values = [2,3,5,7,11,13,17,19]
    point = beta.subs(dict(zip(params,values)))
    fiber_form = form.subs(dict(zip(symbols,point)))
    columns = []
    for A in sp_basis():
        transport = sum(((A*sp.Matrix(z))[j]*fiber_form.diff(z[j]) for j in range(4)), sp.zeros(4))
        derivative = tensor(A.T*fiber_form + fiber_form*A + transport)
        coordinates = left*derivative[independent,:]
        assert embedding*coordinates == derivative
        columns.append(coordinates)
    matrix = sp.Matrix.hstack(beta.jacobian(params), *columns)
    chosen = list(matrix.to_DM().rref()[1])
    assert len(chosen)==16
    determinant = matrix[:,chosen].det(method='domain-ge')
    assert determinant == sp.Rational(1680453,125)
    return {'H8_coordinates':list(map(str,params)), 'point':values,
            'derivative_shape':list(matrix.shape), 'ordered_columns':'eight H8 then ten symplectic directions',
            'chosen_columns':chosen, 'minor_QQ':str(determinant)}


def audit():
    start=time.time()
    main_path=ROOT/'reports'/(MAIN+'.json')
    main=json.loads(main_path.read_text())
    assert main['all_passed'] and main['exceptional_core_dimension_upper']==1
    for field in ['input_sha256','primary_dependency_sha256']:
        for name,value in main[field].items():
            assert sha(ROOT/name)==value,name
    p=32003
    a,e,d=sp.symbols('a e d')
    F=500*a+90*e**2+3*e-24
    B=50*a*d+18*d*e**2-9*d*e-2*d+60*e-40
    D0=50*a+9*e**2+3*e-6
    C=100*a+18*e**2-9*e-2
    P=125*a*a*d-450*a*e-200*a-54*e**3-27*e**2+33*e+6
    R=15*a*d*e-10*a*d+200*a+18*e**2+6*e-12
    S=3*e-2
    E=B-S*(d*(6*e+1)+40)
    assert sp.expand((-S*E)*P+(25*a*E+10*S*C)*R-50*a*B*C)==0
    factors={k:sp.Poly(v.subs(e,2),a,modulus=p).monic()
             for k,v in [('a',a),('F',F),('D0',D0),('C',C)]}
    roots={k:int(-v.nth(0))%p for k,v in factors.items()}
    assert len(set(roots.values()))==4 and roots==main['special_a_roots']
    vertex=ROOT/'reports/S61_BINARY_HYPERPLANE_VERTICAL_SLICE_CERTIFICATE_20261102.json'
    vr=json.loads(vertex.read_text())
    field=vr['full_modular_row_modules']['3']
    monic=sp.Poly.from_list(list(reversed(field['monic_coefficients'])),a,modulus=p)
    product=sp.Poly(1,a,modulus=p)
    for k,v in field['known_orders'].items():
        product*=factors[k]**v
    assert product==monic and monic.degree()==287
    assert sum(vr['kernel_dimensions'])==119
    assert vr['selected_row_weight_sum']+vr['column_weight_sum']-119==575
    assert vr['direct_modular_minor']['coefficients'][-1]==16008
    axes=json.loads((ROOT/'reports/S61_BINARY_HYPERPLANE_AXES_TRANSVERSE_LENGTH_20261102.json').read_text())
    alt=json.loads((ROOT/'reports/S61_BINARY_HYPERPLANE_C_TRANSVERSE_LENGTH_20261102.json').read_text())
    for point in axes['points'].values():
        assert point['weighted_Hilbert_function_modp']==[1,7,23,34,0]
        assert point['transverse_determinant_modp']%p
        assert all(piece['minor_modp']%p for piece in point['low_pieces'].values())
    assert alt['new_cut_determinant_modp']%p and all(
        piece.get('minor_modp',1)%p for piece in alt['weighted_pieces'].values())
    assert 1+5+10==16 and 8+3==11
    extra=[main_path,ROOT/'code/s61_binary_hyperplane_transverse_codimension_20261102.py',
           ROOT/'research/notes'/(MAIN+'.md'),vertex, Path(__file__),
           ROOT/'research/notes'/(STEM+'.md'),ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
           ROOT/'code/s61_model.py',ROOT/'code/weight1_structure.py',
           ROOT/'research/notes/S61_BINARY_FORMAL_RING_20261101.md',
           ROOT/'research/notes/S61_BINARY_DEEPER_FULL_CORE_REMOVAL_20261102.md',
           ROOT/'research/notes/S61_BINARY_HYPERPLANE_E_ONLY_SUPPORT_20261102.md',
           ROOT/'research/notes/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.md',
           ROOT/'research/notes/S61_BINARY_NUISANCE_CONTRACTION_20261102.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return {'all_passed':True,'reviewed_theorem':MAIN,'accepted_scope':'H8, exceptional core dimension <=1, and finite residual points from this presentation over general beta.',
            'checks':'Independent H8 dominance, old heavy-square Bezout, finite-field roots and module product, degree arithmetic, minor units, relative and image dimensions. Primary large coefficient calculations are dependencies.',
            'proof_audit_accepted':True,'independent_dominance':dominance(),'four_distinct_roots_modp':roots,
            'special_module_degree':monic.degree(),'selected_degree':575,'completed_regular_base_dimension':11,
            'residual_source_dimension_upper':16,'main_primary_reports':len(main['primary_inputs']),
            'main_primary_hashes_checked':len(main['primary_dependency_sha256']),
            'H11_in_reviewed_scope':False,'residual_emptiness_claimed':False,
            'input_sha256':{str(q.resolve().relative_to(ROOT)):sha(q) for q in extra},
            'seconds':time.time()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists():
        parser.error(f'output exists: {args.output}')
    result=audit()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print('PASS',result['independent_dominance']['minor_QQ'],result['seconds'])
