#!/usr/bin/env python3
"""Affine-linear change of the normalized primitive jet under basepoint motion.

This constructs transport, not new Poisson-ideal constraints. In particular,
constancy of an elliptic-curve invariant under this transport is not assumed.
"""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X
from elliptic_weight1_dominance import projective_affine
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J,COORDS,constraints
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'


def primitive_projection(linear,z):
    tau=linear[0,1]+linear[2,3];alpha=s.Matrix([s.diff(tau,zi)/5 for zi in z]);v=J*alpha;zz=s.Matrix(z)
    return (-J*(linear-2*alpha.dot(zz)*J+zz*v.T-v*zz.T)*J).applyfunc(s.expand)


def transport(model,point):
    weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    q=s.symbols('q0:45');values=list(point)
    for i,v in zip(two,q,strict=True):values[i]=v
    z=s.symbols('z1:5');zero=dict.fromkeys(z,0)
    polynomial=projective_affine(model.reconstruct(values),X,z)
    assert polynomial.subs(zero)==J
    linear=[polynomial.diff(zi).subs(zero) for zi in z]
    z0,u,beta=primitive_beta(model)
    beta=beta.subs(dict(zip(z0,z)),simultaneous=True)
    basis=s.Matrix([s.diff(beta[i,j],z[k]) for (i,j),k in COORDS]).jacobian(u)
    pivots=list(basis.T.rref()[1]);left_inverse=basis[pivots,:].inv()
    constraint=s.Matrix.vstack(*constraints());results=[];constraints_checked=True
    for h in range(4):
        frame=s.Rational(1,2)*linear[h]*J
        assert frame*J+J*frame.T+linear[h]==s.zeros(4)
        variation=[]
        for i in range(4):
            value=polynomial.diff(z[h],z[i]).subs(zero)+frame*linear[i]+linear[i]*frame.T
            value-=sum((frame[k,i]*linear[k] for k in range(4)),s.zeros(4))
            variation.append(value)
        variation_form=primitive_projection(sum((z[i]*variation[i] for i in range(4)),s.zeros(4)),z)
        flat=s.Matrix([s.diff(variation_form[i,j],z[k]) for (i,j),k in COORDS])
        assert (constraint*flat).applyfunc(s.expand)==s.zeros(8,1)
        coordinate=(left_inverse*flat[pivots,0]).applyfunc(s.expand)
        assert (basis*coordinate-flat).applyfunc(s.expand)==s.zeros(24,1)
        assert all(s.Poly(v,*q).total_degree()<=1 for v in coordinate)
        coefficients=coordinate.jacobian(q);constant=coordinate.subs(dict.fromkeys(q,0))
        results.append({'direction':h,'coefficients':coefficients,'constant':constant})
    return one,two,q,basis,results


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    record=next(p for p in json.loads(SOURCE.read_text())['points'] if p['component_id']=='T01')
    point=list(map(s.Rational,record['S61_coordinates_u1_through_u61']));model=S61()
    one,two,q,basis,result=transport(model,point)
    matrix=s.Matrix.vstack(*(entry['coefficients'] for entry in result));rank=matrix.rank()
    report={'id':'S61_JET_BASEPOINT_TRANSPORT_20261101','field':'QQ','all_passed':True,
        'fixed_primitive_u_order':[i+1 for i in one],'fixed_u_values':[str(point[i]) for i in one],
        'q_to_u_indices_zero_based':two,'q_coordinates':'uncentered original weight2 coordinates',
        'primitive_coefficient_basis':[[str(v) for v in row] for row in basis.tolist()],
        'derivative_by_basepoint_direction':[{'direction':e['direction'],'q_coefficient_matrix':[[str(v) for v in row] for row in e['coefficients'].tolist()],'constant':list(map(str,e['constant']))} for e in result],
        'combined_q_rank':rank,
        'formula':'A_h=(1/2)L_h J; deltaL_i=P_hi(0)+A_h L_i+L_i A_h^T-sum_k (A_h)_ki L_k; project to primitive beta.',
        'scope':'Only exact basepoint transport; no invariant constancy or new ideal membership is claimed.',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),SOURCE,ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'combined_q_rank':rank,'output':str(args.output)}))

if __name__=='__main__':main()
