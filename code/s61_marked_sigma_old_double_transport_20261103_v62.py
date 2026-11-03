#!/usr/bin/env python3
"""Exact canonical pure-double-root frame at the old plane's discriminant point."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
from weight1_structure import COORDS,J

ROOT=Path(__file__).resolve().parents[1]


def verify(folder,output):
    start=time.monotonic();L=load_model(folder);K=s.QQ.algebraic_field(s.sqrt(-2),alias='h')
    hh=K.unit;h=K.to_sympy(hh);n=[K.zero]*13;n[9]=hh*hh;n[10]=hh;n[12]=coerce(K,s.Rational(2,3))*hh
    point=L.lift(hh,n,domain=K);assert point['all135_zero'] and not point['t']
    y,b,B=primitive_beta(L.model);y=s.Matrix(y)
    vector=lambda tensor:s.Matrix([tensor[i,j].diff(y[k]) for (i,j),k in COORDS])
    embedding=vector(B).jacobian(b);rows=list(embedding.T.to_DM().rref()[1]);inverse=embedding[rows,:].inv()
    C=s.Matrix([[-h,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,h]])
    D=s.diag(s.eye(2),s.Matrix([[1,-h/2],[2/(3*h),s.Rational(2,3)]]))
    F=(D*C).applyfunc(s.simplify);assert (F*J*F.T-J).applyfunc(s.simplify)==s.zeros(4)
    g=F.inv();old=B.subs(dict(zip(b,[K.to_sympy(value) for value in point['beta']])))
    moved=g.T*old.subs(dict(zip(y,g*y)),simultaneous=True)*g
    beta=[K.from_sympy(s.simplify(value)) for value in inverse*vector(moved)[rows,:]]
    physical=coords()['E']*s.Matrix([K.to_sympy(value) for value in point['q']]);Q=s.zeros(4)
    for k,(i,j) in enumerate(pairs):Q[i,j]=sum(physical[k*10+a]*y[z]*y[w] for a,(z,w) in enumerate(m2))
    Q-=Q.T
    canonical=(-2/h*F*Q.subs(dict(zip(y,F.inv()*y)),simultaneous=True)*F.T).applyfunc(s.simplify)
    expected=s.zeros(4);expected[1,2]=y[2]**2;expected[1,3]=-2*y[2]*y[3];expected-=expected.T
    assert canonical==expected
    coefficients=s.Matrix([s.Poly(canonical[i,j],*y).coeff_monomial(y[z]*y[w])
                           for i,j in pairs for z,w in m2])
    q=[K.from_sympy(value) for value in coords()['Pi']*coefficients]
    values=dict(zip(L.one,beta))|dict(zip(L.two,q))
    for form in L.original:
        value=K.zero
        for mon,numerator,denominator in form:
            if sum(index in L.two for index in mon)!=2:continue
            term=coerce(K,s.Rational(numerator,denominator))
            for index in mon:term*=values[index]
            value+=term
        assert not value
    def expr(value):return str(L.expression(K.from_sympy(value) if isinstance(value,s.Basic) else value,K,s.Symbol('h')))
    A,B0,U,V=[beta[L.one.index(i-1)] for i in [14,49,26,35]]
    assert A==K.one and B0==2*K.one and 2*A-B0==K.zero
    report=dict(id='S61_MARKED_SIGMA_OLD_DOUBLE_TRANSPORT_20261103_v62',all_passed=True,
                field='QQ(h),h^2=-2',source_plane_point={'u':'h','v':'2h/3'},
                new_coordinate_frame=[[expr(value) for value in row] for row in F.tolist()],
                projective_factor='-2/h',canonical_Hamiltonian='z^2*w',
                beta_nonzero={str(index+1):expr(value) for index,value in zip(L.one,beta) if value},
                q_nonzero={str(index):expr(value) for index,value in enumerate(q) if value},
                A=expr(A),B=expr(B0),U=expr(U),V=expr(V),
                old_double_root_table_factor='2A-B=0',all135_highest_equations_zero=True,
                actual_ambient_nilpotence_claimed=False,
                limitations=['This is the old table proper-beta-image stratum, not its generic free3 open.',
                             'No pointwise ambient nilpotence follows from the proper image bound.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
                                        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
                                        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']},
                seconds=time.monotonic()-start)
    with output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    r=verify(args.model_folder,args.output);print('PASS',r['id'],r['seconds'])
