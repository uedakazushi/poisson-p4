#!/usr/bin/env python3
"""Verify the actual fixed-beta P1 and both d=0 pencil charts over Z_(32003)."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_flag_transverse_obstruction_20261101 import specialize
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
from weight1_structure import COORDS,J

ROOT=Path(__file__).resolve().parents[1]
P=32003


def verify():
    start=time.monotonic()
    model=S61(); y,u,tensor=primitive_beta(model); y=s.Matrix(y)
    source=ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    flag=json.loads(source.read_text())['flag']
    H=s.Matrix(flag['beta_kernel']).applyfunc(s.Rational)
    core=s.zeros(15,1)
    core[3],core[12],core[5],core[6],core[14]=s.Rational(4,5),1,1,3,2
    beta0=H*core
    beta_tensor=tensor.subs(dict(zip(u,beta0)))
    vector=lambda f:s.Matrix([f[i,j].diff(y[k]) for (i,j),k in COORDS])
    embedding=vector(tensor).jacobian(u)
    pivot=list(embedding.T.to_DM().rref()[1]);inverse=embedding[pivot,:].inv()
    a=s.Symbol('m')
    exchange=s.eye(4);exchange[2,2]=exchange[3,3]=0
    exchange[2,3]=-1;exchange[3,2]=1
    U=s.Matrix([[-s.Rational(5,18),s.Rational(10,3)],[s.Rational(1,3),0]])*s.Matrix([[1,0],[-a/12,1]])
    sign=s.diag(-1,1);V=sign*U.inv().T*sign
    group=s.zeros(4)
    for indices,block in [([0,3],V),([1,2],U)]:
        for i in range(2):
            for j in range(2):
                group[indices[i],indices[j]]=block[i,j]
    group=exchange*group
    assert (group.T*J*group-J).applyfunc(s.expand)==s.zeros(4)
    for value in list(group)+list(group.inv()):
        for coefficient in s.Poly(value,a,domain=s.QQ).coeffs():
            assert int(coefficient.q)%P
    pulled=group.T*beta_tensor.subs(dict(zip(y,group*y)),simultaneous=True)*group
    beta=(inverse*vector(pulled)[pivot,:]).applyfunc(s.factor)
    expected=s.zeros(16,1)
    expected[0]=-(5*a*a+7)/120;expected[4]=(a+2)/9
    expected[9]=1;expected[13]=(5*a+1)/9;expected[14]=-5
    assert (beta-expected).applyfunc(s.cancel)==s.zeros(16,1)
    q0=s.zeros(4);q0[1,2]=y[0]*y[3];q0-=q0.T
    moving=(group*q0.subs(dict(zip(y,group.inv()*y)),simultaneous=True)*group.T).applyfunc(s.factor)
    physical=s.Matrix([s.Poly(moving[i,j],*y).coeff_monomial(y[k]*y[l])
                       for i,j in pairs for k,l in m2])
    q=(coords()['Pi']*physical).applyfunc(s.expand)
    assert all(s.Poly(value,a).degree()<=1 for value in q)
    line=s.Matrix([[value.coeff(a,0),value.coeff(a,1)] for value in q])
    rows=list(line.T.to_DM().rref()[1]);assert len(rows)==2
    minor=s.Rational(line[rows,:].det());assert int(minor.p)%P and int(minor.q)%P
    top=line[:,1]
    assert [i for i,value in enumerate(top) if value]==[25]
    top_value=s.Rational(top[25]);assert int(top_value.p)%P and int(top_value.q)%P
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    _,forms=specialize(model,dict(zip(one,beta0)),True)
    for form in forms:
        value=sum((coefficient*q[mon[0]]*q[mon[1]] for mon,coefficient in form.items()
                   if len(mon)==2),s.S.Zero)
        assert s.expand(value)==0
    paths=[Path(__file__),source,ROOT/'code/s61_model.py',
           ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
           ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
           ROOT/'code/toric_closed_jet_coordinates_20261102.py',
           ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_BINARY_D_ZERO_GEOMETRY_20261103',all_passed=True,prime=P,
                finite_beta=list(map(str,beta)),fixed_beta=list(map(str,beta0)),
                finite_symplectic_matrix=[[str(value) for value in row] for row in group.tolist()],
                moving_tensor=[str(moving[i,j]) for i,j in pairs],
                moving_q=list(map(str,q)),projective_line_minor=str(minor),
                projective_line_minor_rows=rows,infinity_flag_q_index=25,
                infinity_leading_coefficient=str(top_value),
                fixed_beta_highest_135_identities=True,
                flat_projective_line_over_V=True,
                limitations=['Checks this fixed-beta pencil and its true flag endpoint.',
                             'Whole support-plane factor locus is a separate geometric dependency.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'],report['seconds'])
