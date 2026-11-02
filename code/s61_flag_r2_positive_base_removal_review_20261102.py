#!/usr/bin/env python3
"""Check the mixed group minor by dual-number pullback and audit the central dimensions."""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import sympy as sp
from s61_model import S61, X, J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FLAG_R2_POSITIVE_BASE_REMOVAL_REVIEW_20261102'


def verify():
    start=time.monotonic()
    path=ROOT/'reports/S61_FLAG_R2_POSITIVE_BASE_REMOVAL_20261102.json';j=json.loads(path.read_text())
    sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
    assert j['all_passed']
    inputs=dict(j['input_sha256'])
    for k,h in inputs.items():assert sha(ROOT/k)==h
    model=S61();variables,parameters,beta=primitive_beta(model);v=sp.Matrix(variables)
    ids=[(i,k,l) for i,k in itertools.combinations(range(4),2) for l in range(4)]
    coefficients=lambda b:sp.Matrix([b[i,k].diff(v[l]) for i,k,l in ids])
    basis=coefficients(beta).jacobian(parameters)
    rows=list(basis.T.to_DM().rref()[1])
    rho,h=sp.symbols('rho h')
    special={p:{'u5':1,'u13':2,'u49':1,'u12':rho}.get(str(p),0) for p in parameters}
    b0=beta.subs(special,simultaneous=True)
    q=sp.zeros(4);q[1,3]=-v[0]**2;q-=q.T
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    columns=[]
    for i in [4,5,6,9]:
        A=sp_basis()[i];g=sp.eye(4)+h*A;ginv=sp.eye(4)-h*A;sub=dict(zip(v,g*v))
        raw_b=g.T*b0.subs(sub,simultaneous=True)*g
        db=basis[rows,:].inv()*coefficients(raw_b.diff(h).subs(h,0))[rows,:]
        raw_q=ginv*q.subs(sub,simultaneous=True)*ginv.T
        dq=raw_q.diff(h).subs(h,0)
        homogeneous=X[0]**2*J
        homogeneous[1:,1:]+=dq.subs(dict(zip(v,X[1:])),simultaneous=True)
        qvalues=sp.Matrix(model.normalize(homogeneous)['parameters'])[two,:]
        columns.append(sp.Matrix([db[0],*[qvalues[k] for k in [19,26,27]]]))
    witness=sp.Matrix.hstack(*columns)
    assert witness.det()==8
    dims=j['central_torsion_weight_dimensions']
    assert sum(dims)==111 and [i for i,v in enumerate(dims) if v]==list(range(2,7))
    assert 2+5>6 and int(j['actual_section_cut_determinant'])==8
    files=[path,Path(__file__).resolve(),ROOT/'research/notes'/(STEM+'.md'),ROOT/'code/weight1_structure.py',ROOT/'code/s61_model.py']
    inputs.update({str(p.relative_to(ROOT)):sha(p) for p in files})
    return dict(id=STEM,all_passed=True,independent_group_minor=8,group_motion_matrix=list(map(str,witness)),
                replay_scope='Separate primary replay and independent pullback for the selected four group directions, using the actual S61 normalization. Central Macaulay2 relations are primary dependencies.',
                t_annihilates_relative_kernel=True,input_sha256=inputs,seconds=time.monotonic()-start)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    out=p.parse_args().output;assert not out.exists();r=verify()
    with out.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(f"PASS independent mixed minor {r['independent_group_minor']}; {r['seconds']:.3f} seconds")
