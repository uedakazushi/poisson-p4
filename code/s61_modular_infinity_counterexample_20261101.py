#!/usr/bin/env python3
"""A nonzero-modular point of135 unsaturated infinity with smooth C_beta.

Refutes a proposed modular-stabilizer identity on that leading locus. It says
nothing about saturation or whether this point is approached by affine fibers.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61,X,J as J5
from weight1_structure import COORDS,J,infinitesimal_action,minor_witness,sp_basis

ROOT=Path(__file__).resolve().parents[1]


def run():
    model=S61();x,y,z,w=X[1:];Q=s.zeros(4)
    entries={(0,1):-z*z,(0,2):0,(0,3):-(x*x+8*x*z+8*z*z)/2,
             (1,2):-z*z/2,(1,3):x*x+x*y+2*x*z+2*y*z+3*z*z,(2,3):z*(x+z)}
    for (i,j),value in entries.items():Q[i,j]=value;Q[j,i]=-value
    B=X[0]**2*J5
    B[1:,1:]+=Q
    normalization=model.normalize(B)
    beta_zero_verification=model.verify_point(normalization)
    weights=model.reduced['parameter_weights'];one=[i for i,v in enumerate(weights) if v==1]
    two=[i for i,v in enumerate(weights) if v==2]
    beta_values=list(map(s.Rational,[0,0,-24,784,32,24,18,0,1,'-311/3',2,'80/3',3,4,5,6]))
    point=list(normalization['parameters'])
    for i,value in zip(one,beta_values,strict=True):point[i]=value
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra=json.loads(extra_path.read_text())['universal_new_equation']
    top=[[term for term in f if sum(i in two for i in term[0])==2]
         for f in model.reduced['equations']+[extra]]
    residuals=[model.evaluate(f,point) for f in top]
    assert residuals==[0]*135
    variables=list(X[1:])
    modular=s.Matrix([sum(s.diff(Q[i,j],variables[j]) for j in range(4)) for i in range(4)])
    A=modular.jacobian(variables)
    assert A.rank()==2 and A*A==s.zeros(4) and A.T*J+J*A==s.zeros(4)
    zz,u,beta=primitive_beta(model)
    vector=s.Matrix([s.diff(beta[i,j],zz[k]) for (i,j),k in COORDS])
    primitive=vector.subs(dict(zip(u,beta_values,strict=True)))
    action=infinitesimal_action(primitive,A)
    assert action!=s.zeros(24,1)
    lie=s.Matrix.hstack(*[infinitesimal_action(primitive,aa) for aa in sp_basis()])
    lie_minor=minor_witness(lie);assert lie_minor['rank']==10
    c4,c6,delta,j=curve_invariants_from_beta(primitive)
    assert c4==47302703984349184 and c6==-6657339632224197100699648
    assert delta==35602891117000474561543155905702835990508339200
    assert (c4**3-c6**2)/1728==delta and delta!=0
    inputs=[Path(__file__),ROOT/'code/elliptic_curve_invariants_20261101.py',
            ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/s61_model.py',
            model.source/'symplectic_slice_reduced_equations.json',extra_path]
    return {'id':'S61_MODULAR_INFINITY_COUNTEREXAMPLE_20261101','all_passed':True,'field':'QQ',
        'Q':[[str(v) for v in row] for row in Q.tolist()],
        'source_normal_form':'old nil22 U1 representative','beta_zero_verification':beta_zero_verification,
        'weight_one_u_indices_zero_based':one,'weight_two_u_indices_zero_based':two,
        'beta_values':list(map(str,beta_values)),'infinity_S61_coordinates':list(map(str,point)),
        'all_135_leading_residuals':list(map(str,residuals)),
        'modular_vector_field':list(map(str,modular)),
        'modular_matrix':[[str(v) for v in row] for row in A.tolist()],
        'primitive_beta_coefficients':list(map(str,primitive)),
        'modular_action_on_beta':list(map(str,action)),
        'beta_sp4_lie_stabilizer_dimension':0,'sp4_orbit_minor':lie_minor,
        'c4':str(c4),'c6':str(c6),'Delta':str(delta),'j':str(j),
        'scope':'Counterexample on the UNSATURATED135-leading-equation locus, even Delta!=0. '
                'Does not claim an affine extension, a point of the saturated projective closure, '
                'or a component dominating W16.',
        'proof_note':'research/notes/S61_MODULAR_INFINITY_COUNTEREXAMPLE_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'leading_equations':135,'modular_rank':2,
                      'Delta_nonzero':True,'beta_lie_stabilizer_dimension':0,'output':str(args.output)}))


if __name__=='__main__':main()
