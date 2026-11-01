#!/usr/bin/env python3
"""A smooth canonical curve at a binary infinity point of smaller Jacobian rank."""
import argparse
import hashlib
import json
import time
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1];m=S61()
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    data={4:1,5:-3,25:2,27:3,34:-2,37:-2,48:3,50:-3}
    beta={i:s.Integer(data.get(i+1,0)) for i in one}
    two,F=specialize(m,beta,True);fixed=two.index(29)
    q0=s.zeros(45,1);q0[fixed]=1
    leading=s.zeros(135,45);tcol=s.zeros(135,1)
    for r,f in enumerate(F):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q0[i]*q0[j]
                leading[r,i]+=c*q0[j];leading[r,j]+=c*q0[i]
            elif len(mon)==1:tcol[r,0]+=c*q0[mon[0]]
        assert value==0
    full=leading[:,[i for i in range(45) if i!=fixed]].row_join(tcol)
    ranks=[];minors=[];kernels=[]
    for A in (leading,full):
        _,cols=A.to_DM().rref();cols=list(cols)
        _,rows=A[:,cols].T.to_DM().rref();rows=list(rows)
        determinant=A.extract(rows,cols).det(method='domain-ge');assert determinant
        kernel=s.Matrix.hstack(*A.nullspace());assert A*kernel==s.zeros(A.rows,kernel.cols)
        ranks.append(len(cols));minors.append(dict(rows=rows,columns=cols,value=str(determinant)))
        kernels.append([[str(v) for v in row] for row in kernel.tolist()])
    assert ranks==[32,33]
    z,uu,form=primitive_beta(m);form=form.subs(dict(zip(uu,[beta[i] for i in one])))
    tensor=[form[i,j].diff(z[k]) for i,j in combinations(range(4),2) for k in range(4)]
    invariants=curve_invariants_from_beta(tensor)
    assert invariants[:3]==(22528,-2435968,3182443776)
    a,b,d=s.symbols('a b d');v1=s.Matrix([a,1,0,b]);v2=s.Matrix([-b,0,1,d])
    section=s.Matrix([(v1.T*form.subs(dict(zip(z,v1)))*v2)[0],(v1.T*form.subs(dict(zip(z,v2)))*v2)[0]])
    assert section.subs({a:0,b:0,d:0})==s.zeros(2,1)
    curve_jac=section.jacobian([a,b,d]).subs({a:0,b:0,d:0});assert curve_jac.rank()==2
    paths=[Path(__file__),root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',
           root/'code/primitive_jet_elliptic_curve_20261101.py',root/'code/elliptic_curve_invariants_20261101.py',
           root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_BINARY_RANK_DROP_SMOOTH_CURVE_20261101',all_passed=True,
        field='QQ',first_jet_nonzero_coordinates={f'u{k}':str(v) for k,v in data.items()},
        infinity_normal_form='Q=xw*d_y wedge d_z; weight2 coordinate u30=1',
        all135_leading_values_zero=True,leading_jacobian_rank=32,
        full_homogenized_jacobian_rank=33,full_projective_tangent_dimension=12,
        jacobian_minors=minors,jacobian_kernels=kernels,
        canonical_curve_derivative=[[str(v) for v in row] for row in curve_jac.tolist()],
        c4=str(invariants[0]),c6=str(invariants[1]),discriminant=str(invariants[2]),j=str(invariants[3]),
        conclusion='Smoothness of C_beta alone does not imply the generic binary Jacobian ranks34/35.',
        limitations=['No local-ring, t-saturation, or affine-branch claim at this exceptional point.',
                    'No claim that this rank-drop incidence dominates the first-jet space.'],
        proof_note='research/notes/S61_BINARY_RANK_DROP_SMOOTH_CURVE_20261101.md',
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','leading_jacobian_rank','full_homogenized_jacobian_rank','discriminant','seconds']}))


if __name__=='__main__':main()
