#!/usr/bin/env python3
"""Pointwise conditions for the rank65 finite free binary completed local ring."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1]
    paths={name:root/f'reports/{stem}_20261101.json' for name,stem in [
        ('quadratic','S61_BINARY_QUADRATIC_OBSTRUCTION'),('local','S61_BINARY_LOCAL_DIMENSION'),
        ('transverse','S61_BINARY_TRANSVERSE_ALGEBRA'),('uniform_rank','S61_BINARY_UNIFORM_RANK'),
        ('uniform_quadratic','S61_BINARY_UNIFORM_QUADRATIC'),('uniform_cubic','S61_BINARY_UNIFORM_CUBIC')]}
    reports={name:json.loads(path.read_text()) for name,path in paths.items()}
    assert all(r['all_passed'] for r in reports.values())
    source=reports['quadratic'];local=reports['local'];transverse=reports['transverse']
    assert source['point']==transverse['point']
    assert transverse['transverse_algebra_length']==65 and transverse['transverse_infinity_quotient_length']==64
    assert transverse['homogenizing_class_nonzero'] and transverse['homogenizing_class_in_cubic_socle']
    assert reports['uniform_cubic']['normal_length_lower_bound']==65
    assert reports['uniform_cubic']['cubic_relation_rank_upper_bound']==50
    assert 44 not in reports['uniform_cubic']['lift_support_columns']
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    point=list(map(s.Rational,source['point']));beta={i:point[i] for i in one}
    assert [beta[i] for i in one]==[s.S.Zero if i in (11,12) else s.Integer(j+2) for j,i in enumerate(one)]
    two,F=specialize(m,beta,True);q0=s.Matrix([point[i] for i in two]);fixed=source['chart_fixed_q'];cols=[i for i in range(45) if i!=fixed]
    D=s.zeros(135,45);tc=s.zeros(135,1)
    for r,f in enumerate(F):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q0[i]*q0[j];D[r,i]+=c*q0[j];D[r,j]+=c*q0[i]
            elif len(mon)==1:tc[r,0]+=c*q0[mon[0]]
        assert value==0
    Jac=D[:,cols].row_join(tc)
    assert Jac==s.Matrix(source['jacobian']).applyfunc(s.Rational)
    assert len(D.to_DM().rref()[1])==34 and len(Jac.to_DM().rref()[1])==35
    T=s.Matrix(source['tangent']).applyfunc(s.Rational)
    assert T.shape==(45,10) and Jac*T==s.zeros(135,10) and T[-1,:]==s.zeros(1,10)
    pivots=local['jacobian_pivot_columns'];free=[i for i in range(45) if i not in pivots]
    assert T[free,:]==s.eye(10)
    parameters=local['linear_parameters'];normal=local['remaining_normal_variables']
    assert parameters==[0,1,3] and normal==[2,4,5,6,7,8,9]
    assert transverse['linear_parameters']==parameters and transverse['normal_variables']==normal
    basis=s.Matrix(source['quadratic_obstruction_basis']).applyfunc(s.Rational)
    assert len(basis.to_DM().rref()[1])==5
    mons=[tuple(v) for v in source['quadratic_monomials']]
    hess=[]
    for row in basis.tolist():
        H=s.zeros(10)
        for c,(i,j) in zip(row,mons):H[i,j]+=c;H[j,i]+=c
        hess.extend(H.tolist())
    hess=s.Matrix(hess);assert len(hess.to_DM().rref()[1])==7
    radical=s.Matrix.hstack(*hess.nullspace());assert radical.shape==(10,3)
    # Construct the actual fixed-beta binary family, rather than identifying
    # its tangent space from the Hessian radical alone.
    z,uu,form=primitive_beta(m);form=form.subs(dict(zip(uu,[beta[i] for i in one])))
    a,b,d,r,h=s.symbols('a b d r h')
    v1=s.Matrix([a,1,0,b]);v2=s.Matrix([-b,0,1,d]);origin={a:0,b:0,d:0,r:0,h:0}
    section=s.Matrix([(v1.T*form.subs(dict(zip(z,v1)))*v2)[0],(v1.T*form.subs(dict(zip(z,v2)))*v2)[0]])
    assert section.subs(origin)==s.zeros(2,1)
    curve_jac=section.jacobian([a,b,d]).subs(origin)
    assert len(curve_jac.to_DM().rref()[1])==2
    curve_direction=curve_jac.nullspace()[0]
    x,y,zz,w=X[1:];ell1=x-a*y+b*zz;ell2=w-b*y-d*zz
    family=(ell1**2+r*ell1*ell2+(1+h)*ell2**2)*(v1*v2.T-v2*v1.T)
    base=s.zeros(4);base[1,2]=x*x+w*w;base[2,1]=-base[1,2]
    assert family.subs(origin)==base
    direction_columns=[]
    for variable in [a,b,d,r,h]:
        delta=family.diff(variable).subs(origin)
        P=X[0]**2*J;P[1:,1:]+=delta
        values=m.normalize(P)['parameters'];assert all(values[i]==0 for i in one)
        q=s.Matrix([values[i] for i in two]);q=q-q0*(q[fixed]/q0[fixed])
        direction_columns.append(q[cols,:].col_join(s.zeros(1,1)))
    plane=s.Matrix.hstack(*direction_columns[:3])*curve_direction
    family_tangent=s.Matrix.hstack(plane,direction_columns[3],direction_columns[4])
    assert Jac*family_tangent==s.zeros(135,3)
    coordinates=family_tangent[free,:]
    assert T*coordinates==family_tangent
    assert len(coordinates.to_DM().rref()[1])==3 and hess*coordinates==s.zeros(50,3)
    parameter_minor=coordinates[parameters,:].det();assert parameter_minor!=0
    assert len(radical.row_join(coordinates).to_DM().rref()[1])==3
    inputs=list(paths.values())+[Path(__file__),root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/primitive_jet_elliptic_curve_20261101.py',root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',root/'research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',root/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    return dict(id='S61_BINARY_FORMAL_RING_20261101',all_passed=True,field='QQ; completed local ring over C',point=source['point'],chart_fixed_q=fixed,leading_jacobian_rank=34,full_jacobian_rank=35,tangent_dimension=10,tangent_t_zero=True,quadratic_rank=5,common_hessian_radical_dimension=3,canonical_curve_derivative=enc(curve_jac),canonical_curve_tangent=list(map(str,curve_direction)),actual_binary_family_tangent_in_saved_coordinates=enc(coordinates),parameter_columns=parameters,parameter_minor=str(parameter_minor),same_fixed_beta_smooth_family_dimension=3,transverse_length=65,uniform_normal_length_lower_bound=65,completed_ring_finite_free_over='C[[v0,v1,v3]]',completed_ring_module_rank=65,reduced_completed_ring='C[[v0,v1,v3]]',nilradical_rank=64,nilradical_fourth_power_zero=True,associated_graded_free_ranks=[1,7,23,34],homogenizer_nonzero=True,homogenizer_square_zero=True,homogenizer_in_nilradical_cube=True,homogenizer_annihilates_nilradical=True,completed_infinity_quotient_module_rank=64,local_t_saturation_empty=True,proof_note='research/notes/S61_BINARY_FORMAL_RING_20261101.md',limitations=['No assertion that the rank65 algebra is the scalar extension of its special Artin fiber.','The multiplication among the seven normal generators is not completely classified.','Only the specified point and justified constant-rank open; no all-strata or full affine-fiber exhaustion.'],source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();r=verify()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','parameter_minor','completed_ring_module_rank','homogenizer_nonzero','homogenizer_square_zero','seconds']}))


if __name__=='__main__':main()
