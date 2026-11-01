#!/usr/bin/env python3
"""A flag/binary junction: finite free rank70 and a nonzero square-zero t."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_finite_field_local_length_20261101 import implicit_jets,truncated_lengths
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]


def verify():
    started=time.time();m=S61()
    flag_path=ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    flag=json.loads(flag_path.read_text());assert flag['all_passed']
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    H=s.Matrix(flag['flag']['beta_kernel']).applyfunc(s.Rational)
    coefficients=s.Matrix(range(2,17));coefficients[0]=0;beta_values=H*coefficients
    two,F=specialize(m,dict(zip(one,beta_values)),True)
    point=s.Matrix(flag['flag']['point61']).applyfunc(s.Rational)
    for i,v in zip(one,beta_values):point[i]=v
    q0=point[two,:];fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed]
    x,y,z,w=X[1:];Q0=s.zeros(4);Q0[1,3]=-x*x;Q0-=Q0.T
    def qcoordinates(Q):
        P=X[0]**2*J0;P[1:,1:]+=Q
        return s.Matrix([m.normalize(P)['parameters'][i] for i in two])
    assert qcoordinates(Q0)==q0 and q0[fixed]==-1
    J=s.zeros(135,45)
    for row,polynomial in enumerate(F):
        assert sum(c*s.prod(q0[i] for i in mon) for mon,c in polynomial.items() if len(mon)==2)==0
        for mon,c0 in polynomial.items():
            if len(mon)==2:
                i,j=mon
                if i!=fixed:J[row,cols.index(i)]+=c0*q0[j]
                if j!=fixed:J[row,cols.index(j)]+=c0*q0[i]
            elif len(mon)==1:J[row,44]+=c0*q0[mon[0]]
    T=s.Matrix.hstack(*J.nullspace());assert T.cols==11 and T[-1,:]==s.zeros(1,11)
    piv=list(J.to_DM().rref()[1]);free=[i for i in range(45) if i not in piv]
    assert T[free,:]==s.eye(11) and J*T==s.zeros(135,11)
    # The actual flag fiber is smooth3; this uses the entire flag incidence,
    # not just the tangent space to its b0=0 boundary.
    boundary=flag['frame_change']['boundary'][0]
    bd=s.Matrix(boundary['source_beta_differential']).applyfunc(s.Rational)
    assert bd.rank()==16
    BD=H.row_join(bd[:,-10:]);N=s.Matrix.hstack(*BD.nullspace())
    QD=s.Matrix(flag['flag']['q_orbit_differential']).applyfunc(s.Rational)
    CF=s.zeros(44,15).row_join(QD)*N;cp=list(CF.to_DM().rref()[1]);CF=CF[:,cp].col_join(s.zeros(1,len(cp)))
    assert CF.cols==3 and J*CF==s.zeros(135,3)
    flag_coordinates=CF[free,:];assert T*flag_coordinates==CF
    assert any(v!=0 for v in N[0,:])  # b0 changes along the same fixed-beta flag fiber.
    # The distinct actual binary fiber has the same square point as boundary.
    zz,uu,beta=primitive_beta(m);beta=beta.subs(dict(zip(uu,beta_values)))
    a,c,d,h,k=s.symbols('a c d h k');variables=[a,c,d,h,k];origin=dict.fromkeys(variables,0)
    v1=s.Matrix([a,1,c,0]);v2=s.Matrix([c,0,d,1])
    sec=s.Matrix([(v1.T*beta.subs(dict(zip(zz,v)))*v2)[0] for v in [v1,v2]])
    assert sec.subs(origin)==s.zeros(2,1)
    vertex_jac=sec.jacobian([a,c,d]).subs(origin);assert vertex_jac.rank()==2
    direction=vertex_jac.nullspace()[0]
    ell1=x-a*y-c*w;ell2=z-c*y-d*w
    family=-(ell1**2+h*ell1*ell2+k*ell2**2)*(v1*v2.T-v2*v1.T)
    derivatives=[]
    for v in variables:
        dq=qcoordinates(family.diff(v).subs(origin));dq-=q0*dq[fixed]/q0[fixed]
        derivatives.append(dq[cols,:].col_join(s.zeros(1,1)))
    CB=s.Matrix.hstack(s.Matrix.hstack(*derivatives[:3])*direction,derivatives[3],derivatives[4])
    assert J*CB==s.zeros(135,3)
    binary_coordinates=CB[free,:];assert T*binary_coordinates==CB
    projection=s.zeros(3,11);projection[0,0]=projection[1,1]=projection[1,5]=projection[2,2]=1
    flag_det=(projection*flag_coordinates).det();binary_det=(projection*binary_coordinates).det()
    assert flag_det==-s.Rational(46,7) and binary_det==1
    # Normal slice v0=v2=0, v1+v5=0 is transverse to BOTH actual components.
    normal=list(range(3,11));Tn=T[:,normal];Tn[:,normal.index(5)]-=T[:,1]
    assert projection*Tn[free,:]==s.zeros(3,8)
    # A nonsquare point on this same fixed-beta binary component reaches its
    # already established generic thickening of length65.
    Q1=s.zeros(4);Q1[1,3]=-(x*x+x*z+z*z);Q1-=Q1.T
    p1=point.copy()
    for i,v in zip(two,qcoordinates(Q1)):p1[i]=v
    comparison=calculate(m,True,p1)
    assert (comparison['leading_jacobian_rank'],comparison['homogenized_jacobian_rank'],
            comparison['quadratic_span_rank'],len(comparison['tangent_basis'][0])-comparison['common_hessian_rank'])==(34,35,5,3)
    jets,tjets,integral=implicit_jets(F,J,Tn,cols,maximum=6,prime=32003)
    lengths=truncated_lengths(jets,tjets,8,bounds=(5,6),prime=32003)
    assert lengths['5']['length']==lengths['6']['length']==70
    assert not lengths['6']['homogenizer_zero']
    paths=[Path(__file__),flag_path,ROOT/'code/s61_finite_field_local_length_20261101.py',
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        ROOT/'research/notes/ADIC_MODULAR_LENGTH_BOUND_20261101.md',
        ROOT/'research/notes/S61_FLAG_BINARY_CONNECTION_FORMAL_20261101.md']
    return dict(id='S61_FLAG_BINARY_CONNECTION_FORMAL_20261101',all_passed=True,
        verification_field='QQ geometry and full GF(32003) local-algebra truncations with integral lifting',
        point61=list(map(str,point)),jacobian=enc(J),tangent_basis=enc(T),free_rows=free,
        transverse_forms=['v0','v1+v5','v2'],transverse_tangent=enc(Tn),
        flag_tangent=enc(flag_coordinates),binary_tangent=enc(binary_coordinates),
        projection_determinants=dict(flag=str(flag_det),binary=str(binary_det)),
        fixed_beta_flag_b0_derivative=list(map(str,N[0,:])),binary_vertex_jacobian=enc(vertex_jac),
        binary_same_beta_comparison=comparison,integral_model=integral,truncated_algebras=lengths,
        characteristic_zero_transverse_length=70,component_generic_lengths=[5,65],
        completed_ring_finite_free_rank=70,completed_homogenizer_nonzero=True,
        completed_homogenizer_square_zero=True,relative_parameter_dimension=17,
        nonempty_ambient_open_on_flag_b0_divisor=True,entire_flag_b0_divisor_removed=False,
        proof_note='research/notes/S61_FLAG_BINARY_CONNECTION_FORMAL_20261101.md',
        limitations=['Specified completed ring and a nonempty ambient open on the18D flag boundary.',
            'The whole b0=0 boundary and the distinct b1=0 boundary remain unresolved.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','characteristic_zero_transverse_length','seconds']}))
