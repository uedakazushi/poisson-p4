#!/usr/bin/env python3
"""A rational seven-variable slice for the remaining deeper binary chart."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_binary_deeper_boundary_removal_20261101 import verify_point
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    rank=json.loads(rank_path.read_text());family=json.loads(family_path.read_text())
    assert rank['all_passed'] and family['all_passed']
    one=[3,4,11,12,13,24,25,26,33,34,35,36,47,48,49,50]
    u=s.symbols('u1:62');uu=[u[i] for i in one]
    a,b,c,d,e,f=[u[i] for i in [3,13,34,35,48,49]]
    K=5*c+f;L=5*b-e;Q=600*b*b-255*b*e-63*e*e+100*a*f
    A=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    reps=[A[:,i].jacobian(uu) for i in [0,1,3]]
    assert all(M*N==N*M for M in reps for N in reps)
    t=s.symbols('t0:3');N=sum((v*M for v,M in zip(t,reps)),s.zeros(16))
    assert N**4==s.zeros(16)
    exp=s.eye(16)+N+N**2/2+N**3/6
    value=exp*s.Matrix(uu)
    selected=[one.index(i) for i in [4,33,50]]
    corezero={u[11]:0,u[12]:0,f:-5*c}
    selected_shift=(value-s.Matrix(uu))[selected,:].subs(corezero).applyfunc(s.expand)
    mat=A[selected,[0,1,3]].subs(corezero)
    assert (selected_shift-mat*s.Matrix(t)).applyfunc(s.expand)==s.zeros(3,1)
    qparam={a:(600*b*b-255*b*e-63*e*e)/(500*c)}
    minor=s.factor(mat.det().subs(qparam))
    assert s.expand(minor-c*L*(20*b-37*e)/50)==0
    for i in [3,13,34,35,48,49]:
        assert s.expand((value[one.index(i)]-u[i]).subs(corezero))==0
    # The diagonal torus diag(x,x^-1,z,z^-1) has integral weights obtained
    # from the saved Lie generators. It acts by c -> z*c and L -> L/x.
    x,z=s.symbols('x z',nonzero=True)
    x_weights=[];z_weights=[]
    for i,v in enumerate(uu):
        wx=s.cancel(-2*A[i,4]/v)
        wz=s.cancel((A[i,2]+2*A[i,4])/v)
        assert wx.is_Integer and wz.is_Integer
        x_weights.append(int(wx));z_weights.append(int(wz))
    torus={v:x**wx*z**wz*v for v,wx,wz in zip(uu,x_weights,z_weights)}
    assert s.expand(c.subs(torus)-z*c)==0
    assert s.cancel(L.subs(torus,simultaneous=True)-L/x)==0
    normalized={v:s.cancel(image.subs({x:L,z:1/c},simultaneous=True)) for v,image in torus.items()}
    assert normalized[c]==1 and s.cancel(5*normalized[b]-normalized[e]-1)==0
    assert s.cancel(normalized[d]-d*L/c**2)==0 and s.cancel(normalized[a]-a*c/L**2)==0
    # The only missing unipotent divisor is H=20b-37e.  It is stable under
    # the full five-dimensional fixed-Q stabilizer.
    H=20*b-37*e
    cons=s.Matrix([u[11],u[12],K,Q,H])
    gb=s.groebner(list(cons),*uu,domain=s.QQ)
    assert all(gb.reduce(s.expand(v))[1]==0 for v in cons.jacobian(uu)*A)
    sample=dict(family['sample'])
    ev=s.Rational(4,33);bv=(ev+1)/5;av=3*(1-2*ev)*(15*ev+8)/500
    sample.update({k:str(v) for k,v in zip(['u4','u14','u35','u36','u49','u50'],[av,bv,1,12,ev,-5])})
    point={s.Symbol(k):s.Rational(v) for k,v in sample.items()}
    assert cons.subs(point)==s.zeros(5,1)
    fullaction=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    normal=(cons.jacobian(uu)*fullaction).subs(point)
    cp=list(normal.to_DM().rref()[1]);assert len(cp)==5
    beta_minor=normal[:,cp].det();assert beta_minor
    witness=verify_point(sample)
    E,D,P,Qv,Rv,Hv,Kv=s.symbols('e d p q r h k')
    slice_beta=dict.fromkeys(uu,s.S.Zero)
    slice_beta.update({a:3*(1-2*E)*(15*E+8)/500,b:(E+1)/5,c:1,d:D,e:E,f:-5,
        u[24]:P,u[25]:Qv,u[26]:Rv,u[36]:Hv,u[47]:Kv})
    assert all(s.expand(v.subs(slice_beta))==0 for v in [u[11],u[12],K,Q,u[4],u[33],u[50],c-1,L-1])
    assert s.expand(H.subs(slice_beta)-(4-33*E))==0
    # Five normal conditions give a transverse rational slice to the fixed-Q
    # stabilizer wherever the stated determinant is nonzero.
    normals=s.Matrix([c-1,L-1,u[4],u[33],u[50]])
    slice_det=s.factor((normals.jacobian(uu)*A).subs(slice_beta).det())
    assert slice_det!=0 and s.cancel(slice_det/(4-33*E)).is_Rational
    paths=[Path(__file__),rank_path,family_path,
        ROOT/'code/s61_binary_deeper_boundary_removal_20261101.py',
        ROOT/'reports/S61_BINARY_T2_STRATUM_REMOVAL_20261101.json',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        ROOT/'research/notes/S61_BINARY_DEEPER_SEVEN_SLICE_20261101.md']
    return dict(id='S61_BINARY_DEEPER_SEVEN_SLICE_20261101',all_passed=True,
        unipotent_generators=[0,1,3],unipotent_target_coordinates=['u5','u34','u51'],
        unipotent_normal_matrix=[[str(v) for v in row] for row in mat.tolist()],
        unipotent_normal_determinant=str(minor),nilpotence_order_bound=4,
        torus_matrix=['L','1/L','1/c','c'],torus_x_weights=x_weights,torus_z_weights=z_weights,
        normalized_beta={str(k):str(v) for k,v in slice_beta.items()},
        normalized_variables=list(map(str,[E,D,P,Qv,Rv,Hv,Kv])),
        slice_normal_determinant=str(slice_det),
        missing_divisor_source=list(map(str,cons)),missing_divisor_irreducible_on_c_nonzero=True,
        missing_divisor_beta_minor=str(beta_minor),missing_divisor_beta_minor_columns=cp,
        missing_divisor_joint_dimension=16,missing_divisor_entire_general_beta_removed=True,
        missing_divisor_witness=witness,
        proof_note='research/notes/S61_BINARY_DEEPER_SEVEN_SLICE_20261101.md',
        limitations=['This is a rational coordinate reduction and one16D boundary removal.',
          'It does not prove the seven-variable open has length65 at every point.',
          'No full affine-fiber exhaustion is asserted.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();data=verify()
    with args.output.open('x') as stream:json.dump(data,stream,indent=2);stream.write('\n')
    print(json.dumps({k:data[k] for k in ['all_passed','slice_normal_determinant','missing_divisor_beta_minor','seconds']}))
