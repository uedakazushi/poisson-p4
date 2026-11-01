#!/usr/bin/env python3
"""A universal nine-parameter actual double line at the binary/mixed junction."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J,COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time();model=S61();u=s.symbols('u1:62')
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    b={i:u[i] for i in one}
    b.update({3:0,4:0,11:0,12:0,13:s.Rational(3,5)*u[48],49:-5*u[34],50:5*u[25]})
    two,F=specialize(model,b,True);fixed=two.index(29);chart=[i for i in range(45) if i!=fixed]
    q0=s.zeros(45,1);q0[fixed]=1
    x,y,z,w=X[1:];P=X[0]**2*J0;P[2,3]+=z*z;P[3,2]-=z*z;P[2,4]-=2*z*w;P[4,2]+=2*z*w
    delta=s.Matrix([model.normalize(P)['parameters'][i] for i in two]);assert delta[fixed]==0
    nil=s.zeros(45,1)
    for i,val in {0:5*u[48],19:8*u[25],24:s.Rational(16,3)*u[34],27:s.Rational(16,5)*u[48],34:u[48]}.items():nil[chart[i]]=val
    for f in F:
        for v,ww in [(q0,delta),(q0,nil),(delta,nil)]:
            assert s.expand(sum(c*(v[mon[0]]*ww[mon[1]]+ww[mon[0]]*v[mon[1]])
                for mon,c in f.items() if len(mon)==2))==0
        for v in [q0,delta]:
            assert s.expand(sum(c*v[mon[0]]*v[mon[1]] for mon,c in f.items() if len(mon)==2))==0
    assert delta[chart[32]]==-2 and nil[chart[32]]==0
    assert delta[chart[34]]==0 and nil[chart[34]]==u[48]

    # Normalize the proved fixed-beta point by an actual finite symplectic pullback.
    point_path=ROOT/'reports/S61_BINARY_T2_CONNECTION_20261101.json'
    sample={s.Symbol(k):s.Rational(v) for k,v in json.loads(point_path.read_text())['point']['sample'].items()}
    zz,uu,beta=primitive_beta(model)
    tensor=s.Matrix([beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=tensor.jacobian(uu)
    rows=list(basis.T.to_DM().rref()[1]);left=basis[rows,:].inv()
    h=s.Integer(23);r=s.Integer(15)
    shear=s.Matrix([[1,0,0,0],[0,1,0,-h],[h,0,1,r],[0,0,0,1]])
    assert shear.T*J*shear==J
    pulled=shear.T*beta.subs(sample).subs(dict(zip(zz,shear*s.Matrix(zz))),simultaneous=True)*shear
    vector=s.Matrix([pulled[i,j].diff(zz[k]) for (i,j),k in COORDS]);values=left*vector[rows,:]
    assert basis*values==vector
    point=dict(zip(uu,values));constraints=s.Matrix([u[3],u[4],u[11],u[12],
        5*u[13]-3*u[48],5*u[34]+u[49],u[50]-5*u[25]])
    assert constraints.subs(point)==s.zeros(7,1) and point[u[48]]==1
    action=s.Matrix.hstack(*[(left*infinitesimal_action(tensor,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()])
    normal=(constraints.jacobian(uu)*action).subs(point);columns=list(normal.to_DM().rref()[1])
    assert len(columns)==7
    determinant=normal[:,columns].det();assert determinant!=0
    # The shear also fixes Q0 and sends the old tail's ell to z.
    V=s.Matrix(zz);Q=s.zeros(4);Q[1,2]=zz[0]*zz[3];Q-=Q.T
    assert shear.inv()*Q.subs(dict(zip(zz,shear*V)),simultaneous=True)*shear.inv().T==Q
    ell=zz[2]-15*zz[3]-23*zz[0]
    assert s.expand(ell.subs(dict(zip(zz,shear*V)),simultaneous=True))==zz[2]
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    paths=[Path(__file__),point_path,ROOT/'code/s61_model.py',
        ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_BINARY_T2_RELATIVE_DOUBLE_LINE_20261101.md']
    return dict(id='S61_BINARY_T2_RELATIVE_DOUBLE_LINE_20261101',all_passed=True,
        source_equations=list(map(str,constraints)),source_dimension=9,open_condition='u49 != 0',
        tail_delta=list(map(str,delta)),nil_direction=list(map(str,nil)),
        all_135_actual_double_line_identities=True,flat_relative_double_line=True,
        normalized_fixed_beta={str(k):str(v) for k,v in point.items()},normalizing_shear=enc(shear),
        beta_normal_matrix=enc(normal),beta_minor_columns=columns,beta_minor=str(determinant),
        proof_note='research/notes/S61_BINARY_T2_RELATIVE_DOUBLE_LINE_20261101.md',
        limitations=['This constructs a flat closed quotient D of the relative full scheme.',
            'Relative torsion-freeness of its kernel and ambient t-nilpotence are not proved here.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,beta_minor=report['beta_minor'],seconds=report['seconds'])))
