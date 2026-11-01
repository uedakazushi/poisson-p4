#!/usr/bin/env python3
"""Independent fixed-beta binary threefold and pointwise135 tangent audit."""
import argparse,hashlib,json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]


def run():
    stems=['S61_BINARY_FORMAL_RING','S61_BINARY_QUADRATIC_OBSTRUCTION','S61_BINARY_TRANSVERSE_ALGEBRA',
           'S61_BINARY_UNIFORM_RANK','S61_BINARY_UNIFORM_QUADRATIC','S61_BINARY_UNIFORM_CUBIC']
    paths=[ROOT/'reports'/f'{stem}_20261101.json' for stem in stems]
    reports=[json.loads(p.read_text()) for p in paths];assert all(r['all_passed'] for r in reports)
    main,quadratic,transverse,rank_report,quad_report,cubic=reports
    assert main['point']==quadratic['point']==transverse['point']
    assert transverse['transverse_algebra_length']==65 and transverse['transverse_infinity_quotient_length']==64
    assert transverse['maximal_ideal_fourth_power_zero'] and transverse['homogenizing_class_in_cubic_socle']
    assert transverse['homogenizing_class_nonzero'] and cubic['normal_length_lower_bound']==65
    assert cubic['cubic_relation_rank_upper_bound']==50 and 44 not in cubic['lift_support_columns']
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];index={i:j for j,i in enumerate(two)}
    point=list(map(s.Rational,main['point']));q=s.Matrix([point[i] for i in two])
    fixed=main['chart_fixed_q'];cols=[i for i in range(45) if i!=fixed]
    source=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    polynomials=[]
    for eq in source:
        f={}
        for mon,num,den in eq:
            qm=tuple(index[i] for i in mon if i in index)
            c=s.Rational(num,den)*s.prod(point[i] for i in mon if i in one)
            f[qm]=f.get(qm,s.S.Zero)+c
        polynomials.append(f)
    A=s.zeros(135,45);tc=s.zeros(135,1)
    for row,f in enumerate(polynomials):
        value=s.S.Zero
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q[i]*q[j];A[row,i]+=c*q[j];A[row,j]+=c*q[i]
            elif len(mon)==1:tc[row]+=c*q[mon[0]]
        assert value==0
    Jac=A[:,cols].row_join(tc)
    assert Jac==s.Matrix(quadratic['jacobian']).applyfunc(s.Rational)
    assert A.rank()==34 and Jac.rank()==35
    T=s.Matrix(quadratic['tangent']).applyfunc(s.Rational)
    assert Jac*T==s.zeros(135,10) and T.rank()==10 and T[-1,:]==s.zeros(1,10)
    full=s.zeros(46,10)
    for j,i in enumerate(cols):full[i,:]=T[j,:]
    full[45,:]=T[-1,:]
    mons=list(combinations_with_replacement(range(10),2));H=s.zeros(135,len(mons))
    for row,f in enumerate(polynomials):
        for mon,c in f.items():
            if fixed in mon:continue
            i,j=[*mon,*([45]*(2-len(mon)))]
            for k,(a,b) in enumerate(mons):H[row,k]+=c*(full[i,a]*full[j,b]+(full[i,b]*full[j,a] if a!=b else 0))
    left=s.Matrix.hstack(*Jac.T.nullspace()).T;projected=left*H
    assert projected.rank()==5
    selected=list(projected.T.rref()[1]);forms=projected[selected,:];hess=[]
    for row in forms.tolist():
        B=s.zeros(10)
        for c,(i,j) in zip(row,mons):B[i,j]+=c;B[j,i]+=c
        hess.append(B)
    stacked=s.Matrix.vstack(*hess);assert stacked.rank()==7

    # Recover the old closed first-jet form directly from the reconstructed
    # projective bivector. Its restriction to isotropic planes is the same
    # as the primitive projection, so no primitive-beta helper is needed.
    z=s.symbols('z0:4');B=m.reconstruct(point);subst=dict(zip(X,(1,*z)))
    P=s.Matrix(4,4,lambda i,j:s.expand((B[i+1,j+1]-X[i+1]*B[0,j+1]+X[j+1]*B[0,i+1]).subs(subst)))
    zero=dict.fromkeys(z,0);linear=s.zeros(4)
    for zi in z:linear+=P.diff(zi).subs(zero)*zi
    form=(-J[1:,1:]*linear*J[1:,1:]).applyfunc(s.expand)
    a,b,d,r,h=s.symbols('a b d r h');v1=s.Matrix([a,1,0,b]);v2=s.Matrix([-b,0,1,d])
    assert (v1.T*J[1:,1:]*v2)[0]==0
    section=s.Matrix([(v1.T*form.subs(dict(zip(z,v1)))*v2)[0],
                      (v1.T*form.subs(dict(zip(z,v2)))*v2)[0]])
    origin={a:0,b:0,d:0,r:0,h:0}
    assert section.subs(origin)==s.zeros(2,1)
    derivative=section.jacobian([a,b,d]).subs(origin)
    assert derivative==s.Matrix([[22,-21,2],[12,-27,3]])
    direction=s.Matrix([s.Rational(1,38),s.Rational(7,57),1])
    assert derivative*direction==s.zeros(2,1) and derivative.rank()==2
    x,y,zz,w=X[1:];ell1=x-a*y+b*zz;ell2=w-b*y-d*zz
    family=(ell1**2+r*ell1*ell2+(1+h)*ell2**2)*(v1*v2.T-v2*v1.T)
    directions=[]
    for variable in (a,b,d,r,h):
        delta=family.diff(variable).subs(origin);biv=X[0]**2*J;biv[1:,1:]+=delta
        values=m.normalize(biv)['parameters'];assert all(values[i]==0 for i in one)
        dq=s.Matrix([values[i] for i in two]);dq-=q*dq[fixed]/q[fixed]
        directions.append(dq[cols,:].col_join(s.zeros(1,1)))
    actual=s.Matrix.hstack(s.Matrix.hstack(*directions[:3])*direction,*directions[3:])
    rows=list(T.T.rref()[1]);coordinates=T[rows,:].inv()*actual[rows,:]
    assert T*coordinates==actual and Jac*actual==s.zeros(135,3)
    assert coordinates==s.Matrix(main['actual_binary_family_tangent_in_saved_coordinates']).applyfunc(s.Rational)
    assert coordinates.rank()==3 and stacked*coordinates==s.zeros(50,3)
    minor=coordinates[[0,1,3],:].det();assert minor==-s.Rational(1,19)
    paths += [Path(__file__),ROOT/'code/s61_binary_formal_ring_20261101.py',
              ROOT/'research/notes/S61_BINARY_FORMAL_RING_20261101.md',
              ROOT/'research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',
              ROOT/'research/notes/S61_BINARY_TRANSVERSE_ALGEBRA_20261101.md',
              m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_BINARY_FORMAL_RING_REVIEW_20261101',all_passed=True,
        leading_rank=34,full_rank=35,tangent_dimension=10,quadratic_rank=5,common_radical_dimension=3,
        independent_curve_derivative=[[str(v) for v in row] for row in derivative.tolist()],
        actual_family_parameter_minor=str(minor),same_fixed_beta_smooth_family_verified=True,
        accepted_module_rank=65,accepted_graded_free_ranks=[1,7,23,34],accepted_infinity_quotient_rank=64,
        scope='Independent point Jacobian/Hessians and direct first-jet/family reconstruction; universal rank/cubic and transverse Artin certificates are explicit existing dependencies.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=run()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'output':str(args.output)}))
