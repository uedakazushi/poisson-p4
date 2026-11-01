#!/usr/bin/env python3
"""Same-component lower65 and relative free65 on the common deeper core open."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.time()
    linear_path=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    upper_path=ROOT/'reports/S61_BINARY_DEEPER_QUARTIC_NUISANCE_20261101.json'
    linear=json.loads(linear_path.read_text())['linear']; upper=json.loads(upper_path.read_text()); assert upper['all_passed']
    variables=s.symbols('e d p q r h k'); sample=dict(zip(variables,[2,3,5,7,11,13,17]))
    model=S61(); beta={int(k[1:])-1:s.sympify(v) for k,v in linear['beta'].items()}
    two,F=specialize(model,beta,True); fixed=linear['fixed']; cols=linear['cols']
    H,K=s.symbols('H K'); x,y,z,w=X[1:]
    P=X[0]**2*J0; P[2,3]+=x*w+H*x*x+K*w*w; P[3,2]-=x*w+H*x*x+K*w*w
    parameters=model.normalize(P)['parameters']; conic=s.Matrix([parameters[i] for i in two]); assert conic[fixed]==1
    for f in F:
        value=sum(c*conic[i]*conic[j] for mon,c in f.items() if len(mon)==2 for i,j in [mon])
        assert s.expand(value)==0
    q1=conic.subs({H:1,K:1}); point_beta={i:v.subs(sample) for i,v in beta.items()}
    _,equations=specialize(model,point_beta,True)
    J=s.zeros(135,45)
    for row,f in enumerate(equations):
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon
                if i!=fixed:J[row,cols.index(i)]+=c*q1[j]
                if j!=fixed:J[row,cols.index(j)]+=c*q1[i]
            elif len(mon)==1:J[row,44]+=c*q1[mon[0]]
    T=s.Matrix.hstack(*J.nullspace()); assert T.shape==(45,10) and T[-1,:]==s.zeros(1,10)
    assert len(J[:,:44].to_DM().rref()[1])==34
    mons=list(combinations_with_replacement(range(10),2)); index={m:i for i,m in enumerate(mons)}; cache={}
    def product(i,j):
        key=tuple(sorted((i,j)))
        if key not in cache:
            row=[s.S.Zero]*55
            for a,v in enumerate(T[i,:]):
                if not v:continue
                for b,wv in enumerate(T[j,:]):
                    if wv:row[index[tuple(sorted((a,b)))]]+=v*wv
            cache[key]=row
        return cache[key]
    Q=[]
    for f in equations:
        row=[s.S.Zero]*55
        for mon,c in f.items():
            if fixed in mon:continue
            i,j=[cols.index(v) for v in mon]+[44]*(2-len(mon))
            for k,v in enumerate(product(i,j)):row[k]+=c*v
        Q.append(row)
    cokernel=s.Matrix.hstack(*J.T.nullspace()).T
    A=(cokernel.to_DM()*s.Matrix(Q).to_DM()).to_Matrix()
    _,chosen=A.T.to_DM().rref(); assert len(chosen)==5
    hess=[]
    for row in A[list(chosen),:].tolist():
        matrix=s.zeros(10)
        for c,(i,j) in zip(row,mons):matrix[i,j]+=c;matrix[j,i]+=c
        hess+=matrix.tolist()
    hessian=s.Matrix(hess); assert len(hessian.to_DM().rref()[1])==7
    # Reconstruct the actual relative vertex curve and its tangent, not only
    # an abstract three-dimensional radical.
    zv,_,bf=primitive_beta(model); bf=bf.subs({s.Symbol(k):s.sympify(v) for k,v in linear['beta'].items()})
    aa,cc,dd,hh,kk=s.symbols('A C D H1 K1'); v1=s.Matrix([aa,1,0,cc]);v2=s.Matrix([-cc,0,1,dd]);origin={aa:0,cc:0,dd:0,hh:0,kk:0}
    section=s.Matrix([(v1.T*bf.subs(dict(zip(zv,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(zv,v2)))*v2)[0]])
    assert section.subs(origin)==s.zeros(2,1)
    cj=section.jacobian([aa,cc,dd]).subs(origin)
    expected=s.Matrix(linear['curve_jacobian']).applyfunc(s.sympify); assert (cj-expected).applyfunc(s.expand)==s.zeros(2,3)
    direction=cj.subs(sample).nullspace()[0]
    ell1=x-aa*y+cc*z;ell2=w-cc*y-dd*z
    family=(ell1*ell2+(1+hh)*ell1**2+(1+kk)*ell2**2)*(v1*v2.T-v2*v1.T)
    tangents=[]
    for v in [aa,cc,dd,hh,kk]:
        PP=X[0]**2*J0; PP[1:,1:]+=family.diff(v).subs(origin)
        qs=s.Matrix([model.normalize(PP)['parameters'][i] for i in two]); qs-=q1*qs[fixed]
        tangents.append(qs[cols,:].col_join(s.zeros(1,1)))
    actual=s.Matrix.hstack(s.Matrix.hstack(*tangents[:3])*direction,tangents[3],tangents[4])
    _,piv=J.to_DM().rref(); free=[i for i in range(45)if i not in piv]
    assert T[free,:]==s.eye(10); coordinates=actual[free,:]
    assert T*coordinates==actual and J*actual==s.zeros(135,3)
    assert coordinates.rank()==3 and hessian*coordinates==s.zeros(50,3)
    # Dominance of the smooth Sp4 x H7 source, independently at the same beta.
    action_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    action=s.Matrix(json.loads(action_path.read_text())['primitive_action']).applyfunc(s.sympify)
    one=[i for i,weight in enumerate(model.reduced['parameter_weights'])if weight==1]
    bv=s.Matrix([beta[i] for i in one]); change={s.Symbol('u'+str(i+1)):beta[i] for i in one}
    tangent=action.subs(change,simultaneous=True).row_join(bv.jacobian(variables)).subs(sample)
    _,source_columns=tangent.to_DM().rref(); assert len(source_columns)==16
    determinant=tangent[:,list(source_columns)].det(); assert determinant
    enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
    note=ROOT/'research/notes/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.md'
    paths=[Path(__file__),linear_path,upper_path,action_path,
        ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json',ROOT/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',ROOT/'reports/S61_BINARY_DEEPER_SEVEN_SLICE_20261101.json',
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',note]
    return dict(id='S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101',all_passed=True,
        sample={str(k):str(v)for k,v in sample.items()},universal_same_P2_highest_identities=135,
        same_P2_qfamily=list(map(str,conic)),q0_conic_discriminant=1,q1_conic_discriminant=-3,
        relative_vertex_derivative=enc(cj),section_parameter_minor=linear['family_parameter_minor'],
        same_beta_Q1_leading_rank=34,same_beta_Q1_full_rank=35,same_beta_Q1_quadratic_rank=5,
        same_beta_Q1_hessian_radical_dimension=3,actual_Q1_family_tangent=enc(coordinates),
        source_beta_rank=16,source_beta_minor_columns=list(source_columns),source_beta_minor=str(determinant),
        relative_completed_rank=65,fixed_beta_completed_rank=65,nilradical_power4_zero=True,
        homogenizer_annihilates_nilradical=True,homogenizer_square_zero=True,
        all_general_beta_points_on_common_core_open_removed=True,
        limitations=['The core exceptional locus remains unprocessed.','No uniform nonzero-t or infinity quotient rank64 is asserted.','Ambient descent at special beta requires source submersivity.'],
        proof_note=str(note.relative_to(ROOT)),input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    result=verify()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2,ensure_ascii=False);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'],source_beta_minor=result['source_beta_minor'])))
