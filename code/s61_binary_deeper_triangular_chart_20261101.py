#!/usr/bin/env python3
"""Exact triangular linear elimination and two removed chart artifacts."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_binary_deeper_boundary_removal_20261101 import verify_point
ROOT=Path(__file__).resolve().parents[1]


def build():
    start=time.time()
    source=json.loads((ROOT/'reports/S61_BINARY_DEEPER_SEVEN_SLICE_20261101.json').read_text())
    old={'variables':source['normalized_variables'],'beta':source['normalized_beta']}
    m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in old['beta'].items()};two,F=specialize(m,beta,True)
    fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45)
    for row,f in enumerate(F):
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon
                if i!=fixed and j==fixed:J[row,cols.index(i)]+=c
                if j!=fixed and i==fixed:J[row,cols.index(j)]+=c
            elif len(mon)==1 and mon[0]==fixed:J[row,44]+=c
    old.update(fixed=fixed,cols=cols,J=[[str(v) for v in row] for row in J.tolist()])
    data=json.load(open(ROOT/'reports/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.json'))
    K=s.Matrix(data['constant_first70_kernel']).applyfunc(s.sympify);fixed=old['fixed'];cols=old['cols'];rr=data['residual_block_rows'];cc=data['residual_block_columns'];qremove=next(j for j in range(29) if K[fixed,j]);assert K[fixed,qremove]==1
    piv=list(J[:70,:44].to_DM().rref()[1]);pr=list(J[:70,piv].T.to_DM().rref()[1]);assert len(piv)==16
    P=s.eye(45)[:,piv];ker=K[cols,[j for j in range(29) if j!=qremove]].col_join(s.zeros(1,28));W=P.row_join(ker).row_join(s.eye(45)[:,44]);assert W.det()!=0
    kcols=[j for j in range(29) if j!=qremove];Acol=[0,3,7,9];Bcol=[0,1,2,3,4,5,7];Brow=[0,1,2,3,5,6,10];Crow=[0,1,2,5,8]
    kp=[cc[0][i] for i in Acol]+[cc[1][i] for i in Bcol]+[cc[2][i] for i in range(5)]
    pivotcols=list(range(16))+[44]+[16+kcols.index(i) for i in kp]
    pivotrows=pr+[15]+[70+i for i in rr[0]]+[70+rr[1][i] for i in Brow]+[70+rr[2][i] for i in Crow]
    JJ=(J.to_DM()*W.to_DM()).to_Matrix();block=JJ[pivotrows,pivotcols].to_DM().to_field();det=block.det();print('det',s.factor(det.as_expr()),time.time()-start,flush=True)
    inv=block.inv();rest=[i for i in range(45) if i not in pivotcols];Tc=s.zeros(45,12);Tc[rest[0],0]=1
    rhs=JJ[pivotrows,rest].to_DM().to_field();sol=-(inv*rhs).to_Matrix()
    for i,j in enumerate(rest):Tc[j,i]=1
    for i,j in enumerate(pivotcols):Tc[j,:]=sol[i,:]
    T=(W.to_DM()*Tc.to_DM()).to_Matrix();assert (J.to_DM()*T.to_DM()).is_zero_matrix
    assert T[-1,:]==s.zeros(1,12)
    Uc=s.zeros(45,33);ii=inv.to_Matrix()
    for i,j in enumerate(pivotcols):Uc[j,:]=ii[i,:]
    U=(W.to_DM()*Uc.to_DM()).to_Matrix();assert (J[pivotrows,:].to_DM()*U.to_DM()).to_Matrix().applyfunc(s.cancel)==s.eye(33)
    # Vertex/conic family in actual q chart.
    m=S61();z,uu,bf=primitive_beta(m);bf=bf.subs({s.Symbol(k):s.sympify(v) for k,v in old['beta'].items()});a,c,d,h,k=s.symbols('A C D H K');v1=s.Matrix([a,1,0,c]);v2=s.Matrix([-c,0,1,d]);origin={a:0,c:0,d:0,h:0,k:0}
    sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]]);cj=sec.jacobian([a,c,d]).subs(origin);direction=s.Matrix(cj[0,:]).cross(s.Matrix(cj[1,:]));direction=s.Matrix(direction).reshape(3,1)
    x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell1*ell2+h*ell1**2+k*ell2**2)*(v1*v2.T-v2*v1.T);ds=[];two=[i for i,weight in enumerate(m.reduced['parameter_weights']) if weight==2];q0=s.zeros(45,1);q0[fixed]=1
    for v in [a,c,d,h,k]:
     P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);q=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);q-=q0*q[fixed];ds.append(q[cols,:].col_join(s.zeros(1,1)))
    family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);fc=(W.inv()*family)[rest,:];params=list(fc.T.to_DM().rref()[1]);print('params',params,'minor',s.factor(fc[params,:].det()),time.time()-start,flush=True)
    assert (J.to_DM()*family.to_DM()).is_zero_matrix
    assert (T.to_DM()*fc.to_DM()-family.to_DM()).is_zero_matrix
    assert len(params)==3
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    out=old|dict(T=enc(T),U=enc(U),pivot_rows=pivotrows,pivot_columns=pivotcols,chart_basis=enc(W),free=rest,parameters=params,normal=[i for i in range(12) if i not in params],curve_jacobian=enc(cj),family_coordinates=enc(fc),family_parameter_minor=str(s.factor(fc[params,:].det())),jacobian_minor=str(s.factor(det.as_expr())),seconds=time.time()-start)
    return out


def verify():
    start=time.time();linear=build()
    source_path=ROOT/'reports/S61_BINARY_DEEPER_SEVEN_SLICE_20261101.json'
    source=json.loads(source_path.read_text())
    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    family=json.loads(family_path.read_text())
    rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    rank=json.loads(rank_path.read_text())
    one=[3,4,11,12,13,24,25,26,33,34,35,36,47,48,49,50]
    u=s.symbols('u1:62');uu=[u[i] for i in one]
    a,b,c,d,e,f=[u[i] for i in [3,13,34,35,48,49]]
    K=5*c+f;Q=600*b*b-255*b*e-63*e*e+100*a*f
    action=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    stab=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    cases=[]
    for label,ev,H in [('A_C_pivot',s.Rational(-8,15),40*b+7*e),('binary_projection',s.Rational(7,24),35*b-31*e)]:
        cons=s.Matrix([u[11],u[12],K,Q,H]);gb=s.groebner(list(cons),*uu,domain=s.QQ)
        assert all(gb.reduce(s.expand(v))[1]==0 for v in cons.jacobian(uu)*stab)
        sample=dict(family['sample']);bv=(ev+1)/5;av=3*(1-2*ev)*(15*ev+8)/500
        sample.update({k:str(v) for k,v in zip(['u4','u14','u35','u36','u49','u50'],[av,bv,1,12,ev,-5])})
        point={s.Symbol(k):s.Rational(v) for k,v in sample.items()}
        normal=(cons.jacobian(uu)*action).subs(point);cp=list(normal.to_DM().rref()[1]);assert len(cp)==5
        det=normal[:,cp].det();assert det
        witness=verify_point(sample)
        cases.append(dict(label=label,source_equations=list(map(str,cons)),source_dimension=11,
            beta_minor=str(det),beta_minor_columns=cp,joint_dimension=16,
            entire_general_beta_stratum_removed=True,witness=witness))
    paths=[Path(__file__),source_path,family_path,rank_path,
        ROOT/'reports/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.json',
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'code/s61_binary_deeper_boundary_removal_20261101.py',
        ROOT/'research/notes/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.md']
    return dict(id='S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101',all_passed=True,
        linear=linear,removed_artifact_divisors=cases,
        proof_note='research/notes/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.md',
        limitations=['The remaining open is only reduced to a controlled linear chart.',
            'No universal higher-order colength65 or full affine exhaustion is asserted.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();data=verify()
    with args.output.open('x') as stream:json.dump(data,stream,indent=2);stream.write('\n')
    print(json.dumps({k:data[k] for k in ['all_passed','seconds']}))
