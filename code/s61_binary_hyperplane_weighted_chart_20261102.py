#!/usr/bin/env python3
"""Exact eight-normal-variable positive-weight chart on the K=0 hypersurface."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,BASIS,X
from s61_model import J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
root=Path(__file__).resolve().parents[1]
NOTE=root/'research/notes/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.md'

def verify():
    start=time.time()
    data=json.loads((root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    model=S61();a,e,d=s.symbols('a e d')
    beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()}
    for name in ['p','q','r','h','k']:
        beta={i:v.subs(s.Symbol(name),0) for i,v in beta.items()}
    beta[3]=a
    two,F=specialize(model,beta,True)
    fixed=data['fixed'];cols=data['cols'];J=s.zeros(135,45);H=[]
    cw=json.loads((root/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json').read_text())['local_coordinate_weights']
    for row,f in enumerate(F):
        terms=[]
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon
                if i!=fixed and j==fixed:J[row,cols.index(i)]+=c
                if j!=fixed and i==fixed:J[row,cols.index(j)]+=c
            elif len(mon)==1 and mon[0]==fixed:J[row,44]+=c
            if fixed not in mon:
                ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
                if c:terms.append((ij,c))
        H.append(terms)
    rw=[]
    for i in range(135):
        weights={cw[j] for j,v in enumerate(J[i,:]) if v}
        weights.update(sum(cw[j] for j in ij) for ij,c in H[i] if c)
        if len(weights)>1:print('BADWEIGHT',i,weights,flush=True)
        rw.append(list(weights))
    assert all(len(v)==1 for v in rw)
    print('ALL_WEIGHTED',time.time()-start,flush=True)
    W=s.Matrix(data['chart_basis']).applyfunc(s.sympify)
    JJ=(J.to_DM()*W.to_DM()).to_Matrix();pr=data['pivot_rows'];pc=data['pivot_columns']
    A=JJ[pr,pc].to_DM().to_field()
    print('BEFORE_INV',A.shape,time.time()-start,flush=True)
    Ai=A.inv();print('INV',time.time()-start,flush=True)
    rest=data['free'];T=s.zeros(45,12)
    for i,j in enumerate(rest):T[j,i]=1
    sol=-(Ai*JJ[pr,rest].to_DM().to_field()).to_Matrix()
    for i,j in enumerate(pc):T[j,:]=sol[i,:]
    TT=(W.to_DM()*T.to_DM()).to_Matrix()
    R=(J.to_DM()*TT.to_DM()).to_Matrix().applyfunc(s.cancel)
    print('RESIDUAL_RANK',R.to_DM().rref()[1],time.time()-start,flush=True)
    nonzero=[(i,j,s.factor(v)) for i in range(R.rows) for j,v in enumerate(R[i,:]) if v]
    print('RESIDUAL_NONZERO',len(nonzero),nonzero[:12],flush=True)
    det=s.factor(A.det().as_expr());print('DETERMINANT',det,flush=True)
    out=dict(all_weights_preserved=True,selected33_minor=str(det),residual_nonzero=[(i,j,str(v)) for i,j,v in nonzero],seconds=time.time()-start)
    pr34=pr+[103];pc34=pc+[rest[10]];rest34=[i for i in range(45) if i not in pc34]
    A34=JJ[pr34,pc34].to_DM().to_field();Ai34=A34.inv()
    T34=s.zeros(45,11);U34=s.zeros(45,34)
    for i,j in enumerate(rest34):T34[j,i]=1
    sol34=-(Ai34*JJ[pr34,rest34].to_DM().to_field()).to_Matrix()
    for i,j in enumerate(pc34):T34[j,:]=sol34[i,:];U34[j,:]=Ai34.to_Matrix()[i,:]
    T34=(W.to_DM()*T34.to_DM()).to_Matrix()
    U34=(W.to_DM()*U34.to_DM()).to_Matrix()
    assert (J.to_DM()*T34.to_DM()).is_zero_matrix
    assert (J[pr34,:].to_DM()*U34.to_DM()).to_Matrix()==s.eye(34)
    cuts=[rest.index(i) for i in rest34 if rest.index(i) in data['parameters']]
    normal34=[i for i,j in enumerate(rest34) if rest.index(j) not in data['parameters']]
    nw34=[json.loads((root/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json').read_text())['normal_weights'][data['normal'].index(rest.index(rest34[i]))] for i in normal34]
    assert len(normal34)==8
    assert all(not T34[i,j] or cw[i]==nw34[k] for i in range(45) for k,j in enumerate(normal34))
    assert all(not U34[i,j] or cw[i]==rw[pr34[j]][0] for i in range(45) for j in range(34))
    det34=s.factor(A34.det().as_expr());print('DETERMINANT34',det34,'NORMAL_WEIGHTS',nw34,flush=True)
    out.update(selected34_minor=str(det34),selected34_rows=pr34,selected34_columns=pc34,free34=rest34,normal34=normal34,normal_weights34=nw34,
               T34=[[str(v) for v in row] for row in T34.tolist()],U34=[[str(v) for v in row] for row in U34.tolist()])
    z,uu,bf=primitive_beta(model);bf=bf.subs({s.Symbol(f'u{i+1}'):v for i,v in beta.items()})
    pa,pc,pd,ph,pk=s.symbols('planeA planeC planeD conicH conicK')
    v1=s.Matrix([pa,1,0,pc]);v2=s.Matrix([-pc,0,1,pd]);origin={pa:0,pc:0,pd:0,ph:0,pk:0}
    sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]])
    cj=sec.jacobian([pa,pc,pd]).subs(origin)
    direction=s.Matrix(cj[0,:]).cross(s.Matrix(cj[1,:])).reshape(3,1)
    x,y,zz,w=X[1:];ell1=x-pa*y+pc*zz;ell2=w-pc*y-pd*zz
    fam=(ell1*ell2+ph*ell1**2+pk*ell2**2)*(v1*v2.T-v2*v1.T)
    ds=[];q0=s.zeros(45,1);q0[fixed]=1
    for v in [pa,pc,pd,ph,pk]:
        P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin)
        q=s.Matrix([model.normalize(P)['parameters'][i] for i in two]);q-=q0*q[fixed]
        ds.append(q[cols,:].col_join(s.zeros(1,1)))
    family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4])
    fc=(W.inv()*family)[rest34,:]
    params34=[i for i in range(11) if i not in normal34]
    fdet=s.factor(fc[params34,:].det())
    assert fdet and (J.to_DM()*family.to_DM()).is_zero_matrix
    assert (T34.to_DM()*fc.to_DM()-family.to_DM()).is_zero_matrix
    print('CURVE_JACOBIAN',cj,'FAMILY_MINOR',fdet,'SECONDS',time.time()-start,flush=True)
    out.update(curve_jacobian=[[str(v) for v in row] for row in cj.tolist()],family_parameter_minor=str(fdet),family_parameters=params34)

    assert max(rw[i][0] for i in pr34)<=4
    assert nw34==[1,1,1,1,1,2,1,1]
    Fcore=500*a+90*e**2+3*e-24
    Bcore=50*a*d+18*d*e**2-9*d*e-2*d+60*e-40
    Ccore=100*a+18*e**2-9*e-2
    assert s.cancel(det34-a**3*Fcore*Bcore/s.Integer(1220703125))==0
    assert s.cancel(fdet+Ccore/25)==0
    # All denominators are controlled by the selected Jacobian determinant.
    allowed=[a,Fcore,Bcore]
    for value in list(T34)+list(U34):
        den=s.fraction(s.cancel(value))[1]
        for factor in allowed:
            while s.rem(den,factor,a,e,d)==0:
                den=s.div(den,factor,a,e,d)[0]
        assert not den.free_symbols,den
    paths=[Path(__file__),NOTE,root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',root/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json',root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/primitive_jet_elliptic_curve_20261101.py']
    out.update(id='S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102',all_passed=True,
        core_variables=['a','e','d'],nuisance_parameters_zero=True,
        beta={f'u{i+1}':str(v) for i,v in beta.items()},local_coordinate_weights=cw,
        equation_weights=[v[0] for v in rw],maximum_selected_equation_weight=4,
        full_linear_kernel_checked=True,selected_right_inverse_checked=True,
        denominators_supported_on_a_F_B=True,
        normal_weight_counts={'1':7,'2':1},weighted_column_counts={'2':29,'3':91,'4':239},
        implicit_solution_exact_degree_bound=4,implicit_polynomials_expanded=False,
        full_remaining_equation_weight_bound=max(v[0] for v in rw),
        proof_note=str(NOTE.relative_to(root)),
        limitations=['Only nuisance-zero on the explicit a*F*B*C open.',
                     'No low-weight or quartic rank coverage or colength65 is asserted.',
                     'No full K=0 source coverage or new ambient saturation removal.'],
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    result=verify()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','implicit_solution_exact_degree_bound','seconds']}))
