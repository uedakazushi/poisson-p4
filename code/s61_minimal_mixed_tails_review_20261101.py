#!/usr/bin/env python3
"""Independent all-B and projective-endpoint audit of four mixed minimal tails."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def run():
    started=time.time();r=json.loads((ROOT/'reports/S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101.json').read_text());model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};equations=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq] for eq in equations]
    x,y,z,w=X[1:];spatial=s.Matrix([x,y,z,w]);B=s.Symbol('B');b=s.symbols('b0:7');field=s.QQ.frac_field(B,*b);a,bb,c,d,h,rho=s.symbols('a b c d h rho')
    def tensor(kind,a,bb,c,d,h,axis=6):
        G=a*x*x*z+bb*x*x*w+c*x*z*w+d*x*(z*z if axis==6 else w*w)+h*(z*z*w if axis==6 else z*w*w);Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-int(kind=='M1')*x*z;Q[1,3]=s.diff(G,z);return Q-Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(model.normalize(P)['parameters'])[two,:]
    zz,bu,be=primitive_beta(model);vec=s.Matrix([be[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);kr=list(basis.T.rref()[1]);iv=basis[kr,:].inv();actions=s.Matrix.hstack(*[iv*infinitesimal_action(vec,A)[kr,:] for A in sp_basis()])
    def proper(data):
        beta=M(data['beta']);variables=sorted(set().union(*(v.free_symbols for v in beta)),key=lambda v:int(str(v)[1:]));H=beta.jacobian(variables);ar=actions.xreplace(dict(zip(bu,beta)));D=H.row_join(ar);L=M(data['left_kernel']);assert (D.T.to_DM()*L.to_DM()).is_zero_matrix;sample=dict(zip(variables,range(1,len(variables)+1)));assert L.subs(sample).rank()==L.cols;rank=D.subs(sample).rank();assert rank==16-L.cols==data['image_dimension_upper'];return rank
    cases=[]
    for rec in r['cases']:
        kind=rec['kind'];Q=tensor(kind,1,B,0,0,1);q=coords(Q);src=rec['universal'];assert q==M(src['q_point']);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];assert cols==src['chart_columns'];H=M(rec['kernel']['kernel']);beta=H*s.Matrix(b);A=s.zeros(135,44);tc=s.zeros(135,1);BM=s.zeros(134,16);products={}
        for row,eq in enumerate(encoded):
            total=0
            for mm,bm,co in eq:
                if bm not in products:products[bm]=s.prod(beta[i] for i in bm)
                coefficient=co*products[bm]
                if len(mm)==2:
                    i,j=mm;total+=coefficient*q[i]*q[j]
                    if i!=fixed:A[row,cols.index(i)]+=coefficient*q[j]
                    if j!=fixed:A[row,cols.index(j)]+=coefficient*q[i]
                    if 70<=row<134:
                        assert len(bm)==1;BM[row,bm[0]]+=co*q[i]*q[j]
                elif len(mm)==1:tc[row]+=coefficient*q[mm[0]]
            assert s.expand(total)==0
        A=A.applyfunc(s.expand);tc=tc.applyfunc(s.expand);BM=BM.applyfunc(s.expand);assert (A-M(src['leading_jacobian'])).applyfunc(s.expand)==s.zeros(135,44) and (tc-M(src['t_column'])).applyfunc(s.expand)==s.zeros(135,1) and (BM-M(rec['kernel']['matrix'])).applyfunc(s.expand)==s.zeros(134,16);assert (BM*H).applyfunc(s.expand)==s.zeros(134,7)
        kd=s.factor(H[rec['kernel']['kernel_minor_rows'],:].det());assert kd==s.sympify(rec['kernel']['kernel_minor']) and not kd.has(B) and kd!=0
        dets=[]
        for item in rec['kernel']['rank_minor_cover']:
            determinant=s.factor(BM[item['rows'],item['columns']].det());assert s.expand(determinant-s.sympify(item['determinant']))==0;dets.append(determinant)
        assert s.expand(sum(s.sympify(v)*de for v,de in zip(rec['kernel']['bezout_coefficients'],dets)))==1
        pc=src['first_minor_columns'];pr=src['first_minor_rows'];first=A[pr,pc];d0=s.factor(first.det());assert s.expand(d0-s.sympify(src['first_minor']))==0;C=M(src['first_kernel']);free=[i for i in range(44) if i not in pc];assert C[free,:]==s.eye(12) and (A[:70,:].to_DM()*C.to_DM()).is_zero_matrix
        U=s.zeros(44,32);inverse=first.inv(method='DM')
        for j,i in enumerate(pc):U[i,:]=inverse[j,:]
        EE=(A*C).row_join(tc-A*U*tc[pr,:]).applyfunc(s.cancel);assert (EE.to_DM()-M(src['residual_full_matrix']).to_DM()).is_zero_matrix;assert EE[:70,-1]!=s.zeros(70,1)
        selected=EE[src['full_minor_rows'],src['full_minor_columns']];denom=s.prod(map(s.sympify,src['row_denominators']));cleared=selected.copy()
        for i,den in enumerate(map(s.sympify,src['row_denominators'])):
            for j in range(cleared.cols):cleared[i,j]=s.cancel(cleared[i,j]*den)
        dm=cleared.to_DM();det=s.factor(dm.domain.to_sympy(dm.det())/denom);assert s.cancel(det-s.sympify(src['full_minor']))==0;assert selected.shape==(12,12)
        allowed={s.Poly(f,B,*b).monic().as_expr() for v in [d0,s.numer(det)] for f,e in s.factor_list(v)[1]};expected={b[6],B,B+1} if kind=='M0' else {b[6],B,B+s.Rational(1,2),B-s.Rational(1,4)};assert allowed==expected
        od=[]
        for A0 in sp_basis():
            dq=coords(sum(((A0*spatial)[i]*Q.diff(spatial[i]) for i in range(4)),s.zeros(4))-A0*Q-Q*A0.T);dq-=q*dq[fixed]/q[fixed];od.append(dq)
        OD=s.Matrix.hstack(*od);ar=actions.xreplace(dict(zip(bu,beta)));st=M(rec['geometry']['generic']['stabilizer']);assert (OD*st).applyfunc(s.cancel)==s.zeros(45,1);assert (H*M(rec['geometry']['generic']['stabilizer_action'])-ar*st).applyfunc(s.cancel)==s.zeros(16,1)
        parameters=(*b,B);BD=beta.jacobian(parameters).row_join(ar);joint=BD.col_join(q.jacobian(parameters).row_join(OD));sample={B:2,**dict(zip(b,range(1,8)))};assert OD.subs(B,2).rank()==9 and BD.subs(sample).rank()==16 and joint.subs(sample).rank()==17
        exceptions=[]
        for item in rec['geometry']['fixed_B']:
            val=s.Rational(item['B']);Hb=H.subs(B,val);arb=ar.subs(B,val);Ob=OD.subs(B,val);assert Ob.rank()==9;stab=M(item['stabilizer']);assert Ob*stab==s.zeros(45,1);assert (Hb*M(item['stabilizer_action'])-arb*stab).applyfunc(s.expand)==s.zeros(16,1);D=Hb.row_join(arb);rank=D.subs(sample).rank();assert rank==item['beta_rank']
            if 'proper_image' in item:assert proper(item['proper_image'])==15
            else:
                loc=item['local_data'];pt=list(map(s.Rational,loc['point61']));assert s.Matrix([pt[i] for i in one])==beta.subs(B,val).subs(sample) and s.Matrix([pt[i] for i in two])==q.subs(B,val);JJ=A.row_join(tc).subs({**sample,B:val});assert JJ==M(loc['jacobian']) and JJ.rank()==44
            exceptions.append(dict(B=str(val),beta_image_rank=rank,full_rank_44_sample_verified='local_data' in item))
        assert proper(rec['geometry']['lambda_zero'])==9
        cases.append(dict(kind=kind,all135_matrix_reconstructed=True,all_B_kernel_rank=7,rank9_minor_bezout_cover=len(dets),first_minor=str(d0),residual_minor=str(det),all_bad_factors_exhausted=sorted(map(str,allowed)),nonzero_top70_t_compatibility_retained=True,beta_rank=16,joint_dimension=17,exceptions=exceptions,lambda_zero_beta_image_dimension=9));print(kind,'raw135/allB/minors/geometry PASS',time.time()-started,flush=True)
    coverage=r['coverage']
    for axis in ['6','7']:assert proper(coverage['normalizer_zero_exclusions'][axis])==15
    for kind in ['M0','M1']:
        normal=coverage['normalizations'][kind];T=M(normal['shear']);R=M(normal['torus']);E=M(normal['exchange']);aa=s.sympify(normal['a_reduced']);bt=s.sympify(normal['b_reduced'])
        for C in [T,R,E]:assert (C.T*J[1:,1:]*C-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
        pull=T.inv()*tensor(kind,a,bb,c,d,h).subs(dict(zip(spatial,T*spatial)),simultaneous=True)*T.inv().T;assert (pull-tensor(kind,aa,bt,0,0,h)).applyfunc(s.cancel)==s.zeros(4)
        pull=R.inv()*tensor(kind,aa,bt,0,0,h).subs(dict(zip(spatial,R*spatial)),simultaneous=True)*R.inv().T/rho**2;assert (pull-tensor(kind,aa*rho**2/h,bt*h,0,0,1)).applyfunc(s.cancel)==s.zeros(4)
        paired=-c-int(kind=='M1');pull=E.inv()*tensor(kind,a,bb,c,d,h,axis=7).subs(dict(zip(spatial,E*spatial)),simultaneous=True)*E.inv().T;assert (pull-tensor(kind,bb,-a,paired,d,-h)).applyfunc(s.expand)==s.zeros(4)
    old=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());k=s.Symbol('k');K,L,U,V=s.symbols('K L U V');covered=[]
    for family in old['families']:
        kind=family['label'];normal=coverage['normalizations'][kind];aa=s.sympify(normal['a_reduced']);bt=s.sympify(normal['b_reduced'])
        for tail in family['rational_tails']:
            axis=tail['axis']
            if axis not in [6,7]:continue
            parameter=s.Symbol(tail['parameter']);cf={name:s.sympify(expression) for name,expression in tail['coordinate_formulas'].items()};cf[tail['parameter']]=parameter;assert len(cf)==9
            for name,expr in cf.items():assert s.cancel(expr.subs(parameter,0)-s.sympify(family['baseline_solution'][name]).subs(k,s.Rational(tail['k'])))==0
            sub={a:cf['g0'],bb:cf['g1'],c:cf['g3'],d:cf['g2'],h:parameter} if axis==6 else {a:cf['g1'],bb:-cf['g0'],c:-cf['g3']-int(kind=='M1'),d:cf['g4'],h:-parameter};factor=K-6*V if axis==6 else L-6*U;numerator=s.factor((h*aa).subs(sub,simultaneous=True));ratio=s.cancel(numerator/(parameter*factor));assert ratio!=0 and not ratio.free_symbols
            G=sum(cf[f'g{i}']*mon for i,mon in enumerate([x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]));assert s.degree(G,parameter)==1;Gi=s.diff(G,parameter);rec=next(c0 for c0 in coverage['covered_presentations'] if c0['kind']==kind and c0['axis']==axis);assert s.expand(Gi-s.sympify(rec['infinity_cubic']))==0;endpoint=M(rec['infinity_symplectic_matrix']);assert endpoint.T*J[1:,1:]*endpoint==J[1:,1:]
            Qinf=s.zeros(4);Qinf[1,2]=-s.diff(Gi,w);Qinf[1,3]=s.diff(Gi,z);Qinf-=Qinf.T;pure=z*z*w if axis==6 else z*w*w;target=s.zeros(4);target[1,2]=-s.diff(pure,w);target[1,3]=s.diff(pure,z);target-=target.T;pull=endpoint.inv()*Qinf.subs(dict(zip(spatial,endpoint*spatial)),simultaneous=True)*endpoint.inv().T;assert (pull-target).applyfunc(s.expand)==s.zeros(4)
            covered.append(dict(kind=kind,axis=axis,k=tail['k'],zero_parameter_all9_verified=True,normalizer_zero_factor=str(factor),projective_endpoint_full4x4_verified=True))
    assert len(covered)==4
    paths=[Path(__file__),ROOT/'code/s61_minimal_mixed_tails_removal_20261101.py',ROOT/'reports/S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101.json',ROOT/'research/notes/S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101.md',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101.json',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101.md',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_MIXED_TAILS_REVIEW_20261101',all_passed=True,cases=cases,all_B_and_lambda_covered=True,normalizer_zero_images_at_most15=True,full_rational_normalizations_verified=True,covered_presentations=covered,accepted_all_four_projective_parameter_presentations_removed=True,scope='The four stated mixed g6/g7 tails, including zero and infinity. Other minimal main exceptions and binary/flag strata, and full affine exhaustion remain separate.',seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}))
