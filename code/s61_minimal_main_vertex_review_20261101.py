#!/usr/bin/env python3
"""Independent raw135 rational main-formula and rank-one M1 local audit."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    started=time.time();path=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';src=json.loads(path.read_text());assert src['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];index={i:j for j,i in enumerate(two)};params=list(map(s.sympify,src['source_variables']));bv=M(src['source_beta']);sample={s.sympify(k):s.Rational(v) for k,v in src['sample'].items()};A,B,P,R,U,V,K,L=params[:8];den=s.sympify(src['denominator_factor']);disc=s.sympify(src['vertex_discriminant']);assert s.expand(4*den-25*disc+81*B**2)==0
    eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];poly_ring,*unused=ring(','.join(map(str,params)),QQ);beta_poly=[poly_ring.from_expr(v) for v in bv]
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();actions=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);vertex=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());out=[]
    for kind,record in src['cases'].items():
        qr=M(record['q_coordinates']);common=6*B*den;qn=[poly_ring.from_expr(s.cancel(common*v)) for v in qr];allvals=[beta_poly[one.index(i)] if i in one else qn[index[i]] for i in range(61)]
        for row,eq in enumerate(eqs):
            residual=poly_ring.zero
            for mon,a,b in eq:
                if sum(i in two for i in mon)!=2:continue
                if any(not allvals[i] for i in mon):continue
                term=poly_ring.ground_new(QQ(a,b))
                for i in mon:term*=allvals[i]
                residual+=term
            assert residual==0,(kind,row)
        # Reconstruct the actual tensor and confirm all45 saved coordinates.
        x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);H=s.sympify(record['H']);G=s.sympify(record['solution_a'])*x*x*z+s.sympify(record['solution_b'])*x*x*w+x*H
        Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T;Biv=X[0]**2*J;Biv[1:,1:]+=Q;check=m.normalize(Biv)['parameters'];assert all(s.cancel(check[i]-qr[j])==0 for j,i in enumerate(two))
        assert all(s.denom(s.cancel(qv)).subs(sample)!=0 for qv in qr)
        for qv in qr:
            for fac,power in s.factor_list(s.denom(s.cancel(qv)))[1]:assert s.rem(s.Poly(B*den,*params),s.Poly(fac,*params)).is_zero
        k=s.Symbol('k');diag={P:0,R:0,A:k,B:2*k+1};block=next(v for v in vertex['families'] if v['label']==kind);baseline={s.sympify(i):s.sympify(v) for i,v in block['baseline_solution'].items()};gs=s.symbols('g0:9');assert s.cancel(s.sympify(record['solution_a']).subs(diag)-baseline[gs[0]])==0 and s.cancel(s.sympify(record['solution_b']).subs(diag)-baseline[gs[1]])==0
        rankone={P:-1,R:0,A:1,B:2,K:0,L:0};rankone_saved={s.sympify(i):s.sympify(v).subs(gs[8],0) for i,v in block['rankone_solution'].items()};rankone_saved[gs[8]]=0;mons=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3];assert s.cancel(G.subs(rankone)-sum(rankone_saved[i]*mm for i,mm in zip(gs,mons)))==0
        ld=record['local_data'];point=list(map(s.Rational,ld['point61']));q0=qr.xreplace(sample);beta=bv.xreplace(sample);assert q0==s.Matrix([point[i] for i in two]) and beta==s.Matrix([point[i] for i in one]);cols=ld['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
        F=[]
        for eq in eqs:
            f={}
            for mon,a,b in eq:
                qm=tuple(index[i] for i in mon if i in index);cc=s.Rational(a,b)*s.prod(point[i] for i in mon if i in one);f[qm]=f.get(qm,s.S.Zero)+cc
            F.append(f)
        Jac=s.zeros(135,45)
        for row,f in enumerate(F):
            assert sum(cc*s.prod(q0[i] for i in mm) for mm,cc in f.items() if len(mm)==2)==0
            for mm,cc in f.items():
                if len(mm)==2:
                    for j,i in enumerate(cols):Jac[row,j]+=cc*sum(q0[mm[1-h]] for h in range(2) if mm[h]==i)
                elif len(mm)==1:Jac[row,-1]+=cc*q0[mm[0]]
        assert Jac==M(ld['jacobian']);T=M(ld['tangent_basis']);assert Jac*T==s.zeros(135,T.cols)
        BD=bv.jacobian(params).row_join(actions.subs(dict(zip(bu,beta))));assert BD==M(record['beta_differential'])
        dQ=qr.jacobian(params).xreplace(sample);dQ-=q0*dQ[fixed,:]/q0[fixed];Q0=Q.xreplace(sample);orbit=[]
        for HH in sp_basis():
            delta=sum(((HH*Y)[i]*Q0.diff(Y[i]) for i in range(4)),s.zeros(4))-HH*Q0-Q0*HH.T;PP=X[0]**2*J;PP[1:,1:]+=delta;vv=m.normalize(PP)['parameters'];dq=s.Matrix([vv[i] for i in two]);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[cols,:])
        QD=dQ[cols,:].row_join(s.Matrix.hstack(*orbit));assert QD==M(record['q_differential']);assert BD.rank()==16 and BD.col_join(QD).rank()==17
        dd=s.Matrix([[s.diff(disc,v).xreplace(sample) for v in params]+[0]*10]);assert [i for i,v in enumerate(dd) if v]==[3]
        restricted=[i for i in range(24) if i!=3];assert BD[:,restricted].rank()==16 and BD[:,restricted].col_join(QD[:,restricted]).rank()==16
        direction=M(record['source_direction']);assert BD*direction==s.zeros(16,1) and (dd*direction)[0]==s.Rational(record['discriminant_derivative'])!=0
        tangent=M(record['main_branch_tangent']);assert T*tangent==(QD*direction).col_join(s.zeros(1,1))
        data=dict(kind=kind,all135_universal_identities=True,denominators_only_B_and_E=True,diagonal_and_rankone_matches=True,source_beta_rank=16,joint_rank=17,rankone_joint_rank=16,discriminant_derivative=str((dd*direction)[0]),leading_rank=Jac[:,:44].to_DM().rank(),full_rank=Jac.to_DM().rank())
        if kind=='M1':
            assert data['leading_rank']==42 and data['full_rank']==43 and T.cols==2 and T[-1,:]==s.zeros(1,2)
            ET=s.zeros(46,2)
            for j,i in enumerate(cols):ET[i,:]=T[j,:]
            ET[45,:]=T[-1,:];HQ=s.zeros(135,3)
            for row,f in enumerate(F):
                for mm,cc in f.items():
                    if fixed in mm:continue
                    i,j=(*mm,*([45]*(2-len(mm))));HQ[row,0]+=cc*ET[i,0]*ET[j,0];HQ[row,1]+=cc*(ET[i,0]*ET[j,1]+ET[i,1]*ET[j,0]);HQ[row,2]+=cc*ET[i,1]*ET[j,1]
            LF=s.Matrix.hstack(*Jac.T.nullspace()).T;OP=LF*HQ;assert OP.to_DM().rank()==1;form=next(row for row in OP.tolist() if any(row));hess=s.Matrix([[2*form[0],form[1]],[form[1],2*form[2]]]);assert hess.rank()==1 and hess*tangent==s.zeros(2,1);data['M1_square_obstruction_verified']=True
        out.append(data)
    paths=[Path(__file__),path,ROOT/'code/s61_minimal_main_vertex_formula_20261101.py',ROOT/'research/notes/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.md',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_MINIMAL_MAIN_VERTEX_REVIEW_20261101',all_passed=True,cases=out,accepted_M1_rankone_point_ring='C[[u,e]]/(e^2), t=0',accepted_M1_generic_rankone_stratum_removed=True,scope='Both universal main formulas; only M1 rank-one saturation, not M0 A1 junction.',seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
