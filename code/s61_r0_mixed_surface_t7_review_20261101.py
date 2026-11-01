#!/usr/bin/env python3
"""Independent raw135 three-jet and quartic audit of the T7 surface theorem."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();r=json.loads((ROOT/'reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json').read_text());o=json.loads((ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text())['orbits']['T7'];u=r['universal'];p=r['point'];ld=p['local'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq] for eq in eqs]
    x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=3*w*w;Q[1,3]=-x*x;Q-=Q.T
    def coords(q):
        P=X[0]**2*J;P[1:,1:]+=q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    q=coords(Q);assert q==s.Matrix([s.Rational(o['point61'][i]) for i in two]);fixed=next(i for i,c in enumerate(q) if c);cols=[i for i in range(45) if i!=fixed];H=M(u['beta_kernel']);b=s.symbols('b0:11');beta=H*s.Matrix(b);field=s.QQ.frac_field(*b);zero=field.zero;oneK=field.one
    def dm(mat):return mat.to_DM().convert_to(field)
    coeff=[];products={}
    for eq in encoded:
        ff={}
        for mm,bm,c in eq:
            if bm not in products:
                value=oneK
                for k in bm:value*=field.from_sympy(beta[k])
                products[bm]=value
            value=field.from_sympy(c)*products[bm]
            if value:ff[mm]=ff.get(mm,zero)+value
        coeff.append(ff)
    JD=[[zero]*45 for _ in range(135)]
    qk=[field.from_sympy(c) for c in q]
    for row,ff in enumerate(coeff):
        point=zero
        for mm,c in ff.items():
            if len(mm)==2:
                i,j=mm;point+=c*qk[i]*qk[j]
                if i!=fixed:JD[row][cols.index(i)]+=c*qk[j]
                if j!=fixed:JD[row][cols.index(j)]+=c*qk[i]
            elif len(mm)==1:JD[row][-1]+=c*qk[mm[0]]
        assert point==zero
    JD=DM(JD,(135,45),field);assert (JD-dm(M(u['jacobian']))).is_zero_matrix;T=dm(M(u['tangent']));W2=dm(M(u['quadratic_correction']));W3=dm(M(u['cubic_correction']));assert (JD*T).is_zero_matrix and T.extract(u['free_columns'],range(3))==DM.eye(3,field);assert T.extract([44],range(3)).is_zero_matrix and W2.extract([44],range(6)).is_zero_matrix and W3.extract([44],range(10)).is_zero_matrix
    # Use graded46-coordinate arrays and original polynomial products.
    E={};mon={d:list(combinations_with_replacement(range(3),d)) for d in [1,2,3]}
    for d,dd in [(1,T),(2,W2),(3,W3)]:
        arr=[[zero]*len(mon[d]) for _ in range(46)];lst=dd.to_list()
        for j,i in enumerate(cols):arr[i]=lst[j]
        arr[45]=lst[-1];E[d]=arr
    for degree,W in [(2,W2),(3,W3)]:
        ix={mm:j for j,mm in enumerate(mon[degree])};raw=[[zero]*len(ix) for _ in range(135)]
        for row,ff in enumerate(coeff):
            for mm,c in ff.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if fixed in (i,j):continue
                for a in range(1,degree):
                    bb=degree-a
                    for ia,ca in enumerate(E[a][i]):
                        if not ca:continue
                        for ib,cb in enumerate(E[bb][j]):
                            if cb:raw[row][ix[tuple(sorted(mon[a][ia]+mon[bb][ib]))]]+=c*ca*cb
        raw=DM(raw,(135,len(ix)),field);assert (JD*W+raw).is_zero_matrix;print('universal raw135 degree',degree,'PASS',time.time()-start,flush=True)
    sample=dict(zip(b,range(2,13)));denoms=set()
    for mat in [M(u['tangent']),M(u['quadratic_correction']),M(u['cubic_correction'])]:
        for expr in mat:
            de=s.denom(expr);assert de.xreplace(sample)!=0;denoms.add(str(de))
    J0=M(u['jacobian']).xreplace(sample);T0=M(u['tangent']).subs(sample);assert J0==M(ld['jacobian']) and J0.rank()==42 and T0==M(ld['tangent_basis'])
    # Rebuild both source differentials without reading their saved coefficients.
    zz,bu,be=primitive_beta(m);vec=s.Matrix([be[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,A)[rr,:] for A in sp_basis()]);BD=H.row_join(acts.xreplace(dict(zip(bu,beta.xreplace(sample)))));assert BD==M(o['source_beta_differential']) and BD.rank()==16
    OD=[]
    for A in sp_basis():
        dq=coords(sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T);dq-=q*dq[fixed]/q[fixed];OD.append(dq[cols,:])
    OD=s.Matrix.hstack(*OD);assert OD==M(o['q_orbit_differential']);curve=s.zeros(44,11).row_join(OD)*s.Matrix.hstack(*BD.nullspace());assert curve.rank()==2;cp=list(curve.to_DM().rref()[1]);curve=curve[:,cp].col_join(s.zeros(1,2));tr=list(T0.T.to_DM().rref()[1]);ct=T0[tr,:].inv()*curve[tr,:];assert T0*ct==curve and ct==s.Matrix([[0,1],[-2,0],[0,0]])
    beta_minor=BD[:,p['beta_minor_columns']].det();assert beta_minor==s.Rational(-6367977,800)
    # Independently solve the scalar normal arc using raw numeric polynomial coefficients.
    numeric=[];values={str(k):s.Rational(v) for k,v in sample.items()}
    for ff in coeff:
        out={}
        for mm,c in ff.items():out[mm]=s.cancel(field.to_sympy(c)).xreplace(sample)
        numeric.append(out)
    rows=list(J0.T.to_DM().rref()[1]);piv=list(J0.to_DM().rref()[1]);U=s.zeros(45,42);iv=J0[rows,piv].inv(method='DM')
    for i,col in enumerate(piv):U[col,:]=iv[i,:]
    W={1:T0[:,2]};res={}
    for degree in [2,3,4]:
        E={}
        for dd,v in W.items():
            e=s.zeros(46,1)
            for j,i in enumerate(cols):e[i]=v[j]
            e[45]=v[-1];E[dd]=e
        raw=s.zeros(135,1)
        for row,ff in enumerate(numeric):
            for mm,c in ff.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if fixed in (i,j):continue
                for a in range(1,degree):raw[row]+=c*E[a][i]*E[degree-a][j]
        W[degree]=-U*raw[rows,:];res[degree]=raw+J0*W[degree];assert W[degree]==M(p['transverse_corrections'][str(degree)]) and res[degree]==M(p['transverse_residuals'][str(degree)])
    assert res[2]==res[3]==s.zeros(135,1);quartic=res[4][103];assert quartic==s.Rational(-48390411675666648790195,251191872974592);assert all(W[d][-1]==0 for d in [1,2,3,4])
    paths=[Path(__file__),ROOT/'code/s61_r0_mixed_surface_t7_formal_20261101.py',ROOT/'reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json',ROOT/'research/notes/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.md',ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json',ROOT/'reports/S61_R0_MIXED_ORBIT_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_R0_MIXED_SURFACE_T7_REVIEW_20261101',all_passed=True,raw135_universal_jacobian_verified=True,universal_rank=42,tangent_identity_rows=u['free_columns'],raw135_universal_quadratic_and_cubic_corrections_verified=True,universal_t_zero_through_degree3=True,all_denominators_nonzero_at_sample=True,actual_fixed_beta_surface_tangent=[[str(v) for v in row] for row in ct.tolist()],source_beta_minor=str(beta_minor),independent_scalar_quartic_coefficient=str(quartic),completed_ring_accepted='C[[u0,u1,e]]/(e^4), t=0',relative_free_rank4_accepted=True,nonempty_ambient_open_accepted=True,whole_generic_beta_surface_removed=False,scope='Specified point plus nonempty ambient open only; the17D complement can dominate W16.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
