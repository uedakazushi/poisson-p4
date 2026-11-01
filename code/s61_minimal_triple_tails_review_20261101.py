#!/usr/bin/env python3
"""Independent full135 and finite-zero-orbit audit of the M1 triple-tail pair."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_minimal_junction_review_20261101 import standards
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();r=json.loads((ROOT/'reports/S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101.json').read_text());old=json.loads((ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json').read_text())['universal'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq] for eq in eqs];b=s.symbols('b0:10');field=s.QQ.frac_field(*b);zero=field.zero;H=M(old['beta_kernel']);beta=H*s.Matrix(b);x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);aa,bb,dd=s.symbols('a b d')
    def tensor(a,bb,d):
        Q=s.zeros(4);Q[0,1]=Q[2,3]=x*x/2;Q[1,2]=-x*z-bb*x*x;Q[1,3]=a*x*x+2*d*x*z;return Q-Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    def polynomials(beta):
        products={};out=[]
        for eq in encoded:
            f={}
            for mm,bm,c in eq:
                if bm not in products:products[bm]=s.prod(beta[i] for i in bm)
                if products[bm]:f.setdefault(mm,[]).append(c*products[bm])
            out.append({mm:s.expand(s.Add(*v)) for mm,v in f.items()})
        return out
    def raw_jac(F,q,cols):
        A=s.zeros(135,45)
        for row,f in enumerate(F):
            for mm,c in f.items():
                if len(mm)==2:
                    i,j=mm
                    if i in cols:A[row,cols.index(i)]+=c*q[j]
                    if j in cols:A[row,cols.index(j)]+=c*q[i]
                elif len(mm)==1:A[row,-1]+=c*q[mm[0]]
        return A.applyfunc(s.expand)
    q=coords(tensor(1,0,0));cols=old['chart_columns'];fixed=next(i for i in range(45) if i not in cols);F=polynomials(beta);A=raw_jac(F,q,cols);assert A[:,:44]==M(old['leading_jacobian']);T=M(old['tangent']);L=M(old['quadratic_lift']);assert T[[26,36],:]==s.eye(2);AD=A.to_DM().convert_to(field);TD=T.col_join(s.zeros(1,2)).to_DM().convert_to(field);assert (AD*TD).is_zero_matrix
    C=M(old['first_kernel']);pc=list(A[:70,:44].rref()[1]);pr=list(A[:70,pc].T.rref()[1]);U=s.zeros(44,len(pc));iv=A[pr,pc].inv()
    for j,i in enumerate(pc):U[i,:]=iv[j,:]
    R=A[:,:44]*C;ct=(A[:,-1]-A[:,:44]*U*A[pr,-1]).applyfunc(s.cancel);E=R.row_join(ct);src=r['universal'];assert (E.to_DM()-M(src['full_residual']).to_DM()).is_zero_matrix;N=M(src['full_kernel']);assert (E.to_DM()*N.to_DM()).is_zero_matrix and N.rank()==2 and N[-1,:]==s.zeros(1,2);rp=src['full_minor_rows'];cp=src['full_minor_columns'];ed=E[rp,cp].to_DM();assert s.factor(ed.domain.to_sympy(ed.det()))==s.sympify(src['full_minor'])
    allowed={str(s.Poly(f,*b).monic().as_expr()) for f,e in s.factor_list(s.sympify(src['full_minor']),*b)[1]}
    for mat in [T,L,N]:
        for value in set(mat):assert {str(s.Poly(f,*b).monic().as_expr()) for f,e in s.factor_list(s.denom(s.cancel(value)),*b)[1]}<=allowed
    # Raw symmetric quadratic restriction and the full old lift, with t=0.
    Q2=s.zeros(135,3);rawc=[]
    for row,f in enumerate(F):
        for mm,c in f.items():
            if len(mm)!=2 or fixed in mm:continue
            i,j=[cols.index(k) for k in mm]
            for h,(ia,ib) in enumerate([(0,0),(0,1),(1,1)]):Q2[row,h]+=c*(T[i,ia]*T[j,ib]+(T[i,ib]*T[j,ia] if ia!=ib else 0))
    assert (A[:,:44].to_DM()*L.to_DM()-Q2.to_DM()).is_zero_matrix
    assert L[:,0]==s.zeros(44,1)
    # Clear the pure normal first/second corrections before polynomial products.
    from sympy.polys.rings import ring
    PR,*gg=ring(','.join(map(str,b)),s.QQ)
    vv=T[:,1];ww=-L[:,2];dv=s.factor(s.lcm([s.denom(s.cancel(e)) for e in vv]));dw=s.factor(s.lcm([s.denom(s.cancel(e)) for e in ww]));vp=[PR.from_expr(s.cancel(e*dv)) for e in vv];wp=[PR.from_expr(s.cancel(e*dw)) for e in ww];raw=[]
    for f in F:
        value=PR.zero
        for mm,c in f.items():
            if len(mm)==2 and fixed not in mm:
                i,j=[cols.index(k) for k in mm];value+=PR.from_expr(c)*(vp[i]*wp[j]+vp[j]*wp[i])
        raw.append(value)
    AU=(A[:,:44]*U).applyfunc(s.expand);saved=src['cubic'][0];assert s.factor(dv*dw-s.sympify(saved['den']))==0;ZZ=M(saved['z_numerator'])
    for row in range(135):
        value=raw[row]-sum((PR.from_expr(AU[row,j])*raw[k] for j,k in enumerate(pr)),PR.zero);assert value==PR.from_expr(ZZ[row])
    aug=E[:,saved['cols']].row_join(ZZ);md=aug[saved['rows'],:].to_DM();assert s.factor(md.domain.to_sympy(md.det()))==s.sympify(saved['minor'])
    for expression in [s.sympify(saved['den']),s.sympify(saved['minor'])]:assert {str(s.Poly(f,*b).monic().as_expr()) for f,e in s.factor_list(expression,*b)[1]}<=allowed
    print('universal raw135 minors and cubic membership PASS',time.time()-start,flush=True)
    zz,bu,be=primitive_beta(m);vec=s.Matrix([be[i,j].diff(zz[k]) for (i,j),k in COORDS]);base=vec.jacobian(bu);rr=list(base.T.rref()[1]);left=base[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,B)[rr,:] for B in sp_basis()]);sample=dict(zip(b,range(1,11)));geometry=[]
    def geom(H,a0):
        Q=tensor(a0,0,0);q=coords(Q);fi=next(i for i,c in enumerate(q) if c);orb=[]
        for B in sp_basis():
            dq=coords(sum(((B*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-B*Q-Q*B.T);dq-=q*dq[fi]/q[fi];orb.append(dq)
        OD=s.Matrix.hstack(*orb);st=s.Matrix.hstack(*OD.nullspace());bv=H*s.Matrix(b);ar=acts.xreplace(dict(zip(bu,bv)));hr=list(H.T.rref()[1]);hi=H[hr,:].inv();ba=(hi*(ar*st)[hr,:]).applyfunc(s.expand);assert (H*ba-ar*st).applyfunc(s.expand)==s.zeros(16,st.cols);D=H.row_join(ar);joint=D.col_join(s.zeros(45,10).row_join(OD));assert D.xreplace(sample).rank()==16 and joint.xreplace(sample).rank()==16+a0;assert OD.rank()==6+a0;return bv,ar,ba,hi,hr
    beta,ar,ba,hi,hr=geom(H,1);factorout=[]
    for f in r['factors']:
        fac=s.sympify(f['factor']);var=s.sympify(f['solved_variable']);val=s.sympify(f['solved_value']);assert s.degree(fac,var)==1 and s.gcd(s.diff(fac,var),fac.subs(var,0))==1 and s.cancel(fac.subs(var,val))==0
        for j in range(3):assert not s.cancel(sum(s.diff(fac,bb)*ba[i,j] for i,bb in enumerate(b))/fac).free_symbols.intersection(b)
        bp=beta.subs(var,val);pars=[bb for bb in b if bb!=var];D=bp.jacobian(pars).row_join(ar.subs(var,val));rank=D.subs(sample).rank();out=dict(factor=str(fac),beta_rank=rank,type=f['type'])
        if f['type']=='proper_beta_image':
            NN=M(f['certificate']['left_kernel']);assert (D.T.to_DM()*NN.to_DM()).is_zero_matrix and NN.subs(sample).rank()==16-rank and rank==f['certificate']['image_dimension_upper']
        elif f['type']=='cubic_artifact_stratum':
            ld=f['local_data'];beta0=s.Matrix([s.Rational(ld['point61'][i]) for i in one]);QQ=s.Matrix([s.Rational(ld['point61'][i]) for i in two]);ff=polynomials(beta0);J0=raw_jac(ff,QQ,ld['chart_columns']);assert J0==M(ld['jacobian']) and J0.rank()==43 and J0[:,:44].rank()==42;out['independent_jet']=point_jet(ff,J0,M(ld['tangent_basis']),ld['chart_columns'],f['initial']);assert out['independent_jet']['localized_length']==3
        else:assert fac==b[0]
        factorout.append(out)
    z0=r['zero_orbit'];H0=M(z0['beta_kernel']);bv0,ar0,ba0,hi0,hr0=geom(H0,0);assert bv0[:9,:]==beta[:9,:] and bv0[10:,:]==beta[10:,:] and bv0[9]==0
    # Both actual fixed-beta lines and the complete closure family.
    for a0,HH,cc in [(1,H,r['same_beta_curve']),(0,H0,z0['same_beta_curve'])]:
        betaA=HH*s.Matrix(b);F0=polynomials(betaA);qc=coords(tensor(a0+2*b[9]*dd/(3*b[7]),0,dd));assert all(s.cancel(sum(co*s.prod(qc[i] for i in mm) for mm,co in f.items() if len(mm)==2))==0 for f in F0)
    bf=list(bv0);bf[9]=3*aa*b[7];qc=coords(tensor(aa,0,0));assert all(s.expand(sum(co*s.prod(qc[i] for i in mm) for mm,co in f.items() if len(mm)==2))==0 for f in polynomials(bf))
    ld=z0['local_data'];beta0=s.Matrix([s.Rational(ld['point61'][i]) for i in one]);QQ=s.Matrix([s.Rational(ld['point61'][i]) for i in two]);ff=polynomials(beta0);J0=raw_jac(ff,QQ,ld['chart_columns']);assert J0==M(ld['jacobian']) and J0.rank()==43 and J0[:,:44].rank()==42;zjet=point_jet(ff,J0,M(ld['tangent_basis']),ld['chart_columns'],z0['initial']);assert zjet['zero_section_length']==4 and zjet['localized_length']==3
    qc=coords(tensor(2*b[9]*dd/(3*b[7]),0,dd)).subs(sample);dq=qc.diff(dd);fi=next(i for i in range(45) if i not in ld['chart_columns']);dq-=QQ*dq[fi]/QQ[fi];assert M(ld['tangent_basis'])*s.Matrix([2,0])==dq[ld['chart_columns'],:].col_join(s.zeros(1,1))
    for norm,target in [(M(r['normalization']['finite_shear']),tensor(aa-2*bb*dd,0,0)),(M(r['normalization']['nonzero_orbit']['matrix']),tensor(1,0,0))]:
        assert (norm.T*J[1:,1:]*norm-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4);pull=norm.inv()*tensor(aa,bb,dd).subs(dict(zip(V,norm*V)),simultaneous=True)*norm.inv().T;assert (pull-target).applyfunc(s.cancel)==s.zeros(4)
    main=json.loads((ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json').read_text())['cases']['M1'];P,RR,AA,BB,UU,VV,KK,LL=s.symbols('P R A B U V K L');ss={P:beta[1],RR:beta[3],AA:beta[4],BB:beta[13],UU:beta[6],VV:beta[9],KK:beta[14]+5*beta[9],LL:5*beta[6]-beta[15]};ma=s.cancel(s.sympify(main['solution_a']).subs(ss,simultaneous=True));mb=s.cancel(s.sympify(main['solution_b']).subs(ss,simultaneous=True));assert ma.subs(b[0],0)==1 and mb==0
    paths=[Path(__file__),ROOT/'code/s61_minimal_triple_tails_removal_20261101.py',ROOT/'reports/S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101.json',ROOT/'research/notes/S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101.md',ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json',ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_REVIEW_20261101.json',ROOT/'code/s61_minimal_two_variable_initial_20261101.py',ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/s61_minimal_junction_review_20261101.py']
    return dict(id='S61_MINIMAL_TRIPLE_TAILS_REVIEW_20261101',all_passed=True,universal_raw135_jacobian_quadratic_lift_and_two_pure_cubics_verified=True,first_pure_cubic_membership_checked_over_function_field=True,all_open_minor_denominator_factors_covered=True,factor_cases=factorout,zero_orbit=dict(orbit_dimension=6,joint_dimension=16,actual_tangent=[2,0],jet=zjet,full135_closure_family_verified=True),finite_normalizations_verified=True,main_junction_verified=True,accepted_entire_finite_triple_tail_pair=True,scope='Finite M1 g2/g4 pair including normalization-zero orbit. Zero orbit only t nilpotence, no flatness/full t=0. Binary infinity endpoints and other tails separate.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def point_jet(F,J0,T,cols,saved):
    n=2;fixed=next(i for i in range(45) if i not in cols);depth=saved['initial_degree'];rows=list(J0.T.to_DM().rref()[1]);piv=list(J0.to_DM().rref()[1]);U=s.zeros(45,43);iv=J0[rows,piv].inv(method='DM')
    for j,i in enumerate(piv):U[i,:]=iv[j,:]
    mons={d:list(combinations_with_replacement(range(n),d)) for d in range(1,depth+1)};E={1:s.zeros(46,2)}
    for j,i in enumerate(cols):E[1][i,:]=T[j,:]
    E[1][45,:]=T[-1,:];res={}
    for degree in range(2,depth+1):
        ix={mm:i for i,mm in enumerate(mons[degree])};raw=s.zeros(135,len(ix))
        for row,ff in enumerate(F):
            for mm,c in ff.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if fixed in (i,j):continue
                for a in range(1,degree):
                    bb=degree-a
                    for ia,ca in enumerate(E[a][i,:]):
                        if not ca:continue
                        for ib,cb in enumerate(E[bb][j,:]):
                            if cb:raw[row,ix[tuple(sorted(mons[a][ia]+mons[bb][ib]))]]+=c*ca*cb
        lift=-U*raw[rows,:];res[degree]=raw+J0*lift;assert res[degree]==M(saved['residuals'][str(degree)]);E[degree]=s.zeros(46,len(ix))
        for j,i in enumerate(cols):E[degree][i,:]=lift[j,:]
        E[degree][45,:]=lift[-1,:]
    assert all(res[d]==s.zeros(135,len(mons[d])) for d in range(2,depth));rr,pv=res[depth].to_DM().rref();rr=rr.to_Matrix();v=s.symbols('v0:2');forms=[sum(rr[i,j]*s.prod(v[k] for k in mm) for j,mm in enumerate(mons[depth])) for i in range(len(pv))];assert forms==list(map(s.sympify,saved['initial_forms']));g0=s.groebner([f.subs(v[0],0) for f in forms],v[1]);g1=s.groebner([f.subs(v[0],1) for f in forms],v[1]);st0,hf0=standards(g0);st1,hf1=standards(g1)
    return dict(raw135_all_residuals_verified_through_degree=depth,initial_forms=list(map(str,forms)),first_nonzero=saved['first_nonzero'],zero_section_length=len(st0),localized_length=len(st1))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
