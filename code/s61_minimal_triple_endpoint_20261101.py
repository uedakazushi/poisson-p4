#!/usr/bin/env python3
"""Remove the triple minimal-tail binary endpoint using good reduction."""
import argparse,hashlib,json,os,shlex,shutil,subprocess,tempfile,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import numpy as np
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_minimal_baseline_thickenings_20261101 import coords,enc
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_infinity_local_data_20261101 import calculate
from weight1_structure import sp_basis
from truncated_local_ideal_io_20261101 import write_input,read_initial
ROOT=Path(__file__).resolve().parents[1]
PRIME=32003
def point_jets():
    sample=dict(zip(('u4','u5','u12','u13','u14','u25','u26','u27','u34','u35','u36','u37','u48','u49','u50','u51'),map(str,(0,0,0,0,6,2,2,3,5,3,7,11,13,3,-1,5))));expected_parameters=[0,3,8]
    start=time.time();p=32003;m=S61()
    r={'sample':dict(sample)}
    b={int(k[1:])-1:s.Rational(v) for k,v in r['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];q0=s.zeros(45,1);q0[fixed]=1
    J=s.zeros(135,45)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon
       if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
       if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
      elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
    T=s.Matrix.hstack(*J.nullspace());n=T.cols;assert n==10 and J*T==s.zeros(135,10) and T[-1,:]==s.zeros(1,10)
    _,piv=J.to_DM().rref();free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(n)
    z,uu,bf=primitive_beta(m);bf=bf.subs({s.Symbol(k):s.Rational(v) for k,v in r['sample'].items()});a,c,d,h,k=s.symbols('a c d h k');v1=s.Matrix([a,1,0,c]);v2=s.Matrix([-c,0,1,d]);origin={a:0,c:0,d:0,h:0,k:0}
    sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]])
    assert sec.subs(origin)==s.zeros(2,1)
    curve_jac=sec.jacobian([a,c,d]).subs(origin);assert curve_jac.rank()==2
    direction=curve_jac.nullspace()[0]
    x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell1*ell2+h*ell1**2+k*ell2**2)*(v1*v2.T-v2*v1.T)
    ds=[]
    for v in [a,c,d,h,k]:
     P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);qs=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);qs-=q0*qs[fixed];ds.append(qs[cols,:].col_join(s.zeros(1,1)))
    family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);coordinates=family[free,:];assert T*coordinates==family and J*family==s.zeros(135,3)
    params=expected_parameters;assert coordinates[params,:].det()!=0;normal=[i for i in range(n) if i not in params];assert params==expected_parameters and len(normal)==7
    def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p
    def arr(M):return np.array([[mod(v) for v in row] for row in M.tolist()],dtype=np.int64)
    def rref(A):
     A=A.copy()%p;row=0;piv=[]
     for col in range(A.shape[1]):
      choices=np.flatnonzero(A[row:,col])
      if len(choices)==0:continue
      j=row+int(choices[0]);A[[row,j]]=A[[j,row]];A[row]=A[row]*pow(int(A[row,col]),-1,p)%p
      mask=np.flatnonzero(A[:,col]);mask=mask[mask!=row]
      if len(mask):A[mask]=(A[mask]-A[mask,col,None]*A[row,None,:])%p
      piv.append(col);row+=1
      if row==A.shape[0]:break
     return A,piv
    j=arr(J);t=arr(T);_,C=rref(j);_,R=rref(j[:,C].T);rank=len(C);rr,_=rref(np.concatenate([j[np.ix_(R,C)],np.eye(rank,dtype=np.int64)],axis=1));U=np.zeros((45,rank),dtype=np.int64);U[C,:]=rr[:,rank:];JU=j@U%p
    assert len(C)==len(R)==35 and np.array_equal(j[R,:]@U%p,np.eye(35,dtype=np.int64))
    H=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      i,l=[cols.index(v) for v in mon]+[44]*(2-len(mon));v=mod(c)
      if i==l:H[row,i,l]+=v
      else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p

    normal=list(range(T.cols))
    H%=p;n=len(normal);W={1:t[:,normal]};A={};ms={d:list(combinations_with_replacement(range(n),d)) for d in range(1,8)};idx={d:{m:i for i,m in enumerate(ms[d])} for d in range(1,8)}
    for degree in range(2,8):
     Qd=np.zeros((135,len(ms[degree])),dtype=np.int64)
     for da in range(1,degree//2+1):
      db=degree-da;fac=1 if da==db else 2
      indices=np.array([idx[degree][tuple(sorted(ma+mb))] for ma in ms[da] for mb in ms[db]])
      for row in range(135):
       cross=(W[da].T@(H[row]@W[db]%p))%p
       acc=np.bincount(indices,weights=cross.ravel(),minlength=len(ms[degree])).astype(np.int64)
       Qd[row]=(Qd[row]+fac*acc)%p
     W[degree]=-U@Qd[R,:]%p;A[degree]=(Qd-JU@Qd[R,:])%p
     print('degree',degree,'computed',time.time()-start,flush=True)
    return dict(sample=sample,J=J,T=T,F=F,two=two,fixed=fixed,chart=cols,free=free,rows=R,columns=C,parameters=params,family=coordinates,curve_jacobian=curve_jac,A=A,tseries={d:W[d][-1,:].tolist() for d in A},monomials=ms)


def universal_thickening(m):
    """Pull the previously proved quadratic nilpotent deformation to infinity."""
    from weight1_structure import COORDS
    from s61_minimal_cubic_tail_20261101 import tensor
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ
    start=time.time();old=json.loads((ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json').read_text())['universal']
    K,L,V,d,z=s.symbols('K L V d s');ns=s.symbols('n0:6');source=s.Matrix([0,0,0,0,-s.Rational(1,3),ns[0],0,ns[1],ns[2],V,ns[3],ns[4],ns[5],s.Rational(1,3),K-5*V,-L])
    zz,u,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(u);rows=list(basis.T.to_DM().rref()[1]);iv=basis[rows,:].inv()
    a=V-3*L*d;b=-L;g=a-2*b*d;T1=s.eye(4);T1[3,2]=-2*d;T2=s.eye(4);T2[1,3]=b;T2[2,0]=-b;Tn=T1*T2*s.diag(1,1,1/g,g)
    To=s.Matrix([[1,0,0,0],[0,1,-3*L/2,0],[3*L/2,0,0,-1],[0,0,1,0]])
    assert (Tn*J0[1:,1:]*Tn.T-J0[1:,1:]).applyfunc(s.cancel)==s.zeros(4) and To*J0[1:,1:]*To.T==J0[1:,1:]
    BB=B.subs(dict(zip(u,source)))
    def pull_beta(T):
        out=T.T*BB.subs(dict(zip(zz,T*s.Matrix(zz))),simultaneous=True)*T
        return (iv*s.Matrix([out[i,j].diff(zz[k])for(i,j),k in COORDS])[rows,:]).applyfunc(s.factor)
    canonical=pull_beta(Tn);beta=pull_beta(To);M=s.Matrix(old['beta_kernel']);rr=list(M.T.to_DM().rref()[1]);bb=(M[rr,:].inv()*canonical[rr,:]).applyfunc(s.factor)
    assert(M*bb-canonical).applyfunc(s.cancel)==s.zeros(16,1)
    sub=dict(zip(s.symbols('b0:10'),bb));nc=s.Matrix(old['tangent'])[:,1].applyfunc(s.sympify);wc=-s.Matrix(old['quadratic_lift'])[:,2].applyfunc(s.sympify)
    nc=nc.applyfunc(lambda v:s.factor(s.cancel(v.subs(sub,simultaneous=True))));wc=wc.applyfunc(lambda v:s.factor(s.cancel(v.subs(sub,simultaneous=True))));print('canonical nilpotent coefficients',time.time()-start,flush=True)
    S=(Tn.inv()*To).applyfunc(s.cancel);Si=S.inv().applyfunc(s.cancel);xv=s.Matrix(X[1:]);two=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==2]
    outn=s.zeros(45,1);outw=s.zeros(45,1)
    for j,i in enumerate(old['chart_columns']):
        if nc[j]==0 and wc[j]==0:continue
        par=[s.S.Zero]*61;par[two[i]]=1;P=m.reconstruct(par);Q=s.zeros(4)
        for h in range(4):
            for k in range(4):Q[h,k]=s.expand((P[h+1,k+1]-X[h+1]*P[0,k+1]+X[k+1]*P[0,h+1]).subs(X[0],1)-J0[h+1,k+1])
        assert s.Matrix(coords(m,Q))[two,:]==s.eye(45)[:,i]
        transformed=(Si*Q.subs(dict(zip(xv,S*xv)),simultaneous=True)*Si.T).applyfunc(s.cancel);qp=s.Matrix(coords(m,transformed))[two,:]
        if nc[j]!=0:outn+=qp*nc[j]
        if wc[j]!=0:outw+=qp*wc[j]
    outn=outn.applyfunc(s.cancel);outw=outw.applyfunc(s.cancel);print('nilpotent pullback',time.time()-start,flush=True)
    Qt=(To.inv()*tensor(a,b,d).subs(dict(zip(xv,To*xv)),simultaneous=True)*To.inv().T).applyfunc(s.expand);q=s.Matrix(coords(m,Qt))[two,:];fi=two.index(29);q0=-q.diff(d)/2;qc=q.subs(d,0)
    assert q0[fi]==1 and qc[fi]==0;delta=qc;qbase=q0+z*delta
    nn=outn.applyfunc(lambda v:s.cancel(v.subs(d,-1/(2*z))));ww=outw.applyfunc(lambda v:s.cancel(v.subs(d,-1/(2*z))))
    # Projective chart division, followed by the nilpotent rescaling epsilon -> s epsilon.
    n=(z*(z*nn-z*nn[fi]*qbase)).applyfunc(s.cancel)
    w=(z*z*(z*ww-z*z*nn[fi]*nn+((z*nn[fi])**2-z*ww[fi])*qbase)).applyfunc(s.cancel)
    dn=80*(-23*L+10*V*z);dw=19200*V*(-23*L+10*V*z)**3
    assert all(not s.denom(s.cancel(v*dn)).free_symbols for v in n) and all(not s.denom(s.cancel(v*dw)).free_symbols for v in w)
    assert n[fi]==0 and w[fi]==0
    n0=n.subs(z,0);assert n0[24]==-31*L/460 and delta[0]==s.Rational(1,2)
    assert delta.row_join(n0)[[0,24],:].det()==-31*L/920
    variables=[K,L,V]+list(ns)+[z];R,*_=ring(variables,QQ);qp=[R.from_expr(v)for v in qbase];np=[R.from_expr(s.cancel(v*dn))for v in n];wp=[R.from_expr(s.cancel(v*dw))for v in w];DN=R.from_expr(dn);DW=R.from_expr(dw)
    one=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1];_,F=specialize(m,dict(zip(one,beta)),True)
    for f in F:
        values=[R.zero,R.zero,R.zero]
        for mon,c in f.items():
            if len(mon)!=2:continue
            i,j=mon;C=R.from_expr(c);values[0]+=C*qp[i]*qp[j];values[1]+=C*(qp[i]*np[j]+np[i]*qp[j]);values[2]+=C*((qp[i]*wp[j]+wp[i]*qp[j])*DN**2+np[i]*np[j]*DW)
        assert all(v==0 for v in values)
    print('universal flat3 identities',time.time()-start,flush=True)
    return dict(beta=list(map(str,beta)),source_beta=list(map(str,source)),q0=list(map(str,q0)),delta=list(map(str,delta)),n=list(map(str,n)),w=list(map(str,w)),n_den=str(dn),w_den=str(dw),n0=list(map(str,n0)),embedding_minor='-31*L/920',open_condition='L*V != 0',epsilon_order=3,all405_coefficient_identities=True,source_normalization=enc(To),finite_canonical_normalization=enc(Tn),beta_scale_homogeneities=[0,1,2],seconds=time.time()-start)

def source_geometry(m,univ):
    b=s.symbols('b0:10');e=b[0];beta=s.Matrix([0,0,0,0,2*e,b[1],b[7],b[2],b[3],b[8],b[4],b[5],b[6],e,-b[8]/3,b[9]])
    two=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==2];q0=s.Matrix(univ['q0']).applyfunc(s.sympify);fi=two.index(29);K,L,V=s.symbols('K L V');ns=s.symbols('n0:6')
    src=s.Matrix(univ['beta']).applyfunc(s.sympify);sol={K:(5*b[7]-b[9])/(3*e),L:4*b[8]/(9*e),V:b[7]/(3*e)}
    ss=s.solve(list((3*e*src.subs(sol)-beta).applyfunc(s.cancel)),ns,dict=True);assert len(ss)==1;sol.update(ss[0]);assert (3*e*src.subs(sol)-beta).applyfunc(s.cancel)==s.zeros(16,1)
    x,y,z,w=X[1:];vv=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=x*w;Q-=Q.T;orbit=[]
    for A in sp_basis():
        dq=s.Matrix(coords(m,sum(((A*vv)[i]*Q.diff(vv[i])for i in range(4)),s.zeros(4))-A*Q-Q*A.T))[two,:];dq-=q0*dq[fi];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*O.nullspace());assert O.rank()==5 and stab.shape==(10,5)
    u,act=beta_actions(m);H=beta.jacobian(b);ann=s.Matrix.hstack(*H.T.nullspace()).T;action=act.subs(dict(zip(u,beta)));N=(ann*action*stab).applyfunc(s.factor);assert N.to_DM().rank()==1
    sample=dict(zip(b,[3,2,3,5,7,11,13,2,3,5]));D=H.row_join(action);ds=D.subs(sample);pc=list(ds.to_DM().rref()[1]);assert len(pc)==16;det=ds[:,pc].det();joint=D.col_join(s.zeros(45,10).row_join(O));assert joint.subs(sample).rank()==16
    # Match the whole source tail, rather than merely one of its points.
    old=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());fam=next(f for f in old['families']if f['label']=='M1');tail=next(f for f in fam['rational_tails']if f['axis']==2 and f['k']=='-1/3')
    gg=s.symbols('g0:9');d=s.Symbol('d');U=s.Symbol('U');G=sum(s.sympify(tail['coordinate_formulas'].get(str(g),str(g)))*mon for g,mon in zip(gg,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3])).subs({gg[2]:d,U:0})
    assert s.expand(G-((V-3*L*d)*x*x*z-L*x*x*w+d*x*z*z))==0
    # All original U, not only the source section U=0, have endpoints in Z.
    source_all=s.Matrix([0,0,0,0,-s.Rational(1,3),ns[0],U,ns[1],ns[2],V,ns[3],ns[4],ns[5],s.Rational(1,3),K-5*V,5*U-L])
    aa=-3*L+14*U;Toall=s.Matrix([[1,0,0,0],[0,1,aa/2,0],[-aa/2,0,0,-1],[0,0,1,0]])
    assert Toall*J0[1:,1:]*Toall.T==J0[1:,1:]
    Ginf=x*z*z+aa*x*x*z;Pinf=s.zeros(4);Pinf[1,3]=s.diff(Ginf,z);Pinf-=Pinf.T
    assert (Toall.inv()*Pinf.subs(dict(zip(vv,Toall*vv)),simultaneous=True)*Toall.inv().T+2*Q).applyfunc(s.expand)==s.zeros(4)
    zz,uu,BB=primitive_beta(m);from weight1_structure import COORDS
    vec=s.Matrix([BB[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(uu);rr=list(basis.T.to_DM().rref()[1]);iv=basis[rr,:].inv()
    Bold=BB.subs(dict(zip(uu,source_all)));Bnew=Toall.T*Bold.subs(dict(zip(zz,Toall*s.Matrix(zz))),simultaneous=True)*Toall
    allbeta=(iv*s.Matrix([Bnew[i,j].diff(zz[k])for(i,j),k in COORDS])[rr,:]).applyfunc(s.factor)
    assert (ann*allbeta).applyfunc(s.cancel)==s.zeros(6,1)
    return dict(full_source_normalization=enc(Toall),full_source_normalized_beta=list(map(str,allbeta)),beta=list(map(str,beta)),parameters=list(map(str,b)),open_condition='b0*b7*b8 != 0',source_recovery={str(k):str(v)for k,v in sol.items()},source_beta_scale=str(3*e),projective_orbit_rank=5,stabilizer=enc(stab),source_annihilator=enc(ann),normal_stabilizer_matrix=enc(N),normal_stabilizer_rank=1,beta_sample={str(k):str(v)for k,v in sample.items()},beta_minor_columns=pc,beta_minor=str(det),joint_rank=16,image_dimension=16,covered_tails=['M1 g2 k=-1/3','Weyl partner M1 g4 k=-2/3'])

PRIMARY_PROGRAM='''
assert(dim(R/J)==3 and degree(R/J)==65);
B=saturate(J,ideal(v_8));D=saturate(J,ideal(v_7));
assert(D==actualD and intersect(B,D)==J and dim(R/D)==1 and degree(R/D)==3 and degree(R/B)==65);
C=ideal(v_9,v_7,v_6,v_5+(17/18)*v_8,v_4,v_2-(3/2)*v_8,v_1-v_8);
assert((gens B)%C==0);assert all(flatten entries gens C,g->g^4%B==0);
assert(dim(R/(J+ideal(v_0,v_3,v_8)))==0 and degree(R/(J+ideal(v_0,v_3,v_8)))==66);
S=ZZ/32003[w_1,w_2,w_4,w_5,w_6,w_7,w_9,w_0,w_3,MonomialOrder=>{GRevLex=>7,GRevLex=>2}];
ff=map(S,R,{w_0,w_1,w_2,w_3,w_4,w_5,w_6,w_7,1_S,w_9});GG=gens gb(ff B);cc={};
scan(flatten entries GG,g->(ee:=take((exponents leadTerm g)#0,7);mm:=product(7,j->S_j^(ee#j));cf:=(sum(select(terms g,tt->take((exponents tt)#0,7)==ee)))//mm;if first degree cf>0 then cc=append(cc,cf/(leadCoefficient cf))));
assert(cc=={});NI=ideal apply(flatten entries GG,g->product(7,j->S_j^((exponents leadTerm g)#0#j)));assert(dim(S/NI)==2 and degree(S/NI)==65);
print("B65 primary intersect actual flat triple line; transverse66");print("CERTIFICATE PASS");stdio<<flush;
'''


def run(m2_command=None):
    start=time.time();m=S61();univ=universal_thickening(m);geometry=source_geometry(m,univ);r=point_jets();pt=s.zeros(61,1)
    for k,v in r['sample'].items():pt[int(k[1:])-1]=s.Rational(v)
    pt[29]=1;local=calculate(m,True,pt);assert(local['leading_jacobian_rank'],local['homogenized_jacobian_rank'],local['quadratic_span_rank'],local['common_hessian_rank'])==(34,35,5,7)
    def mod(v):
        q=s.Rational(v);assert int(q.q)%PRIME;return int(q.p)*pow(int(q.q),-1,PRIME)%PRIME
    def rank(M):
        A=np.array([[mod(v)for v in row]for row in s.Matrix(M).tolist()],dtype=np.int64);j=0
        for c in range(A.shape[1]):
            ix=np.flatnonzero(A[j:,c])
            if not len(ix):continue
            a=j+int(ix[0]);A[[j,a]]=A[[a,j]];A[j]=A[j]*pow(int(A[j,c]),-1,PRIME)%PRIME;ii=np.flatnonzero(A[:,c]);ii=ii[ii!=j];A[ii]=(A[ii]-A[ii,c,None]*A[j,None,:])%PRIME;j+=1
            if j==len(A):break
        return j
    assert rank(r['J'][:,:44])==34 and rank(r['J'])==35 and rank(local['quadratic_obstruction_basis'])==5 and rank(r['curve_jacobian'])==2
    hs=[]
    for row in s.Matrix(local['quadratic_obstruction_basis']).tolist():
        H=s.zeros(10)
        for c,(i,j)in zip(row,local['quadratic_monomials']):H[i,j]+=s.Rational(c);H[j,i]+=s.Rational(c)
        hs+=H.tolist()
    assert rank(s.Matrix(hs))==7 and s.Matrix(hs)*r['family']==s.zeros(50,3) and mod(r['family'][r['parameters'],:].det())
    env={s.Symbol(k):s.sympify(v).subs({s.Symbol(a):s.Rational(b)for a,b in geometry['beta_sample'].items()})for k,v in geometry['source_recovery'].items()};z=s.Symbol('s')
    delta=s.Matrix(univ['delta']).applyfunc(s.sympify).subs(env);n0=s.Matrix(univ['n0']).applyfunc(s.sympify).subs(env)
    frame=delta[r['chart'],:].col_join(s.zeros(1,1)).row_join(n0[r['chart'],:].col_join(s.zeros(1,1)));fc=frame[r['free'],:];assert r['T']*fc==frame and fc.rank()==2 and fc[3,0]==0 and mod(fc[3,1])!=0
    ann=s.Matrix.hstack(*fc.T.nullspace()).T
    for label in ['n','w']:
        for value in univ[label]:
            num,den=s.fraction(s.cancel(s.sympify(value).subs(env)));assert mod(den.subs(z,0))!=0
            for c in s.Poly(num,z).all_coeffs()+s.Poly(den,z).all_coeffs():mod(c)
    with tempfile.TemporaryDirectory(prefix='poisson-triple-endpoint-') as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';basis=tmp/'basis.txt';binary=tmp/'closure'
        info=write_input(inp,PRIME,10,{d:A.tolist()for d,A in r['A'].items()},r['tseries'],maximum_degree=7)
        cc=subprocess.run(['c++','-O3','-std=c++17','-Wno-unqualified-std-cast-call',str(ROOT/'code/truncated_local_ideal_20261101.cpp'),'-o',str(binary)],capture_output=True,text=True,timeout=30);assert cc.returncode==0,cc.stderr
        cr=subprocess.run([str(binary),str(inp),str(basis)],capture_output=True,text=True,timeout=240);assert cr.returncode==0 and 'rank 16253 length 3195' in cr.stdout,cr.stdout[:2000]+cr.stderr
        allinitial=read_initial(basis,10,7);selected=[0,1,2,3,4,10,11,13,14,15,16,19,21,22,23,24,26,31,33,35,489,490,6047];chosen=[allinitial[i]for i in selected]
        assert [sum(len(next(iter(f)))==d for f in chosen)for d in [2,3,4,5,6,7]]==[5,15,0,2,0,1]
        expr=lambda f:'+'.join(str(c)+'*'+'*'.join('v_'+str(i)for i in mon)for mon,c in f.items())
        linear=lambda row:'+'.join(str(mod(v))+'*v_'+str(i)for i,v in enumerate(row)if v)
        program='R=ZZ/32003[v_0..v_9];J=ideal('+','.join(expr(f)for f in chosen)+');\nactualD=ideal('+','.join(linear(row)for row in ann.tolist())+',v_3^3);\n'+PRIMARY_PROGRAM
        source=tmp/'certificate.m2';source.write_text(program);envp=os.environ.copy();local_binary=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary')
        if m2_command:cmd=shlex.split(m2_command)
        elif local_binary.exists():cmd=[str(local_binary)];envp['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        else:cmd=[shutil.which('M2')or'M2']
        mr=subprocess.run(cmd+['--script',str(source)],env=envp,capture_output=True,text=True,timeout=60);output=mr.stdout+mr.stderr;assert mr.returncode==0 and 'CERTIFICATE PASS'in output and 'error:'not in output,output
    prime=json.loads((ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json').read_text());assert prime['all_passed']and prime['prime']==PRIME
    deps=['code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/s61_minimal_baseline_thickenings_20261101.py','code/s61_minimal_rankone_A1_removal_20261101.py','code/s61_minimal_cubic_tail_20261101.py','code/s61_minimal_infinity_local_data_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py','code/truncated_local_ideal_20261101.cpp','code/truncated_local_ideal_io_20261101.py','reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json','reports/S61_MINIMAL_CUBIC_TAIL_20261101.json','reports/S61_BINARY_PRIME_SPREADING_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md','research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md']
    paths=[Path(__file__)]+[ROOT/f for f in deps]
    return dict(id='S61_MINIMAL_TRIPLE_ENDPOINT_20261101',all_passed=True,prime=PRIME,sample=r['sample'],universal_thickening=univ,source_geometry=geometry,point_local_data=local,curve_jacobian=enc(r['curve_jacobian']),actual_binary_family_tangent=enc(r['family']),binary_parameters=r['parameters'],actual_triple_line_tangent=enc(fc),actual_triple_line_initial_annihilator=enc(ann),truncated_input=info,truncated_closure_output=cr.stdout.splitlines()[:9],initial_selected_indices=selected,initial_generators=[[[list(mon),int(c)]for mon,c in f.items()]for f in chosen],initial_generator_degrees=[len(next(iter(f)))for f in chosen],macaulay2_program=program,macaulay2_output=output,initial_decomposition='B65 intersect actual flat triple line',transverse_initial_length=66,conclusion='Specified ambient t^2=0; the entire sufficiently general-beta fiber of the dominant16D triple-tail endpoint is removed.',proof_note='research/notes/S61_MINIMAL_TRIPLE_ENDPOINT_20261101.md',limitations=['The degree7 relation is essential: degree5 and degree6 leave a spurious plane.','Uses a flat triple-line quotient and good-reduction/DVR proof, not finite-field sample extrapolation.','Other flag/binary endpoints and full affine exhaustion remain unresolved.'],input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in paths},seconds=time.time()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--m2-command');args=parser.parse_args();report=run(args.m2_command)
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
