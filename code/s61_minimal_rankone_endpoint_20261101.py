#!/usr/bin/env python3
"""Remove the rank-one A1 binary endpoint using good reduction."""
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
    sample=dict(zip(('u4','u5','u12','u13','u14','u25','u26','u27','u34','u35','u36','u37','u48','u49','u50','u51'),map(str,(0,3,0,0,2,2,-s.Rational(39,2),3,5,9,7,11,13,4,45,5))));expected_parameters=[1,3,4]
    start=time.time();p=32003;m=S61()
    r={'sample':dict(sample)}
    b={int(k[1:])-1:s.Rational(v) for k,v in r['sample'].items()};two,F=specialize(m,b,True);x,y,zz,w=X[1:];Q0=s.zeros(4);Q0[1,2]=w*w+15*x*x;Q0-=Q0.T;q0=s.Matrix(coords(m,Q0))[two,:];fixed=next(i for i,v in enumerate(q0)if v);cols=[i for i in range(45)if i!=fixed]
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
    x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell2**2+15*ell1**2+h*ell1*ell2+k*ell1**2)*(v1*v2.T-v2*v1.T)
    ds=[]
    for v in [a,c,d,h,k]:
     P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);qs=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);qs-=q0*qs[fixed]/q0[fixed];ds.append(qs[cols,:].col_join(s.zeros(1,1)))
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
    H%=p;n=len(normal);W={1:t[:,normal]};A={};ms={d:list(combinations_with_replacement(range(n),d)) for d in range(1,7)};idx={d:{m:i for i,m in enumerate(ms[d])} for d in range(1,7)}
    for degree in range(2,7):
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
    return dict(sample=sample,q0=q0,J=J,T=T,F=F,two=two,fixed=fixed,chart=cols,free=free,rows=R,columns=C,parameters=params,family=coordinates,curve_jacobian=curve_jac,A=A,tseries={d:W[d][-1,:].tolist() for d in A},monomials=ms)


def universal_double_line(m):
    """Transport the canonical A1 tangent over a quadratic extension, then descend."""
    from weight1_structure import COORDS
    from s61_minimal_rankone_junction_20261101 import tensor
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ
    st=time.time();old=json.loads((ROOT/'reports/S61_MINIMAL_RANKONE_JUNCTION_20261101.json').read_text())['A1_universal'];BK=s.Matrix(old['beta_kernel']);b=s.symbols('b0:9');h,la,z=s.symbols('h la s');source=BK*s.Matrix(b);zz,u,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(zz[k])for(i,j),k in COORDS]);M=vec.jacobian(u);rows=list(M.T.to_DM().rref()[1]);iv=M[rows,:].inv()
    dh=1+b[7]*(h-1)/(2*b[6]);shape=2*b[7]**2*(h-1)/(9*b[6]**2);ro=h*la**3;sh=s.Matrix([[1,0,0,0],[0,1,0,-shape],[shape,0,1,0],[0,0,0,1]]);Tn=sh*s.diag(ro,1/ro,1/la,la);hh=b[7]/(18*b[6]);aa=3*hh;To=s.Matrix([[1/hh,0,0,0],[0,hh,-aa,0],[0,0,1,0],[-aa/hh,0,0,1]])
    assert(Tn*J0[1:,1:]*Tn.T-J0[1:,1:]).applyfunc(s.cancel)==s.zeros(4)and(To*J0[1:,1:]*To.T-J0[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
    def even(f):
        f=s.cancel(f)
        if not f:return s.S.Zero
        nn,dd=s.fraction(f);pn=s.Poly(nn,la);pd=s.Poly(dd,la);kn=pn.monoms()[0][0]%2;kd=pd.monoms()[0][0]%2;assert kn==kd
        assert all(k%2==kn for(k,)in pn.monoms())and all(k%2==kd for(k,)in pd.monoms())
        an=sum(v*dh**(-(k-kn)//2)for(k,),v in pn.terms());ad=sum(v*dh**(-(k-kd)//2)for(k,),v in pd.terms());return s.factor(s.cancel(an/ad))
    BB=B.subs(dict(zip(u,source)))
    def pullbeta(T):
        P=T.T*BB.subs(dict(zip(zz,T*s.Matrix(zz))),simultaneous=True)*T;return iv*s.Matrix([P[i,j].diff(zz[k])for(i,j),k in COORDS])[rows,:]
    canonical=(pullbeta(Tn)/la).applyfunc(even);beta=pullbeta(To).applyfunc(s.factor);br=list(BK.T.to_DM().rref()[1]);bc=(BK[br,:].inv()*canonical[br,:]).applyfunc(s.factor);assert(BK*bc-canonical).applyfunc(s.factor)==s.zeros(16,1)
    xv=s.Matrix(X[1:]);Q=tensor(0,shape,dh,h);test=(Tn.inv()*Q.subs(dict(zip(xv,Tn*xv)),simultaneous=True)*Tn.inv().T/ro**2).applyfunc(even);assert(test-tensor(0,0,1,1)).applyfunc(s.factor)==s.zeros(4)
    clean=(s.Matrix(old['first_kernel'])*s.Matrix(old['residual_kernel'])).applyfunc(s.factor);nc=clean[:,1].applyfunc(lambda v:s.factor(v.subs(dict(zip(b,bc)),simultaneous=True)));S=(Tn.inv()*To).applyfunc(s.factor);Si=S.inv().applyfunc(s.factor);two=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==2];out=s.zeros(45,1)
    for j,i in enumerate(old['chart_columns']):
        if nc[j]==0:continue
        par=[s.S.Zero]*61;par[two[i]]=1;P=m.reconstruct(par);D=s.zeros(4)
        for a in range(4):
            for k in range(4):D[a,k]=s.expand((P[a+1,k+1]-X[a+1]*P[0,k+1]+X[k+1]*P[0,a+1]).subs(X[0],1)-J0[a+1,k+1])
        assert s.Matrix(coords(m,D))[two,:]==s.eye(45)[:,i]
        tr=(ro**2*Si*D.subs(dict(zip(xv,S*xv)),simultaneous=True)*Si.T).applyfunc(s.expand);qp=s.Matrix(coords(m,tr))[two,:].applyfunc(even);out+=qp*nc[j]
    out=out.applyfunc(s.factor);Qt=(To.inv()*Q.subs(dict(zip(xv,To*xv)),simultaneous=True)*To.inv().T).applyfunc(s.factor);q=s.Matrix(coords(m,Qt))[two,:];Q0=s.zeros(4);Q0[1,2]=X[4]**2+15*X[1]**2;Q0-=Q0.T;q0=s.Matrix(coords(m,Q0))[two,:];fi=next(i for i,v in enumerate(q0)if v);kk=s.factor(q.diff(h)[fi]/q0[fi]);assert(q.diff(h)-kk*q0).applyfunc(s.factor)==s.zeros(45,1)
    qc=q.subs(h,0);delta=(qc-q0*qc[fi]/q0[fi]).applyfunc(s.factor);hs=s.factor((q0[fi]/z-qc[fi])/(kk*q0[fi]));qb=q0+z*delta;assert(z*q.subs(h,hs)-qb).applyfunc(s.factor)==s.zeros(45,1);nn=out.applyfunc(lambda v:s.factor(v.subs(h,hs)));n=(z*nn-z*nn[fi]*qb/q0[fi]).applyfunc(s.cancel);n0=n.subs(z,0).applyfunc(s.factor);minor=s.factor(delta.row_join(n0)[[0,34],:].det());assert minor==2*b[6]**3*(216*b[3]*b[6]-73*b[7]**2)/(b[7]**4*(2*b[6]-b[7]))
    dq=5*b[7]**2;dn=24300*b[3]*b[7]**7*(2*b[6]-b[7])*(1188*b[6]**2*z-162*b[6]*b[7]*z-5*b[7]**2)**2;db=96*b[6]*b[7]**2
    R,*_=ring(list(b)+[z],QQ);qp=[R.from_expr(s.cancel(v*dq))for v in qb];np=[R.from_expr(s.cancel(v*dn))for v in n];one=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1];index={i:j for j,i in enumerate(two)};bp={i:R.from_expr(s.cancel(v*db))for i,v in zip(one,beta)};DB=R.from_expr(db);eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    for f in eqs:
        values=[R.zero,R.zero]
        for mon,nu,de in f:
            inds=[index[i]for i in mon if i in index]
            if len(inds)!=2:continue
            bs=[i for i in mon if i in bp];assert len(bs)<=2;C=R(nu)/de*DB**(2-len(bs))
            for i in bs:C*=bp[i]
            i,j=inds;values[0]+=C*qp[i]*qp[j];values[1]+=C*(qp[i]*np[j]+np[i]*qp[j])
        assert not any(values)
    print('universal double line',time.time()-st,flush=True)
    return dict(source_beta=list(map(str,source)),beta=list(map(str,beta)),q0=list(map(str,q0)),delta=list(map(str,delta)),n=list(map(str,n)),n0=list(map(str,n0)),qbase=list(map(str,qb)),n_den=str(dn),q_den=str(dq),beta_den=str(db),embedding_minor=str(minor),open_condition='b3*b6*b7*(2*b6-b7)*(216*b3*b6-73*b7^2) != 0',finite_normalization=enc(Tn),endpoint_normalization=enc(To),quadratic_extension='la^2=1/dh',h_substitution=str(hs),all270_coefficient_identities=True,epsilon_order=2,seconds=time.time()-st)
def source_geometry(m,univ):
    """H9 fills a dominant16D locus; the finite A1 normalization covers its endpoints."""
    from weight1_structure import minor_witness
    a,p=s.symbols('a p');v=s.symbols('c0:7');params=[a,p]+list(v);beta=s.Matrix([0,p,0,0,a,v[0],3*p/2-12*a,v[1],v[2],9*a/2,v[3],v[4],v[5],2*a,45*a/2,v[6]]);H=beta.jacobian(params);ann=s.Matrix.hstack(*H.T.nullspace()).T;assert H.rank()==9
    b=s.symbols('b0:9');nv=s.Matrix(univ['beta']);sub={b[7]:36*a,b[6]:54*a*a/(p+3*a)};unknown=list(b[:6])+[b[8]];solutions=s.solve(list((nv.subs(sub)-beta).applyfunc(s.factor)),unknown);assert len(solutions)==7;inverse=sub|solutions;assert(nv.subs(inverse)-beta).applyfunc(s.factor)==s.zeros(16,1)
    two=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==2];q0=s.Matrix(univ['q0']);fi=next(i for i,vv in enumerate(q0)if vv);x,y,z,w=X[1:];Q=s.zeros(4);Q[1,2]=w*w+15*x*x;Q-=Q.T;xv=s.Matrix(X[1:]);orbit=[]
    for A in sp_basis():
        D=sum(((A*xv)[i]*Q.diff(xv[i])for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=s.Matrix(coords(m,D))[two,:];dq-=q0*dq[fi]/q0[fi];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);assert O.rank()==5;stab=s.Matrix.hstack(*O.nullspace());u,act=beta_actions(m);AB=act.subs(dict(zip(u,beta)));normal=(ann*AB*stab).applyfunc(s.factor);assert normal.to_DM().rank()==2
    sample={a:2,p:3,**dict(zip(v,[2,3,5,7,11,13,5]))};D=H.row_join(AB);wit=minor_witness(D.subs(sample));assert wit['rank']==16;join=D.col_join(s.zeros(45,9).row_join(O));assert join.subs(sample).rank()==16
    return dict(beta_kernel=enc(H),beta=list(map(str,beta)),parameter_order=list(map(str,params)),source_inverse={str(k):str(vv)for k,vv in inverse.items()},source_inverse_open='a*(p+3*a) != 0',projective_orbit=enc(O),projective_stabilizer=enc(stab),stabilizer_normal_matrix=enc(normal),stabilizer_normal_rank=2,projective_orbit_dimension=5,beta_differential=enc(D),beta_minor=wit,joint_rank=16,closure_dimension=16,endpoint_coverage='Every finite A1 curve has a point with h*(d-3ha) nonzero. The existing finite normalization gives tensor(0,0,1,1); its same-beta A1 curve has the stated endpoint. The old b3,b6,b7,2b6-b7 bad-image bounds remain dependencies.')

PRIMARY_PROGRAM='''
B=saturate(J,ideal(v_6));D=saturate(J,ideal(v_4-15*v_6));
assert(dim(R/J)==3 and degree(R/J)==65);assert(dim(R/B)==3 and degree(R/B)==65);assert(dim(R/D)==1 and degree(R/D)==2);
assert(intersect(B,D)==J);assert(D==actualD);assert(dim(R/(J+ideal(v_1,v_3,v_4)))==0 and degree(R/(J+ideal(v_1,v_3,v_4)))==66);
C=ideal(v_9,v_8,v_7,v_5,v_4-15*v_6,v_2,v_0);assert((gens B)%C==0);assert all(flatten entries gens C,g->g^4%B==0);
S=ZZ/32003[w_0,w_2,w_4,w_5,w_7,w_8,w_9,w_1,w_3,MonomialOrder=>{GRevLex=>7,GRevLex=>2}];ff=map(S,R,{w_0,w_1,w_2,w_3,w_4,w_5,1_S,w_7,w_8,w_9});GG=gens gb(ff B);
scan(flatten entries GG,g->(ee:=take((exponents leadTerm g)#0,7);mm:=product(7,j->S_j^(ee#j));cf:=(sum(select(terms g,tt->take((exponents tt)#0,7)==ee)))//mm;assert(first degree cf==0)));
NI=ideal apply(flatten entries GG,g->product(7,j->S_j^((exponents leadTerm g)#0#j)));assert(dim(S/NI)==2 and degree(S/NI)==65);
print("J = B65 intersect actual double line");print("B: constant leading coefficients, 65 normal standard monomials");print("transverse length66");print("CERTIFICATE PASS");stdio<<flush;
'''

def run(m2_command=None):
    start=time.time();m=S61();univ=universal_double_line(m);geometry=source_geometry(m,univ);r=point_jets();one=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1];point=s.zeros(61,1)
    for i in one:point[i]=s.Rational(r['sample']['u'+str(i+1)])
    for i,q in zip(r['two'],r['q0']):point[i]=q
    local=calculate(m,True,point);assert local['leading_jacobian_rank']==34 and local['homogenized_jacobian_rank']==35 and local['quadratic_span_rank']==5 and local['common_hessian_rank']==7
    b=s.symbols('b0:9');z=s.Symbol('s');env=dict(zip(b,[s.Rational(11,36),s.Rational(1,2),-s.Rational(605,72),-s.Rational(155,6),-s.Rational(119,12),s.Rational(341,12),24,72,s.Rational(55,2)]));assert s.Matrix(univ['beta']).subs(env)==point[one,:];delta=s.Matrix(univ['delta']).subs(env).applyfunc(s.cancel);n0=s.Matrix(univ['n0']).subs(env).applyfunc(s.cancel);curve=delta.row_join(n0)[r['chart'],:].col_join(s.zeros(1,2));fc=curve[r['free'],:];assert r['T']*fc==curve and fc.rank()==2;assert fc[6,0]==0 and fc[6,1]!=0;ann=s.Matrix.hstack(*fc.T.nullspace()).T
    def mod(v):
        v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,PRIME)%PRIME
    def rank(M):
        A=np.array([[mod(v)for v in row]for row in s.Matrix(M).tolist()],dtype=np.int64);j=0
        for c in range(A.shape[1]):
            ix=np.flatnonzero(A[j:,c])
            if not len(ix):continue
            a=j+int(ix[0]);A[[j,a]]=A[[a,j]];A[j]=A[j]*pow(int(A[j,c]),-1,PRIME)%PRIME;ii=np.flatnonzero(A[:,c]);ii=ii[ii!=j];A[ii]=(A[ii]-A[ii,c,None]*A[j,None,:])%PRIME;j+=1
            if j==len(A):break
        return j
    assert rank(r['J'])==35 and rank(r['J'][:,:44])==34 and rank(local['quadratic_obstruction_basis'])==5 and rank(r['curve_jacobian'])==2
    hs=[]
    for row in s.Matrix(local['quadratic_obstruction_basis']).tolist():
        H=s.zeros(10)
        for c,(i,j)in zip(row,local['quadratic_monomials']):H[i,j]+=s.Rational(c);H[j,i]+=s.Rational(c)
        hs+=H.tolist()
    assert rank(s.Matrix(hs))==7 and s.Matrix(hs)*r['family']==s.zeros(50,3) and mod(r['family'][r['parameters'],:].det())
    for value in univ['n']+univ['qbase']:
        num,den=s.fraction(s.cancel(s.sympify(value).subs(env)));assert mod(den.subs(z,0))!=0
        for c in s.Poly(num,z).all_coeffs()+s.Poly(den,z).all_coeffs():mod(c)
    assert mod(fc[[1,6],:].det())!=0
    with tempfile.TemporaryDirectory(prefix='poisson-rankone-endpoint-')as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';basis=tmp/'basis.txt';binary=tmp/'closure';info=write_input(inp,PRIME,10,{d:A.tolist()for d,A in r['A'].items()},r['tseries'],maximum_degree=6)
        cc=subprocess.run(['c++','-O3','-std=c++17','-Wno-unqualified-std-cast-call',str(ROOT/'code/truncated_local_ideal_20261101.cpp'),'-o',str(binary)],capture_output=True,text=True,timeout=30);assert cc.returncode==0,cc.stderr
        cr=subprocess.run([str(binary),str(inp),str(basis)],capture_output=True,text=True,timeout=120);assert cr.returncode==0 and 'rank 6042 length 1966'in cr.stdout,cr.stdout[:1000]+cr.stderr
        allinitial=read_initial(basis,10,6);selected=[0,1,2,3,4,5,7,9,10,11,12,13,21,22,23,24,26,31,33,35,697,698,2620];chosen=[allinitial[i]for i in selected];assert [sum(len(next(iter(f)))==d for f in chosen)for d in [2,3,4,5,6]]==[5,15,0,2,1]
        expr=lambda f:'+'.join(str(c)+'*'+'*'.join('v_'+str(i)for i in mon)for mon,c in f.items());linear=lambda row:'+'.join(str(mod(v))+'*v_'+str(i)for i,v in enumerate(row)if v)
        program='R=ZZ/32003[v_0..v_9];J=ideal('+','.join(expr(f)for f in chosen)+');\nactualD=ideal('+','.join(linear(row)for row in ann.tolist())+',v_6^2);\n'+PRIMARY_PROGRAM;source=tmp/'certificate.m2';source.write_text(program);en=os.environ.copy();lb=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary')
        if m2_command:cmd=shlex.split(m2_command)
        elif lb.exists():cmd=[str(lb)];en['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        else:cmd=[shutil.which('M2')or'M2']
        mr=subprocess.run(cmd+['--script',str(source)],env=en,capture_output=True,text=True,timeout=60);output=mr.stdout+mr.stderr;assert mr.returncode==0 and 'CERTIFICATE PASS'in output and 'error:'not in output,output
    prime=json.loads((ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json').read_text());assert prime['all_passed']and prime['prime']==PRIME
    deps=['code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/s61_minimal_baseline_thickenings_20261101.py','code/s61_minimal_rankone_A1_removal_20261101.py','code/s61_minimal_rankone_junction_20261101.py','code/s61_minimal_infinity_local_data_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py','code/truncated_local_ideal_20261101.cpp','code/truncated_local_ideal_io_20261101.py','reports/S61_MINIMAL_RANKONE_JUNCTION_20261101.json','reports/S61_MINIMAL_RANKONE_A1_REMOVAL_20261101.json','reports/S61_BINARY_PRIME_SPREADING_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md','research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md'];paths=[Path(__file__)]+[ROOT/f for f in deps]
    return dict(id='S61_MINIMAL_RANKONE_ENDPOINT_20261101',all_passed=True,prime=PRIME,sample=r['sample'],source_sample={str(k):str(v)for k,v in env.items()},universal_double_line=univ,source_geometry=geometry,point_local_data=local,curve_jacobian=enc(r['curve_jacobian']),actual_binary_family_tangent=enc(r['family']),binary_parameters=r['parameters'],actual_double_line_tangent=enc(fc),actual_double_line_initial_annihilator=enc(ann),truncated_input=info,truncated_closure_output=cr.stdout.splitlines()[:8],initial_selected_indices=selected,initial_generators=[[[list(mon),int(c)]for mon,c in f.items()]for f in chosen],initial_generator_degrees=[len(next(iter(f)))for f in chosen],macaulay2_program=program,macaulay2_output=output,initial_decomposition='B65 intersect actual flat double line',transverse_initial_length=66,conclusion='Specified ambient t^2=0; the entire sufficiently general-beta fiber of the dominant16D rank-one A1 endpoint is removed.',limitations=['The true sextic relation is essential; degree5 leaves an extra plane.','Good reduction and a flat nonreduced quotient are essential; no extrapolation from one finite-field sample.','Other endpoints, flags, deeper binary, and full affine exhaustion remain separate.'],proof_note='research/notes/S61_MINIMAL_RANKONE_ENDPOINT_20261101.md',input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in paths},seconds=time.time()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--m2-command');args=parser.parse_args();report=run(args.m2_command)
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
