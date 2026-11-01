#!/usr/bin/env python3
"""An M1 resonance tail with completed local ring C[[u,e]]/(e^3), t=0."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
ROOT=Path(__file__).resolve().parents[1]
def enc(A):return [[str(v) for v in row] for row in A.tolist()]
def tensor(a,b,d):
    x,y,z,w=X[1:];Q=s.zeros(4);Q[0,1]=x*x/2;Q[1,2]=-x*z-b*x*x;Q[1,3]=a*x*x+2*d*x*z;Q[2,3]=x*x/2
    return Q-Q.T
def coords(m,Q):
    P=X[0]**2*J;P[1:,1:]+=Q;return m.normalize(P)['parameters']

def universal(m):
    Q=tensor(1,0,0);p=coords(m,Q);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];K=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert K.shape==(16,10)
    b=s.symbols('b0:10');two,polys=specialize(m,dict(zip(one,K*s.Matrix(b))),True);q=s.Matrix([p[i] for i in two]);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];D=s.zeros(135,45)
    for row,f in enumerate(polys):
        val=0
        for mon,c in f.items():
            if len(mon)!=2:continue
            i,j=mon;val+=c*q[i]*q[j];D[row,i]+=c*q[j];D[row,j]+=c*q[i]
        assert s.expand(val)==0
    D=D[:,cols];D0=D[:70,:];cpiv=list(D0.rref()[1]);rows=list(D0[:,cpiv].T.rref()[1]);C=s.Matrix.hstack(*D0.nullspace());R=D[70:,:]*C;rr,piv=R.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(R.cols) if i not in piv];N=s.zeros(R.cols,len(free))
    for j,i in enumerate(free):
        N[i,j]=1
        for a,pc in enumerate(piv):N[pc,j]=-rr[a,i]
    assert len(cpiv)==25 and len(cpiv)+len(piv)==42 and len(free)==2
    assert (R.to_DM()*N.to_DM()).is_zero_matrix;T=C*N
    U=s.zeros(44,len(cpiv));iv=D0.extract(rows,cpiv).inv()
    for j,i in enumerate(cpiv):U[i,:]=iv[j,:]
    H=s.zeros(135,3);tf={a:T[j,:] for j,a in enumerate(cols)};tf[fixed]=s.zeros(1,2)
    for row,f in enumerate(polys):
        for mon,c in f.items():
            if len(mon)!=2:continue
            i,j=mon;H[row,0]+=c*tf[i][0]*tf[j][0];H[row,1]+=c*(tf[i][0]*tf[j][1]+tf[i][1]*tf[j][0]);H[row,2]+=c*tf[i][1]*tf[j][1]
    H=H.applyfunc(s.cancel);assert (H[:70,:]-D0*U*H[rows,:]).applyfunc(s.cancel)==s.zeros(70,3)
    Z=(H[70:,:]-D[70:,:]*U*H[rows,:]).applyfunc(s.cancel);ar,apiv=R.row_join(Z).to_DM().rref();assert len(apiv)==len(piv);ar=ar.to_Matrix();W=s.zeros(R.cols,3)
    for j in range(3):
        for row,pc in enumerate(piv):W[pc,j]=ar[row,R.cols+j]
    lift=(U*H[rows,:]+C*W).applyfunc(s.cancel);assert (D.to_DM()*lift.to_DM()-H.to_DM()).is_zero_matrix
    assert (D.to_DM()*T.to_DM()).is_zero_matrix
    return dict(beta_kernel=enc(K),q_point=list(map(str,q)),chart_columns=cols,leading_jacobian=enc(D),first_rank=25,first_kernel=enc(C),residual_matrix=enc(R),residual_kernel=enc(N),free_columns=free,tangent=enc(T),quadratic_matrix=enc(H),quadratic_lift=enc(lift),leading_rank=42,quadratic_obstruction_zero=True)

def local(m):
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];p=coords(m,tensor(16,3,1));bv=[0,0,0,0,-s.Rational(1,3),2,2,3,5,3,7,11,13,s.Rational(1,3),4,5]
    for i,a in zip(one,bv):p[i]=s.sympify(a)
    r=calculate(m,True,p);assert r['leading_jacobian_rank']==42 and r['homogenized_jacobian_rank']==43 and r['quadratic_span_rank']==0
    _,F=specialize(m,dict(zip(one,bv)),True);q=s.Matrix([p[i] for i in two]);cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols);Jac=s.Matrix(r['jacobian']).applyfunc(s.Rational);T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational)
    z=s.Symbol('z');family=coords(m,tensor(3+13*z,3,z));qq=s.Matrix([family[i] for i in two])
    assert all(s.expand(sum(v*s.prod(qq[i] for i in mon) for mon,v in f.items() if len(mon)==2))==0 for f in F)
    dq=qq.diff(z).subs(z,1);dq-=q*dq[fixed]/q[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));rT=list(T.T.rref()[1]);vt=T[rT,:].inv()*tv[rT,:];assert T*vt==tv
    uindex=next(i for i,a in enumerate(vt) if a);normal=1-uindex;rows=list(Jac.T.to_DM().rref()[1]);piv=list(Jac.to_DM().rref()[1]);U=s.zeros(45,43);iv=Jac.extract(rows,piv).inv(method='DM')
    for i,col in enumerate(piv):U[col,:]=iv[i,:]
    assert Jac[rows,:]*U==s.eye(43);JU=Jac*U;terms=[]
    for f in F:
        ts=[]
        for mon,c in f.items():
            if fixed in mon:continue
            i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
        terms.append(ts)
    series={1:T[:,normal]};out=[]
    for d in [2,3]:
        raw=s.zeros(135,1)
        for row,ts in enumerate(terms):
            for i,j,c in ts:raw[row]+=c*sum(series[a][i]*series[d-a][j] for a in range(1,d))
        res=raw-JU*raw[rows,:];wd=-U*raw[rows,:];series[d]=wd
        if d==2:assert res==s.zeros(135,1)
        else:assert res!=s.zeros(135,1)
        assert wd[-1]==0
        out.append(dict(degree=d,raw=list(map(str,raw)),residual=list(map(str,res)),lift=list(map(str,wd))))
    return dict(local_data=r,whole_tail_leading_zero=True,tail_q_coordinates=list(map(str,qq)),tail_parameter='z',tail_tangent_in_kernel=list(map(str,vt)),parameter_kernel_coordinate=uindex,normal_kernel_coordinate=normal,jacobian_pivot_rows=rows,jacobian_pivot_columns=piv,transverse_expansion=out,transverse_length=3)

def normalization():
    a,b,d=s.symbols('a b d');x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);T1=s.eye(4);T1[3,2]=-2*d;T2=s.eye(4);T2[1,3]=b;T2[2,0]=-b;T3=s.diag(1,1,1/(a-2*b*d),a-2*b*d);T=T1*T2*T3
    assert (T*J[1:,1:]*T.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
    result=T.inv()*tensor(a,b,d).subs(dict(zip(V,T*V)),simultaneous=True)*T.inv().T
    assert (result-tensor(1,0,0)).applyfunc(s.cancel)==s.zeros(4)
    return dict(matrix=enc(T),open_condition='a-2*b*d != 0',tail_value_at_point='10',verified=True)

def run():
    start=time.time();m=S61();univ=universal(m);point=local(m);change=normalization()
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_CUBIC_TAIL_20261101',all_passed=True,field='QQ and QQ(b0,...,b9)',universal=univ,point=point,normalization=change,conclusion='At the specified M1 resonance-tail point and on the justified open, the completed full135 homogenized local ring is C[[u,e]]/(e^3), with t=0.',limitations=['Does not remove exceptional points or endpoints of the full tail family.','Does not prove full affine-fiber finiteness or exhaustion.'],proof_note='research/notes/S61_MINIMAL_CUBIC_TAIL_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,transverse_length=3,seconds=r['seconds'])))
if __name__=='__main__':main()
