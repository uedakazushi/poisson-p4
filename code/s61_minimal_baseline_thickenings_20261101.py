#!/usr/bin/env python3
"""Local square-zero rings on the corrected free baseline resonance lines."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
ROOT=Path(__file__).resolve().parents[1]
def enc(A):return [[str(v) for v in row] for row in A.tolist()]
def tensor(kind,a,b,c):
    x,y,z,w=X[1:];Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-b*x*x-c*x*z-(x*z if kind=='M1' else 0);Q[1,3]=a*x*x+c*x*w
    return Q-Q.T
def coords(m,Q):
    P=X[0]**2*J;P[1:,1:]+=Q;return m.normalize(P)['parameters']
def universal(m,kind):
    c=2 if kind=='M0' else 1;Q=tensor(kind,0,0,c);pt=coords(m,Q);BK=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert BK.shape==(16,10);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];b=s.symbols('b0:10');two,F=specialize(m,dict(zip(one,BK*s.Matrix(b))),True);q=s.Matrix([pt[i] for i in two]);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];A=s.zeros(135,45);tc=s.zeros(135,1)
    for row,f in enumerate(F):
        val=0
        for mon,v in f.items():
            if len(mon)==2:
                i,j=mon;val+=v*q[i]*q[j];A[row,i]+=v*q[j];A[row,j]+=v*q[i]
            elif len(mon)==1:tc[row]+=v*q[mon[0]]
        assert s.expand(val)==0
    A=A[:,cols];A0=A[:70,:];cc=list(A0.rref()[1]);rr=list(A0[:,cc].T.rref()[1]);C=s.Matrix.hstack(*A0.nullspace());U=s.zeros(44,len(cc));iv=A0[rr,cc].inv()
    for j,i in enumerate(cc):U[i,:]=iv[j,:]
    R=(A*C).applyfunc(s.cancel);ct=(tc-A*U*tc[rr,:]).applyfunc(s.cancel);full=R.row_join(ct);rf,piv=full.to_DM().rref();rf=rf.to_Matrix();free=[i for i in range(full.cols) if i not in piv];N=s.zeros(full.cols,len(free))
    for j,i in enumerate(free):
        N[i,j]=1
        for r,pc in enumerate(piv):N[pc,j]=-rf[r,i]
    assert (full.to_DM()*N.to_DM()).is_zero_matrix;_,lpiv=R.to_DM().rref();actual=(len(cc)+len(lpiv),len(cc)+len(piv),N[-1,:]==s.zeros(1,N.cols));assert actual==((41,41,False) if kind=='M0' else (41,42,True))
    # M1's ten-variable kernel introduces a genuine top70 compatibility
    # obstruction to a t direction; do not discard these residual rows.
    assert (ct[:70,:]==s.zeros(70,1))==(kind=='M0')
    return dict(beta_kernel=enc(BK),q_point=list(map(str,q)),chart_columns=cols,leading_jacobian=enc(A),t_column=list(map(str,tc)),first_rank=len(cc),first_pivot_rows=rr,first_pivot_columns=cc,first_kernel=enc(C),residual_matrix=enc(R),residual_t_column=list(map(str,ct)),residual_kernel=enc(N),free_columns=free,leading_rank=actual[0],full_rank=actual[1],homogenizer_zero_on_tangent=actual[2],top70_t_compatibility_zero=kind=='M0')
def local(m,kind):
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];k=-s.Rational(2,3) if kind=='M0' else -s.Rational(4,9);z=s.Symbol('z');Q=tensor(kind,20,z,2) if kind=='M0' else tensor(kind,z,-1,1);family=coords(m,Q);pt=[s.sympify(v).subs(z,1) for v in family];bv=[0,0,0,0,k,2,2,3,5,3,7,11,13,2*k+1,4,5]
    for i,v in zip(one,bv):pt[i]=s.sympify(v)
    _,F=specialize(m,dict(zip(one,bv)),True);qq=s.Matrix([family[i] for i in two]);assert all(s.expand(sum(v*s.prod(qq[i] for i in mon) for mon,v in f.items() if len(mon)==2))==0 for f in F)
    data=calculate(m,True,pt);n=4 if kind=='M0' else 3;assert data['homogenized_jacobian_rank']==45-n and data['leading_jacobian_rank']==41 and data['quadratic_span_rank']==n*(n-1)//2 and data['common_hessian_rank']==n-1
    T=s.Matrix(data['tangent_basis']).applyfunc(s.Rational);rad=s.Matrix(data['quadratic_common_radical']).applyfunc(s.Rational);q=qq.subs(z,1);cols=data['chart_columns'];fixed=next(i for i in range(45) if i not in cols);dq=qq.diff(z).subs(z,1);dq-=q*dq[fixed]/q[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));rows=list(T.T.rref()[1]);vt=T[rows,:].inv()*tv[rows,:];assert T*vt==tv and s.Matrix.hstack(vt,rad).rank()==1 and vt!=s.zeros(n,1)
    assert data['homogenizer_zero_on_tangent']==(kind=='M1')
    return dict(local_data=data,tail_q_coordinates=list(map(str,qq)),tail_tangent_in_kernel=list(map(str,vt)),whole_tail_leading_zero=True,transverse_length_upper_bound=n)
def normalization(kind):
    a,b=s.symbols('a b');c=2 if kind=='M0' else 1;A=-a/(c+1) if kind=='M0' else -a/c;B=-b/(c-1) if kind=='M0' else -b/(c+1);T=s.Matrix([[1,0,0,0],[0,1,A,-B],[B,0,1,0],[A,0,0,1]]);V=s.Matrix(X[1:]);Q=tensor(kind,a,b,c)
    assert T*J[1:,1:]*T.T==J[1:,1:];assert (T.inv()*Q.subs(dict(zip(V,T*V)),simultaneous=True)*T.inv().T-tensor(kind,0,0,c)).applyfunc(s.expand)==s.zeros(4)
    E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);paired=-c if kind=='M0' else -c-1;Qc=tensor(kind,0,0,c)
    assert E*J[1:,1:]*E.T==J[1:,1:];assert (E.inv()*Qc.subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T-tensor(kind,0,0,paired)).applyfunc(s.expand)==s.zeros(4)
    return dict(shear=enc(T),canonical_source_c=c,paired_source_c=paired,paired_exchange=enc(E),verified=True)
def run():
    start=time.time();m=S61();cases={}
    for kind in ['M0','M1']:
        cases[kind]=dict(universal=universal(m,kind),point=local(m,kind),normalization=normalization(kind));print(json.dumps({'case':kind,'passed':True}),flush=True)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_BASELINE_THICKENINGS_20261101',all_passed=True,field='QQ and QQ(b0,...,b9)',cases=cases,conclusion='Specified M0 completion is C[[u,e1,e2,e3]]/(e1,e2,e3)^2 with nonzero t and t^2=0; specified M1 completion is C[[u,e1,e2]]/(e1,e2)^2 with t=0.',limitations=['Specified points and the justified open on the paired baseline lines; not their intersection points or endpoints.','No full affine-fiber finiteness or exhaustion conclusion.'],proof_note='research/notes/S61_MINIMAL_BASELINE_THICKENINGS_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
if __name__=='__main__':main()
