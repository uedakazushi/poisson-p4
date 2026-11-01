#!/usr/bin/env python3
"""Two R0 resonance intersections: transverse length4 and homogenizer removal."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import sp_basis,COORDS,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]

def tail_certificate():
    m=S61();x,y,z,w=X[1:];G=x*z*w+z**3;Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);Q[2,1]=-Q[1,2];Q[3,1]=-Q[1,3]
    M=beta_matrix(m,Q);K=s.Matrix.hstack(*M.nullspace());assert K.cols==9
    b=s.symbols(f'b0:{K.cols}');betav=K*s.Matrix(b);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];beta=dict(zip(one,betav));two,F=specialize(m,beta,True)
    P=X[0]**2*J;P[1:,1:]+=Q;point=m.normalize(P)['parameters'];q0=s.Matrix([point[i] for i in two]);D=s.zeros(135,45)
    for row,f in enumerate(F):
     value=0
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon;value+=c*q0[i]*q0[j];D[row,i]+=c*q0[j];D[row,j]+=c*q0[i]
     assert s.expand(value)==0
    first=D[:70,:];C=s.Matrix.hstack(*first.nullspace());res=D[70:,:]*C
    rr,piv=res.to_DM().rref();assert 45-C.cols+len(piv)==42
    rr=rr.to_Matrix();free=[i for i in range(C.cols) if i not in piv];N=s.zeros(C.cols,len(free))
    for j,k in enumerate(free):
     N[k,j]=1
     for i,col in enumerate(piv):N[col,j]=-rr[i,k]
    assert (res*N).applyfunc(s.cancel)==s.zeros(res.rows,N.cols)
    enc=lambda A:[[str(a) for a in row] for row in A.tolist()]
    T=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]])
    V=s.Matrix([x,y,z,w]);Qw=s.zeros(4);Qw[1,2]=x*z+3*w*w;Qw[2,1]=-Qw[1,2];Qw[1,3]=-x*w;Qw[3,1]=x*w
    assert T*J[1:,1:]*T.T==J[1:,1:]
    assert (T.inv()*Q.subs(dict(zip(V,T*V)),simultaneous=True)*T.inv().T+Qw).applyfunc(s.expand)==s.zeros(4)
    return dict(beta_kernel=enc(K),q_point=list(map(str,q0)),constant_kernel=enc(C),residual_matrix=enc(res),residual_kernel=enc(N),free=free,rank=42,symmetric_tail_transport_verified=True)

def branch_certificate():
    m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];two=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==2]
    Q=s.zeros(4);Q[1,2]=x*z;Q[2,1]=-x*z;Q[1,3]=-x*w;Q[3,1]=x*w
    M=beta_matrix(m,Q);H=s.Matrix.hstack(*M.nullspace());zz,bu,B=primitive_beta(m);vector=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vector.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();actions=s.Matrix.hstack(*[(left*infinitesimal_action(vector,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()])
    def coords(Q):
     P=X[0]**2*J;P[1:,1:]+=Q;pp=m.normalize(P)['parameters'];return s.Matrix([pp[i] for i in two])
    q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);chart=[i for i in range(45) if i!=fixed]
    qd=[]
    for A in sp_basis():
     dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T
     dq=coords(dQ);dq-=q0*dq[fixed]/q0[fixed];qd.append(dq[chart,:])
    QD=s.zeros(44,H.cols).row_join(s.Matrix.hstack(*qd));data={};out={}
    for label,a in [('ratio3/5',6),('ratio2/5',4)]:
     point=s.zeros(61,1)
     for j,i in enumerate(one):point[i]=j+2
     for i in [3,4,11,12]:point[i]=0
     point[13]=a;point[48]=10;point[49]=-5*point[34];point[50]=5*point[25];point[40]=-1
     data[label]=calculate(m,True,point)
    extra=json.loads(Path('reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation'];top=[[t for t in eq if sum(i in two for i in t[0])==2] for eq in m.reduced['equations']+[extra]]
    for label,G in [('ratio3/5',z**3),('ratio2/5',w**3)]:
     r=data[label];point=s.Matrix(list(map(s.Rational,r['point61'])));sub=dict(zip(bu,[point[i] for i in one]));BD=H.row_join(actions.subs(sub));kernel=s.Matrix.hstack(*BD.nullspace());C=QD*kernel
     _,piv=C.to_DM().rref();assert BD.rank()==16 and len(piv)==1
     constraint=s.zeros(1,16);constraint[0,one.index(13)]=1;constraint[0,one.index(48)]=-point[13]/point[48]
     HH=H*s.Matrix.hstack(*(constraint*H).nullspace());assert HH.cols==9
     BDres=HH.row_join(actions.subs(sub));Qres=s.zeros(44,HH.cols).row_join(QD[:,H.cols:])
     assert BDres.rank()==16 and BDres.col_join(Qres).rank()==16
     orbit=QD[:,H.cols:];stabilizer=s.Matrix.hstack(*orbit.nullspace());assert orbit.rank()==7 and stabilizer.cols==3
     bb=s.symbols('bb0:9');beta_parameters=HH*s.Matrix(bb);stab_action=(actions*stabilizer).subs(dict(zip(bu,beta_parameters)),simultaneous=True)
     assert (constraint*stab_action).applyfunc(s.expand)==s.zeros(1,3)
     assert (M*stab_action).applyfunc(s.expand)==s.zeros(134,3)
     minor_cols=list(BDres.to_DM().rref()[1]);node_minor=BDres[:,minor_cols].det(method='domain-ge');assert node_minor!=0
     main=C[:,piv[0]];Qt=s.zeros(4);Qt[1,2]=s.diff(G,w);Qt[2,1]=-Qt[1,2];Qt[1,3]=-s.diff(G,z);Qt[3,1]=-Qt[1,3];tail=coords(Qt)[chart,:]
     c=s.Symbol('c');family=point.copy()
     qt=coords(Qt)
     for j,i in enumerate(two):family[i]+=c*qt[j]
     assert all(s.expand(m.evaluate(eq,family))==0 for eq in top)
     T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational);jr=s.Matrix(r['jacobian']).applyfunc(s.Rational)
     assert jr*main.col_join(s.zeros(1,1))==s.zeros(135,1);assert jr*tail.col_join(s.zeros(1,1))==s.zeros(135,1)
     rr=list(T.T.rref()[1]);L=T[rr,:].inv();vm=L*main.col_join(s.zeros(1,1))[rr,:];vt=L*tail.col_join(s.zeros(1,1))[rr,:]
     assert s.Matrix.hstack(vm,vt).rank()==2
     a,b=s.symbols('a b');forms=[]
     for row in r['quadratic_obstruction_basis']:
      forms.append(s.expand(sum(s.Rational(c)*[a,b,-a][i]*[a,b,-a][j] for c,(i,j) in zip(row,r['quadratic_monomials']))))
     assert len(forms)==2 and s.gcd(*forms)==1
     gb=s.groebner(forms,a,b);leading=[p.monoms()[0] for p in gb.polys]
     standard=[(i,j) for i in range(5) for j in range(5) if not any(i>=v and j>=w for v,w in leading)]
     assert len(standard)==4 and all(i+j<=2 for i,j in standard)
     assert vm[0]+vm[2]!=0 and vt[0]+vt[2]!=0
     out[label]=dict(beta_rank=16,node_incidence_dimension=16,node_beta_minor_columns=minor_cols,node_beta_minor=str(node_minor),projective_Q_orbit_dimension=7,restricted_beta_kernel=[[str(v) for v in row] for row in HH.tolist()],restricted_beta_differential=[[str(v) for v in row] for row in BDres.tolist()],connected_Q_stabilizer_preserves_restricted_beta=True,main_tangent=list(map(str,vm)),tail_tangent=list(map(str,vt)),whole_tail_leading_zero=True,
      beta_parameter_differential=[[str(v) for v in row] for row in BD.tolist()],local_data=r,
      transverse_parameter='v0+v2',transverse_initial_quadrics=list(map(str,forms)),transverse_initial_groebner=list(map(str,gb)),
      transverse_standard_monomials=standard,transverse_length_upper_bound=4)
    return out

def verify():
    tail=tail_certificate();branches=branch_certificate()
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_NONBINARY_FORMAL_RING_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_NONBINARY_RESONANCE_NODES_20261101',all_passed=True,field='QQ',tail_universal_rank=tail,cases=branches,
        conclusion='At both specified135 homogenized points, t=0; the completed local ring is a one-dimensional complete intersection, finite free rank4 over C[[v0+v2]], with reduced ordinary node and square-zero nilradical.',
        proof_note='research/notes/S61_NONBINARY_RESONANCE_NODES_20261101.md',
        limitations=['The two specified local node rings and generic-beta removal of their sixteen-dimensional node-incidence closures; not all R0 endpoints or all saturated infinity.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,cases=2,tail_leading_rank=42,transverse_length=4)))
if __name__=='__main__':main()
