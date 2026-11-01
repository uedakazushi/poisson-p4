"""Finite stabilizer slice and actual two-dimensional fixed-jet T5 incidence."""
import sympy as s
from s61_model import X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action


def verify(model,orbit,frame):
    b=s.symbols('b0:12');H=s.Matrix(orbit['beta_kernel']).applyfunc(s.Rational)
    z,u,beta=primitive_beta(model);v=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS])
    basis=v.jacobian(u);rows=list(basis.T.rref()[1]);inverse=basis[rows,:].inv()
    hr=list(H.T.rref()[1]);hi=H[hr,:].inv();bet=beta.xreplace(dict(zip(u,H*s.Matrix(b))))
    x,y,zz,w=X[1:];coordinates=s.Matrix([x,y,zz,w]);Q=s.zeros(4);Q[1,2]=2*x*w;Q[1,3]=-x*x;Q-=Q.T
    A0=s.zeros(4);A0[1,0]=-1
    A1=s.zeros(4);A1[1,3]=-1;A1[2,0]=1
    A3=s.zeros(4);A3[1,2]=-s.Rational(1,2);A3[2,3]=1;A3[3,0]=-s.Rational(1,2)
    parameter=s.Symbol('s');lam=s.Symbol('lam',nonzero=True);L=68*b[0]-23*b[10]-121*b[5]
    finite={}
    def finite_beta(g):
        new=g.T*bet.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g
        flattened=s.Matrix([s.expand(new[i,j]).diff(z[k]) for (i,j),k in COORDS])
        un=inverse*flattened[rows,:];bn=(hi*un[hr,:]).applyfunc(s.expand)
        assert (H*bn-un).applyfunc(s.expand)==s.zeros(16,1)
        return bn
    for label,A in [('0',A0),('1',A1),('3',A3)]:
        assert A**4==s.zeros(4)
        g=s.eye(4)+parameter*A+parameter**2*A**2/2+parameter**3*A**3/6
        assert (g.T*J[1:,1:]*g-J[1:,1:]).applyfunc(s.expand)==s.zeros(4)
        assert (g.inv()*Q.subs(dict(zip(coordinates,g*coordinates)),simultaneous=True)*g.inv().T-Q).applyfunc(s.expand)==s.zeros(4)
        bn=finite_beta(g);assert bn[9]==b[9] and s.expand(68*bn[0]-23*bn[10]-121*bn[5]-L)==0
        finite[label]=dict(matrix=[[str(x) for x in row] for row in g.tolist()],beta_action=list(map(str,bn)))
    g=s.diag(lam**-3,lam**3,lam,lam**-1);assert g.T*J[1:,1:]*g==J[1:,1:]
    assert (g.inv()*Q.subs(dict(zip(coordinates,g*coordinates)),simultaneous=True)*g.inv().T-lam**-8*Q).applyfunc(s.expand)==s.zeros(4)
    bn=finite_beta(g);assert bn[9]==lam**3*b[9] and s.expand(68*bn[0]-23*bn[10]-121*bn[5]-lam*L)==0
    finite['torus']=dict(matrix=[[str(x) for x in row] for row in g.tolist()],beta_action=list(map(str,bn)),q_scale='lam**(-8)')
    assert s.sympify(finite['3']['beta_action'][5])==b[5]+7*parameter*b[9]/20
    assert s.sympify(finite['1']['beta_action'][2])==b[2]-4*parameter*b[9]/5
    assert s.sympify(finite['0']['beta_action'][7])==b[7]+2*parameter*b[9]/5
    # The four constraints are attained successively, with denominators only b9,L.
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    point=s.Matrix(orbit['point61']).applyfunc(s.Rational);original=hi*point[one,:][hr,:]
    assert H*original==point[one,:]
    values=list(original);lv=L.subs(dict(zip(b,values)))
    values=[s.cancel(s.sympify(x).subs(dict(zip(b,values)),simultaneous=True).subs(s.Symbol('lam'),1/lv)) for x in finite['torus']['beta_action']]
    for label,scalar in [('3',lambda v:20*(1-v[5])/(7*v[9])),('1',lambda v:5*v[2]/(4*v[9])),('0',lambda v:-5*v[7]/(2*v[9]))]:
        ss=scalar(values)
        values=[s.cancel(s.sympify(x).subs(dict(zip(b,values)),simultaneous=True).subs(parameter,ss)) for x in finite[label]['beta_action']]
    assert values[2]==values[7]==0 and values[5]==1 and L.subs(dict(zip(b,values)))==1
    # Recover the full normal beta differential and the SAME actual surface tangents.
    action=s.Matrix.hstack(*[inverse*infinitesimal_action(v,A)[rows,:] for A in sp_basis()])
    left=s.Matrix.hstack(*H.T.nullspace()).T
    N=(left*action.xreplace(dict(zip(u,H*s.Matrix(b))))).applyfunc(s.expand)
    nonzero=[4,5,6,7];free=[0,1,2,3,8,9]
    assert N[:,free]==s.zeros(4,6)
    P=1425*b[0]*b[10]*b[5]+1375*b[0]*b[5]**2-1050*b[0]*b[6]*b[9]+105*b[10]*b[2]*b[9]-180*b[11]*b[5]*b[9]+275*b[2]*b[5]*b[9]+210*b[7]*b[9]**2+84*b[8]*b[9]**2
    determinant=s.factor(N[:,nonzero].det());assert s.expand(determinant+2*b[9]*P/375)==0
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    def qcoords(Q):
        A=X[0]**2*J;A[1:,1:]+=Q
        return s.Matrix(model.normalize(A)['parameters'])[two,:]
    q=qcoords(Q);fixed=frame['fixed_q'];chart=frame['chart_columns'];OD=[]
    assert q==point[two,:]
    for A in sp_basis():
        dQ=sum(((A*coordinates)[i]*Q.diff(coordinates[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T
        dq=qcoords(dQ);dq-=q*dq[fixed]/q[fixed];OD.append(dq[chart,:])
    OD=s.Matrix.hstack(*OD);assert OD.rank()==6
    surface=OD[frame['free_columns'],free]
    assert surface[2:,:]==s.zeros(3,6) and surface[:2,:].rank()==2
    assert b[9].subs(dict(zip(b,original)))!=0 and lv!=0 and P.subs(dict(zip(b,original)))!=0
    assert b[9].subs(dict(zip(b,values)))!=0 and P.subs(dict(zip(b,values)))!=0
    return dict(finite_actions=finite,normalization_domain='b9*(68*b0-23*b10-121*b5)!=0',
                normalization_steps=['lam=1/L','s3=20*(1-b5)/(7*b9)','s1=5*b2/(4*b9)','s0=-5*b7/(2*b9)'],
                original_point_b=list(map(str,original)),normalized_point_b=list(map(str,values)),
                beta_normal_action=[[str(x) for x in row] for row in N.tolist()],beta_normal_minor=str(determinant),
                actual_surface_tangent=[[str(x) for x in row] for row in surface.tolist()],
                source_beta_rank=16,source_joint_rank=18,surface_parameter_rows=[0,1],
                relative_actual_surface_on='b9*P!=0',old_point_in_all_required_opens=True)
