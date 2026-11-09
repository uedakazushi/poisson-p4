"""Exact generic logarithmic Poisson extensions on a smooth quadric.

A universal quadratic form on U=span(x0,x1,x2) and a 2-dimensional orthogonal
complement are used for the all-parameter Jacobi/Pfaffian identities. The
strict so(q) stabilizer rank is computed at one exact rational sample and
promoted to a general-position statement by openness and the universal SO(2)
symmetry. Requires SymPy 1.14.0. No external or historical computation inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,default=Path('work/log_quadric.json'))
args=p.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('x0:5')
z=s.symbols('z0:10')
G=s.symbols('g00 g01 g02 g11 g12 g22')
a,b,t=s.symbols('a b t')
qU=G[0]*x[0]**2+2*G[1]*x[0]*x[1]+2*G[2]*x[0]*x[2]+G[3]*x[1]**2+2*G[4]*x[1]*x[2]+G[5]*x[2]**2
q=qU+x[3]**2+x[4]**2
R=s.Matrix([0,0,0,x[4],-x[3]])
E=s.Matrix([*x[:3],0,0])
xi=s.Matrix([a*x[0],b*x[1],-(a+b)*x[2],0,0])
J=s.zeros(5);J[3,4]=1;J[4,3]=-1

def wed(A,B):return A*B.T-B*A.T

def tensor(A):
    g=s.expand(sum(A[i]*s.diff(qU,x[i]) for i in range(3))/2)
    return wed(A,R)+g*J

P0=tensor(2*xi)
K=tensor(E)
assert K==q*J+wed(s.Matrix(x),R)
P=P0+t*K

# In fact the equality Jacobi=0, Pfaffian=0, and q-Casimir holds
# over Q[a,b,g00,g01,...,t,x0,...,x4], not just at samples.
qgrad=s.Matrix([s.diff(q,u) for u in x])
assert all(s.expand(z)==0 for z in P*qgrad)
jac=[]
for i,j,k in itertools.combinations(range(5),3):
    v=sum(P[i,l]*s.diff(P[j,k],x[l])+P[j,l]*s.diff(P[k,i],x[l])+P[k,l]*s.diff(P[i,j],x[l]) for l in range(5))
    jac.append(s.expand(v))
assert jac==[0]*10
pf=[]
for i,j,k,l in itertools.combinations(range(5),4):
    pf.append(s.expand(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
assert pf==[0]*5
assert all(s.expand(z)==0 for z in K.subs({G[0]:1,G[1]:0,G[2]:0,G[3]:1,G[4]:0,G[5]:1})*qgrad.subs({G[0]:1,G[1]:0,G[2]:0,G[3]:1,G[4]:0,G[5]:1}))
assert all(s.expand(P0[i,j].diff(t))==0 for i,j in itertools.combinations(range(5),2))
assert all(s.expand(z)==0 for z in [P0[i,j] - P0[i,j].subs(t,0) for i,j in itertools.combinations(range(5),2)])
div=s.Matrix([sum(s.diff(P[i,j],x[i]) for i in range(5))for j in range(5)])
assert all(s.expand(div[i]-3*t*R[i])==0 for i in range(5))

# A generic log form comes from three independent linear hyperplanes with
# weights of sum zero. The stabilizer on Q equals the strict stabilizer of
# the divergence-free lift P0; for the incidence this uses the kernel K_T.
# Choose a rational orthogonal sample, a= -5/3,b=4/3,c=1/3.
point={G[0]:1,G[1]:0,G[2]:0,G[3]:1,G[4]:0,G[5]:1,a:s.Rational(-5,3),b:s.Rational(4,3)}
Psample=P0.subs(point)
pairs=list(itertools.combinations(range(5),2));A=s.zeros(5)
for (i,j),v in zip(pairs,z):A[i,j]=v;A[j,i]=-v
V=A*s.Matrix(x)
eq=[]
for i,j in pairs:
    L=s.expand(sum(V[k]*s.diff(Psample[i,j],x[k]) for k in range(5))-
               sum(A[i,k]*Psample[k,j]+Psample[i,k]*A[j,k] for k in range(5)))
    eq.extend(s.Poly(L,*x).coeffs())
M=s.linear_eq_to_matrix(eq,z)[0]
assert M.rank()==9
row=list(M.T.rref()[1]);col=list(M.rref()[1])
assert row==[0,1,2,3,4,5,11,13,33]
assert col==list(range(9))
minor=s.factor(M.extract(row,col).det())
assert minor==s.Rational(10240000,243)
null=M.nullspace()
assert len(null)==1 and null[0]==s.eye(10)[:,9]

# A special boundary of the Log source: two weights collide and the third
# residue becomes zero. Its intrinsic Poisson section has a hyperplane factor.
# The strict stabilizer is SO(2) on (x0,x1) times SO(2) on (x3,x4).
special_point={G[0]:1,G[1]:0,G[2]:0,G[3]:1,G[4]:0,G[5]:1,a:1,b:1}
Ps=P0.subs(special_point)
eqs=[]
for i,j in pairs:
    L=s.expand(sum(V[k]*s.diff(Ps[i,j],x[k]) for k in range(5))-
               sum(A[i,k]*Ps[k,j]+Ps[i,k]*A[j,k] for k in range(5)))
    eqs.extend(s.Poly(L,*x).coeffs())
Ms=s.linear_eq_to_matrix(eqs,z)[0]
assert Ms.rank()==8
rs=list(Ms.T.rref()[1]);cs=list(Ms.rref()[1]);msdet=Ms.extract(rs,cs).det()
assert rs==[0,1,2,3,4,5,12,30] and cs==[1,2,3,4,5,6,7,8] and msdet==2304
assert sorted([list(w) for w in Ms.nullspace()]) == sorted([list(s.eye(10)[:,0]),list(s.eye(10)[:,9])])
S=s.Matrix([x[1],-x[0],0,0,0])
JS=s.zeros(5);JS[0,1]=1;JS[1,0]=-1
Es=s.Matrix(x)
# The symplectic normal vector to the constant bivector JS is -2 S.
KS=q.subs(special_point)*JS+wed(Es,S)
KR=K.subs(special_point)
vpar=s.Symbol('v')
Pboundary=Ps+t*KR+vpar*KS
assert all(s.expand(z)==0 for z in Pboundary*qgrad.subs(special_point))
coeff=[]
for i,j,k in itertools.combinations(range(5),3):
    pol=s.Poly(s.expand(sum(Pboundary[i,l]*s.diff(Pboundary[j,k],x[l])+
      Pboundary[j,l]*s.diff(Pboundary[k,i],x[l])+Pboundary[k,l]*s.diff(Pboundary[i,j],x[l])
      for l in range(5))),*x)
    coeff.extend(pol.coeffs())
conditions=s.groebner(coeff,t,vpar)
assert len(conditions.polys)==1
assert s.expand(conditions.polys[0].as_expr()-vpar*(t+2))==0
PF=[]
for missing in range(5):
    i,j,k,l=[n for n in range(5) if n!=missing]
    PF.append((-1)**missing*(Pboundary[i,j]*Pboundary[k,l]-Pboundary[i,k]*Pboundary[j,l]+Pboundary[i,l]*Pboundary[j,k]))
f_boundary=s.factor(s.cancel(PF[0]/(2*x[0])))
assert s.expand(f_boundary-vpar*x[2]*(t*q.subs(special_point)+2*x[0]**2+2*x[1]**2-4*x[2]**2-4*x[3]**2-4*x[4]**2)/2)==0
assert all(s.expand(PF[i]-2*x[i]*f_boundary)==0 for i in range(5))
f_second=s.factor(f_boundary.subs(t,-2))
assert s.expand(f_second+3*vpar*x[2]*(x[2]**2+x[3]**2+x[4]**2))==0
assert s.expand(f_boundary.subs(vpar,0))==0
D=[s.expand(sum(s.diff(Pboundary[i,j],x[i]) for i in range(5))) for j in range(5)]
assert all(s.expand(D[j]-3*t*R[j]-3*vpar*S[j])==0 for j in range(5))

result={
 'date':'2026-11-10','status':'PASS','sympy':s.__version__,
 'scope':'universal rank-two Poisson pencil identity over Q[a,b,G,t]; generic strict stabilizer via one nonzero exact minor and upper-semicontinuity',
 'general_log_form':'alpha=sum(lambda_i * prod(l_j,j != i) dl_i), sum(lambda_i)=0',
 'quadratic_form_U':'g00*x0^2+2g01*x0*x1+2g02*x0*x2+g11*x1^2+2g12*x1*x2+g22*x2^2',
 'pencil':'P_t=(2xi+t E_U) wedge R +(xi(q_U)+t q_U) J',
 'parameter_weights':'xi=(a*x0,b*x1,-(a+b)*x2), R=x4*d3-x3*d4, J=d3 wedge d4',
 'universal_jacobi_equations_checked':len(jac),
 'universal_pfaffian_minors_checked':len(pf),
 'universal_q_Casimir_components_checked':5,
 'universal_modular_field':'3*t*R',
 'strict_so5_sample':{'q':'x0^2+x1^2+x2^2+x3^2+x4^2','xi':'(-5*x0/3,4*x1/3,x2/3)',
                       'matrix_shape':list(M.shape),'rank':M.rank(),'pivot_rows':row,'pivot_columns':col,'nonzero_minor':str(minor),
                       'kernel_generator':'rotation in (x3,x4) plane'},
 'mathematical_scope':'generic Log(1,1,1) has reduced A1 of rank-two lifts; one special logarithmic boundary point has a nodal two-line fibre with rank-four second branch. Other special strata not classified.',
 'special_log_boundary':{'weights':'lambda=(3,-3,0), up to common scale','q':'sum_{i=0}^4 xi^2','xi':'x0*d0+x1*d1-2*x2*d2',
 'strict_stabilizer_dimension':2,'strict_stabilizer_matrix_shape':list(Ms.shape),'strict_stabilizer_rank8_witness':{'rows':rs,'columns':cs,'determinant':str(msdet)},'strict_stabilizer_generators':['rotation (x0,x1)','rotation (x3,x4)'],
 'full_extension_coordinates':'P0 + t*KJ + v*KS','Jacobi_scheme_ideal_generator':'v*(t+2)',
 'scheme':'Spec C[t,v]/(v*(t+2)), two reduced affine lines meeting transversely at (-2,0)',
 'pfaffian_cubic':'v*x2*(t*q+2*x0^2+2*x1^2-4*x2^2-4*x3^2-4*x4^2)/2',
 'rank_four_branch_at_t_minus2':'f=-3*v*x2*(x2^2+x3^2+x4^2), v != 0',
 'rank_four_actual_divisor':'q*x2*(x2^2+x3^2+x4^2), up to nonzero scalar',
 'rank_two_branch':'v=0 for all t'}
}
result['producer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'status':result['status'],'strict_rank':M.rank(),'minor':str(minor),'Jacobi':len(jac),'Pfaffians':len(pf)},ensure_ascii=False))
