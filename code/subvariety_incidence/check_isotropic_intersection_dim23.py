#!/usr/bin/env python3
"""Exact common rank-23 family of Log1112 and long rank-four components.
Independent SymPy-1.14 rational-identity producer; no historical imports.
"""
import argparse,itertools,json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ap=argparse.ArgumentParser()
ap.add_argument('--output',type=Path,default=Path('work/isotropic_intersection_dim23.json'))
args=ap.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('u0 u1 z w0 w1');u0,u1,z,w0,w1=x;E=s.Matrix(x)
pairs=list(itertools.combinations(range(5),2))
triples=list(itertools.combinations(range(5),3))
mon2=list(itertools.combinations_with_replacement(range(5),2))
v,t,b,lam,eps=s.symbols('v t b lam eps')
def sg(a):return (-1)**sum(a[i]>a[j]for i in range(len(a))for j in range(i+1,len(a)))
def biv(q,xi):
  gam={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]}
  B=s.zeros(5)
  for i,j in pairs:
    rest=tuple(k for k in range(5)if k not in(i,j))
    B[i,j]=s.expand(sg((i,j,*rest))*sum((-1)**k*s.diff(q,x[l])*gam.get(tuple(rest[a]for a in range(3)if a!=k),0)for k,l in enumerate(rest)))
    B[j,i]=-B[i,j]
  return B
def jac(B):
  return {I:s.expand(sum(B[i,l]*s.diff(B[j,k],x[l])+B[j,l]*s.diff(B[k,i],x[l])+B[k,l]*s.diff(B[i,j],x[l])for l in range(5)))for I in triples for i,j,k in[I]}
def kernel(A,q):
  G=s.hessian(q,E)/2
  assert A.T*G+G*A==s.zeros(5)
  return (q*A*G.inv()/2+(E*(A*E).T-(A*E)*E.T)/2).applyfunc(s.expand)
def vec(B):return s.Matrix([s.Poly(B[i,j],*x).coeff_monomial(x[a]*x[b])for i,j in pairs for a,b in mon2])
def lie(M,B):
  w=M*E
  return s.Matrix([[s.expand(sum(w[l]*s.diff(B[i,j],x[l])-M[i,l]*B[l,j]-M[j,l]*B[i,l]for l in range(5)))for j in range(5)]for i in range(5)])
q=u0**2+2*u1*w0+2*z*w1+v*w1**2
xi=(u0+t*u1,u1,-2*z)
A1=s.zeros(5);A1[0,1]=1;A1[3,0]=-1
A2=s.zeros(5);A2[3,2]=1;A2[3,4]=v;A2[4,1]=-1
B=(biv(q,xi)+b*kernel(A2,q)+lam*kernel(A1,q)).applyfunc(s.expand)
assert not any(jac(B).values())
# Affine chart h=u1, coordinates x=u0/h, y=-w0/h, z=z/h, w=w1/h.
xx,yy,zz,ww=s.symbols('x y zz w')
coords=(xx,yy,zz,ww)
fs=s.Matrix([u0/u1,-w0/u1,z/u1,w1/u1]);Jf=fs.jacobian(E)
P=(Jf*B*Jf.T).subs({u0:xx,u1:1,z:zz,w0:-yy,w1:ww},simultaneous=True).applyfunc(s.expand)
pf=s.factor(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
assert s.expand(pf+3*lam*zz*(xx**2-2*yy+2*zz*ww+v*ww**2))==0
mod=s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i])for i in range(4)))for j in range(4)])
C=10+b;D=lam**2+C**2*v
assert all(s.cancel(z_)==0 for z_ in mod-s.Matrix([-lam,-lam*xx+C*(zz+v*ww),0,C]))
# A GENUINE AFFINE-LINEAR change, invertible for lambda*D != 0.
rect=s.Matrix([(-lam*xx+C*(zz+v*ww))/D,yy/D,zz,ww+C*xx/lam])
Jr=rect.jacobian(coords)
assert s.factor(Jr.det()+1/(lam*D))==0
assert all(s.cancel(z_)==0 for z_ in Jr*mod-s.Matrix([1,rect[0],0,0]))
assert s.cancel((Jr*P*Jr.T)[2,3]-6*zz)==0
X,Y,Z,W=s.symbols('X Y Z W')
ix=-lam*X+lam*C*(Z+v*W)/D
sub={xx:ix,yy:D*Y,zz:Z,ww:W-C*ix/lam}
Pb=(Jr*P*Jr.T).subs(sub,simultaneous=True).applyfunc(s.cancel)
r=Y-X**2/2
qaff=xx**2-2*yy+2*zz*ww+v*ww**2
K=s.cancel(-qaff.subs(sub,simultaneous=True)/(2*D)-r)
assert not any(s.diff(K,z_)for z_ in (X,Y))
uz=s.cancel(Pb[0,2]);uw=s.cancel(Pb[0,3]);g=6*Z
d=s.cancel(s.diff(uz,Z)+s.diff(uw,W))
assert all(s.diff(f,h)==0 for f in (uz,uw)for h in (X,Y))
assert all(s.Poly(f,Z,W).total_degree()<=1 for f in (uz,uw))
assert s.cancel(uz*s.diff(g,Z)+uw*s.diff(g,W)-d*g)==0
F=s.cancel(uz*s.diff(K,Z)+uw*s.diff(K,W)+(d+1)*K)
forms={(0,1):-F-(1+d)*r,(0,2):uz,(0,3):uw,
(1,2):X*uz+r*s.diff(g,W)+s.diff(g*K,W),
(1,3):X*uw-r*s.diff(g,Z)-s.diff(g*K,Z),(2,3):g}
for (i,j),f in forms.items():assert s.cancel(Pb[i,j]-f)==0,(i,j)
assert s.cancel(Pb[0,1]*Pb[2,3]-Pb[0,2]*Pb[1,3]+Pb[0,3]*Pb[1,2]+g*(r+K))==0
# Exact logarithmic source for epsilon !=0, with its well-defined epsilon=0 limit.
xis=(u0+t*u1,(1+eps)*u1,(-2-eps)*z)
Psem=(biv(q,xis)+b*kernel(A2,q)).applyfunc(s.expand)
hypers=(u0-t*u1/eps,u1,z,q)
rr=[s.Matrix([1,-1,0,0]),s.Matrix([1,0,-1,0]),s.Matrix([2,0,0,-1])]
bs=[]
for a1,a2 in itertools.combinations(rr,2):
  bs.append({(i,j):a1[i]*a2[j]-a1[j]*a2[i]for i,j in itertools.combinations(range(4),2)})
weights=[2+b/2,-3-eps,-eps]
resid={(i,j):sum(weights[k]*bs[k][(i,j)]for k in range(3))for i,j in itertools.combinations(range(4),2)}
def omega(B):
  out={}
  for i,j in pairs:
    z_=0
    for k in range(5):
      if k in(i,j):continue
      others=tuple(h for h in range(5)if h not in(k,i,j))
      z_+=x[k]*sg((others[0],others[1],k,i,j))*B[others[0],others[1]]
    out[i,j]=s.expand(z_)
  return out
om=omega(Psem)
for i,j in pairs:
  predicted=0
  for k,l in itertools.combinations(range(4),2):
    fac=s.prod(hypers[h]for h in range(4)if h not in(k,l))
    predicted+=resid[k,l]*fac*(s.diff(hypers[k],x[i])*s.diff(hypers[l],x[j])-s.diff(hypers[k],x[j])*s.diff(hypers[l],x[i]))
  assert s.cancel(om[i,j]-predicted)==0,(i,j)
assert (Psem-B.subs(lam,0)).subs(eps,0).applyfunc(s.expand)==s.zeros(5)
# Tangent rank of PGL5 saturation of the three-parameter log subfamily.
q0=q.subs(v,0);P0=biv(q0,(u0,u1,-2*z))
mats=[]
for i,j in itertools.product(range(5),repeat=2):
  M=s.zeros(5);M[i,j]=1;mats.append(M)
cols=[vec(lie(M,P0))for M in mats]
cols += [vec(biv(w1**2,(u0,u1,-2*z))),vec(biv(q0,(u1,0,0))),vec(kernel(A2.subs(v,0),q0))]
J=s.Matrix.hstack(*cols)
cpiv=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,19,20,22,23,25,26,27)
rpiv=(1,2,5,6,16,17,30,31,32,34,36,39,45,46,47,84,111,114,135,136,137,140,141,144)
minor=DomainMatrix.from_Matrix(J.extract(list(rpiv),list(cpiv))).convert_to(s.QQ).det()
assert minor==10871635968
assert DomainMatrix.from_Matrix(J).convert_to(s.QQ).rank()==24
assert DomainMatrix.from_Matrix(J[:,:25].row_join(vec(P0))).convert_to(s.QQ).rank()==21
gauge=s.Matrix.hstack(*[vec(E*(M*E).T-(M*E)*E.T)for M in mats])
assert DomainMatrix.from_Matrix(gauge).convert_to(s.QQ).rank()==24
assert DomainMatrix.from_Matrix(gauge.row_join(J)).convert_to(s.QQ).rank()==48
out={'date':'2026-11-10','version':'v14','status':'PASS','sympy':s.__version__,
'family':'P(v,t,b,lambda)=P(q_v,xi_t)+b*K_A2(v)(q_v)+lambda*K_A1(q_v)',
'q_v':'u0**2+2*u1*w0+2*z*w1+v*w1**2','xi_t':['u0+t*u1','u1','-2*z'],
'A1':'u1*d/du0-u0*d/dw0','A2_v':'(z+v*w1)*d/dw0-u1*d/dw1',
'full_strict_Jacobi':'identically zero over QQ(v,t,b,lambda)',
'affine_Pfaffian':'-3*lambda*z*(x**2-2*y+2*z*w+v*w**2)',
'affine_modular':['-lambda','-lambda*x+(10+b)*(z+v*w)','0','10+b'],
'affine_linear_rectification':['X=(-lambda*x+(10+b)*(z+v*w))/D','Y=y/D','Z=z','W=w+(10+b)*x/lambda','D=lambda**2+(10+b)**2*v'],
'rectification_Jacobian':'-1/(lambda*D)','normalized_modular':['1','X','0','0'],
'normalized_g':'6*Z (nonconstant)','source_bracket_match':'all six exact polynomial identities over QQ(v,t,b,lambda)',
'rank2_log_semisimple_residue_basis':['(1,-1,0,0)^(1,0,-1,0)','(1,-1,0,0)^(2,0,0,-1)','(1,0,-1,0)^(2,0,0,-1)'],
'rank2_log_semisimple_residue_coefficients':['2+b/2','-3-eps','-eps'],
'rank2_log_F_source':['u0-(t/eps)*u1','u1','z','q_v'],
'rank2_log_identity':'all ten 2-form coefficients agree identically over QQ(v,t,b,eps)',
'linearized_affine_GL5_and_three_parameters_rank':24,
'projective_PGL5_and_three_parameters_rank':23,'gauge_rank':24,'gauge_plus_family_rank':48,
'24x24_minor_row_indices':list(rpiv),'24x24_minor_column_indices':list(cpiv),
'24x24_minor_determinant':str(minor),'intersection_dimension_bounds':[23,24],
'scope':'Uses prior long-source classification and both component theorems. Does not identify the exact dimension (23 or 24), the scheme intersection, or ambient reducedness.'}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print('PASS: universal Poisson family, logarithmic source, long-source six brackets and rank-23 saturation')
