#!/usr/bin/env python3
"""Independent exact check of the universal A-invariant quadratic Poisson pencil.
SymPy 1.14.0. No network, previous producer, or untracked input.
"""
import argparse,itertools as it,json
from pathlib import Path
import sympy as s
p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,default=Path('work/universal_nilpotent_pencil.json'))
p.add_argument('--check',type=Path)
a=p.parse_args()
assert s.__version__=='1.14.0'
u0,u1,z,w0,w1=s.symbols('u0 u1 z w0 w1')
alpha,beta,gamma,delta,eta,t,lam,eps=s.symbols('alpha beta gamma delta eta t lam eps')
xs=(u0,u1,z,w0,w1);E=s.Matrix(xs)
q=u0**2+2*u1*w0+2*z*w1+alpha*u1**2+beta*u1*z+gamma*u1*w1+delta*z**2+eta*w1**2
pairs=list(it.combinations(range(5),2))
def sg(a):return (-1)**sum(a[i]>a[j]for i in range(len(a))for j in range(i+1,len(a)))
def P(q,xi):
 g={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]};B=s.zeros(5)
 for i,j in pairs:
  rem=tuple(k for k in range(5)if k not in(i,j))
  B[i,j]=s.expand(sg((i,j,*rem))*sum((-1)**r*s.diff(q,xs[k])*g.get(tuple(rem[h]for h in range(3)if h!=r),0)for r,k in enumerate(rem)))
  B[j,i]=-B[i,j]
 return B
def Jac(B):
 return [s.factor(sum(B[i,l]*s.diff(B[j,k],xs[l])+B[j,l]*s.diff(B[k,i],xs[l])+B[k,l]*s.diff(B[i,j],xs[l])for l in range(5)))for i,j,k in it.combinations(range(5),3)]
A=s.zeros(5);A[0,1]=1;A[3,0]=-1
G=s.hessian(q,E)/2
assert s.factor(G.det())==1-delta*eta
assert A.T*G+G*A==s.zeros(5)
xi=(u0+t*u1,u1,-2*z)
B0=P(q,xi)
K=(q*A*G.inv()/2+(E*(A*E).T-(A*E)*E.T)/2).applyfunc(s.cancel)
assert Jac(B0)==[0]*10 and Jac(K)==[0]*10 and Jac(B0+K)==[0]*10
B=B0+lam*K
for omit in range(5):
 i,j,k,l=[h for h in range(5)if h!=omit]
 pf=(-1)**omit*(B[i,j]*B[k,l]-B[i,k]*B[j,l]+B[i,l]*B[j,k])
 assert s.cancel(pf-s.Rational(3,2)*lam*u1**2*z*s.diff(q,xs[omit]))==0
# Explicit logarithmic residue model over Q(alpha,...,eps), all ten coefficients.
P_eps=P(q,(u0+t*u1,(1+eps)*u1,-(2+eps)*z))
r1=s.Matrix([1,-1,0,0]);r2=s.Matrix([1,0,-1,0]);r3=s.Matrix([2,0,0,-1])
res=[]
for v,w in [(r1,r2),(r1,r3),(r2,r3)]:
 res.append({ij:v[ij[0]]*w[ij[1]]-v[ij[1]]*w[ij[0]]for ij in it.combinations(range(4),2)})
c={ij:2*res[0][ij]-(3+eps)*res[1][ij]-eps*res[2][ij]for ij in res[0]}
f=(u0-t*u1/eps,u1,z,q)
for i,j in pairs:
 actual=0
 for k in range(5):
  if k in(i,j):continue
  h,l=tuple(h for h in range(5)if h not in(k,i,j))
  actual+=xs[k]*sg((h,l,k,i,j))*P_eps[h,l]
 predicted=sum(c[h,l]*s.prod(f[k]for k in range(4)if k not in(h,l))*(s.diff(f[h],xs[i])*s.diff(f[l],xs[j])-s.diff(f[h],xs[j])*s.diff(f[l],xs[i]))for h,l in it.combinations(range(4),2))
 assert s.cancel(actual-predicted)==0,(i,j)
assert P_eps.subs(eps,0)==B0
# True affine-linear change to the already proved long nonconstant-g source.
x,y,Z,w=s.symbols('x y Z w')
V=(x,y,Z,w)
chart=s.Matrix([u0/u1,-w0/u1,z/u1,w1/u1])
J=chart.jacobian(E)
Baf=(J*B*J.T).subs({u0:x,u1:1,z:Z,w0:-y,w1:w},simultaneous=True).applyfunc(s.cancel)
mod=[s.factor(sum(s.diff(Baf[i,j],V[i])for i in range(4)))for j in range(4)]
assert mod==[-lam,-lam*x+10*(Z+eta*w)+5*gamma,0,10]
D=lam**2+100*eta
change=s.Matrix([(-lam*x+10*(Z+eta*w)+5*gamma)/D,y/D,Z,w+10*x/lam])
M=change.jacobian(V)
assert s.factor(M.det())==-1/(lam*D)
Xn,Yn,Zn,Wn=s.symbols('Xn Yn Zn Wn')
oldx=-lam*Xn+lam*(10*(Zn+eta*Wn)+5*gamma)/D
inverse={x:oldx,y:D*Yn,Z:Zn,w:Wn-10*oldx/lam}
T=(M*Baf*M.T).subs(inverse,simultaneous=True).applyfunc(s.cancel)
modn=[s.factor(sum(s.diff(T[i,j],(Xn,Yn,Zn,Wn)[i])for i in range(4)))for j in range(4)]
assert modn==[1,Xn,0,0]
r=Yn-Xn**2/2
qaff=x**2-2*y+2*Z*w+alpha+beta*Z+gamma*w+delta*Z**2+eta*w**2
Ksrc=s.factor(-qaff.subs(inverse,simultaneous=True)/(2*D)-r)
assert s.diff(Ksrc,Xn)==s.diff(Ksrc,Yn)==0
UZ=s.factor(T[0,2]);UW=s.factor(T[0,3]);g=6*Zn
div=s.factor(s.diff(UZ,Zn)+s.diff(UW,Wn))
assert s.cancel(UZ+60*eta*Zn/D)==0
assert s.cancel(UW-2*t/lam-60*Zn/D)==0
assert s.cancel(div+60*eta/D)==0
assert s.cancel(UZ*s.diff(g,Zn)+UW*s.diff(g,Wn)-div*g)==0
F=s.factor(UZ*s.diff(Ksrc,Zn)+UW*s.diff(Ksrc,Wn)+(1+div)*Ksrc)
src={(0,1):-F-(1+div)*r,(0,2):UZ,(0,3):UW,
(1,2):Xn*UZ+r*s.diff(g,Wn)+s.diff(g*Ksrc,Wn),
(1,3):Xn*UW-r*s.diff(g,Zn)-s.diff(g*Ksrc,Zn),(2,3):g}
for ij,val in src.items():assert s.cancel(T[ij]-val)==0,ij
# Quadratic normal form: kappa solves eta*kappa**2-2*kappa+delta=0.
kap=s.symbols('kap')
eq=eta*kap**2-2*kap+delta
assert s.expand((1-eta*kap)**2-(1-delta*eta)-eta*eq)==0
result={
 'date':'2026-11-10','version':'v15','status':'PASS','sympy':'1.14.0',
 'A':'u1*d/du0-u0*d/dw0','quadratic_parameters':5,
 'smooth_condition':'1-delta*eta != 0',
 'full_Jacobi':'all ten strict components vanish over Q(alpha,beta,gamma,delta,eta,t,lam)',
 'Pfaffian':'(3*lam/2)*u1**2*z*dq; all five components checked',
 'log_source':'(u0-t*u1/eps,u1,z,q), residues (2,-3-eps,-eps)',
 'log_twoform_coefficients_checked':10,
 'long_affine_change_det':'-1/(lam*(lam**2+100*eta))',
 'long_source_g':'6*Zn','long_source_UZ':'-60*eta*Zn/(lam**2+100*eta)',
 'long_source_UW':'2*t/lam+60*Zn/(lam**2+100*eta)',
 'long_source_six_brackets_checked':6,
 'quadratic_normal_form':'eta_new=eta/(1-delta*eta); eta*kappa**2-2*kappa+delta=0',
 'scope':'Every smooth q in the five-parameter invariant slice has rank-two Log and rank-four long-component specializations. No decision between global intersection dimensions 23 and 24.'}
if a.check is not None:
 assert result==json.loads(a.check.read_text(encoding='utf-8')), 'certificate mismatch'
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('PASS: universal invariant quadratic pencil, 10 log and 6 long-source coefficients')
