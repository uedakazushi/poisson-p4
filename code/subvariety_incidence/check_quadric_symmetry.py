"""Exact modular/Pfaffian and generic rational stabilizer checks; SymPy 1.14."""
from itertools import combinations
import json
import sys
import sympy as s
x=s.symbols('x0:5');E=s.Matrix(x);pairs=list(combinations(range(5),2))
def wedge(a,b):return a*b.T-b*a.T
def sign(t):return (-1)**sum(t[i]>t[j] for i in range(len(t)) for j in range(i+1,len(t)))
def jac3(q,a,b):
 p=s.zeros(5)
 for i,j in pairs:
  r=[k for k in range(5) if k not in (i,j)]
  p[i,j]=sign((i,j,*r))*s.det(s.Matrix([[s.diff(f,x[k]) for k in r] for f in (q,a,b)]));p[j,i]=-p[i,j]
 return p
def div(p):return s.Matrix([sum(s.diff(p[i,j],x[i]) for i in range(5))for j in range(5)])
def checks(p,q,f):
 assert all(s.expand(z)==0 for z in p*s.Matrix([s.diff(q,z)for z in x]))
 assert all(s.expand(sum(p[i,k]*s.diff(p[j,l],x[k])+p[j,k]*s.diff(p[l,i],x[k])+p[l,k]*s.diff(p[i,j],x[k]) for k in range(5)))==0 for i,j,l in combinations(range(5),3))
 beta=[]
 for miss in range(5):
  i,j,k,l=[a for a in range(5)if a!=miss]
  beta.append((-1)**miss*(p[i,j]*p[k,l]-p[i,k]*p[j,l]+p[i,l]*p[j,k]))
 assert all(s.expand(beta[i]-f*s.diff(q,x[i]))==0 for i in range(5))
 assert s.expand(sum(x[i]*beta[i] for i in range(5))-2*q*f)==0
 return list(map(str,div(p)))
q=sum(a*a for a in x);l=x[0]
r=3*x[0]**2+2*x[1]**2+5*x[2]**2+7*x[3]**2+11*x[4]**2+13*x[1]*x[2]+17*x[2]*x[3]+19*x[3]*x[4]+23*x[1]*x[4]
p=jac3(q,l,r);assert checks(p,q,0)==['0']*5
z=s.symbols('z0:10');A=s.zeros(5)
for (i,j),a in zip(pairs,z):A[i,j]=a;A[j,i]=-a
v=A*E;eq=[]
for i,j in pairs:
 t=s.expand(sum(v[k]*s.diff(p[i,j],x[k]) for k in range(5))-sum(A[i,k]*p[k,j]+p[i,k]*A[j,k] for k in range(5)))
 eq.extend(s.Poly(t,*x).coeffs())
M=s.linear_eq_to_matrix(eq,z)[0];assert M.rank()==10
rr=list(M.T.rref()[1]);cc=list(M.rref()[1]);minor=M.extract(rr,cc).det();assert minor==27980515917168640
u,a,b,c,w=x;qs=2*u*w+a*a+b*b+c*c;U1=s.Matrix([0,u,0,0,-a]);U2=s.Matrix([0,0,u,0,-b]);U3=s.Matrix([0,0,0,u,-c]);T=s.zeros(5);T[4,3]=s.Rational(1,2);T[3,4]=-s.Rational(1,2);t=s.Symbol('t')
ps=wedge(U1,U2)+t*(qs*T-wedge(E,U3)/2)
ds=checks(ps,qs,-t*u**3/4);assert ds==list(map(str,-s.Rational(3,2)*t*U3))
qa=u*w-4*a*c+3*b*b;H=s.Matrix([2*u,a,0,-c,-2*w]);W=s.Matrix([0,u,2*a,3*b,4*c]);T0=s.zeros(5);T0[1,4]=-1;T0[4,1]=1;T0[2,3]=s.Rational(1,2);T0[3,2]=-s.Rational(1,2)
pa=wedge(H,W)-(qa*T0-wedge(E,W)/2)
da=checks(pa,qa,-s.Rational(5,4)*u*qa)
res={'status':'PASS','sympy':s.__version__,'rational_strict_so5_stabilizer_rank':M.rank(),'rational_nonzero_minor':str(minor),'rational_rows':rr,'rational_columns':cc,'simple_quadric_f':'-t*x0**3/4','simple_quadric_divergence':ds,'double_quadric_f':'-5*x0*q/4','double_quadric_divergence':da,'scope':'exact 3 examples; universal modular multiplicity theorem proved separately in manuscript'}
assert s.__version__=='1.14.0'
print(json.dumps(res,indent=2))
if len(sys.argv)>1:
 from pathlib import Path
 path=Path(sys.argv[1]);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(res,indent=2)+'\n')
