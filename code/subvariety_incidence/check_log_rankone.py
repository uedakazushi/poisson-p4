"""Exact rank-one weight-block resonant Log boundary checks over Q.
Geometric proof in quadric_log_orthogonality.tex; finite tests are supporting data.
"""
from itertools import combinations
import argparse,json
from pathlib import Path
import sympy as s
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path('work/log_rankone.json'))
args=parser.parse_args();assert s.__version__=='1.14.0'
x=s.symbols('x0:5');xx=s.Matrix(x);u0,u1,z,w0,w1=x;pairs=list(combinations(range(5),2))
xi=s.Matrix([u0,u1,-2*z,0,0]);R=s.Matrix([0,0,0,w1,-w0]);A=s.Matrix([0,-u1,z,0,0]);
J=s.zeros(5);J[3,4]=1;J[4,3]=-1;Ju=s.zeros(5);Ju[1,2]=1;Ju[2,1]=-1
T,V=s.symbols('t v')
def wedge(a,b):return a*b.T-b*a.T
def jacobi(P):
 return [s.expand(sum(P[i,h]*s.diff(P[j,k],x[h])+P[j,h]*s.diff(P[k,i],x[h])+P[k,h]*s.diff(P[i,j],x[h]) for h in range(5)))for i,j,k in combinations(range(5),3)]
def strict_rank(P,q):
 G=s.hessian(q,xx)/2;zz=s.symbols('a0:10');Z=s.zeros(5)
 for (i,j),v in zip(pairs,zz):Z[i,j]=v;Z[j,i]=-v
 B=G.inv()*Z;vel=B*xx;eq=[]
 for i,j in pairs:
  pol=s.expand(sum(vel[h]*s.diff(P[i,j],x[h]) for h in range(5))-
    sum(B[i,h]*P[h,j]+P[i,h]*B[j,h] for h in range(5)))
  eq.extend(s.Poly(pol,*x).coeffs())
 return s.linear_eq_to_matrix(eq,zz)[0].rank()
ans=[]
for c in (s.Integer(0),s.Integer(1)):
 qU=u0*u0+2*u1*z+c*z*z;q=qU+w0*w0+w1*w1
 g=s.expand(sum(xi[i]*s.diff(qU,x[i])for i in range(3)))
 P0=2*wedge(xi,R)+g*J;KR=q*J+wedge(xx,R)
 Pline=P0+T*KR
 assert all(k==0 for k in jacobi(Pline))
 assert all(s.expand(v)==0 for v in Pline*s.Matrix([s.diff(q,a)for a in x]))
 rk=strict_rank(P0,q)
 assert rk == (8 if c==0 else 9)
 entry={'c':str(c),'strict_so5_rank':rk,'strict_stabilizer_dimension':10-rk,'line_jacobi':True}
 if c==0:
  KA= -q*Ju+wedge(xx,A)
  P=Pline+V*KA
  assert all(s.expand(v)==0 for v in P*s.Matrix([s.diff(q,a)for a in x]))
  eq=[]
  for z0 in jacobi(P):eq += s.Poly(z0,*x).coeffs()
  gb=s.groebner(eq,T,V)
  assert len(gb.polys)==1 and s.expand(gb.polys[0].as_expr()-(T-1)*V)==0
  quartic=[]
  for missing in range(5):
   i,j,k,l=[e for e in range(5)if e!=missing]
   quartic.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
  f=-V*u0*(T*q+2*u0*u0-2*u1*z+2*w0*w0+2*w1*w1)/2
  assert all(s.expand(quartic[i]-f*s.diff(q,x[i]))==0 for i in range(5))
  assert s.factor(f.subs(T,1))==-s.Rational(3,2)*V*u0*(u0*u0+w0*w0+w1*w1)
  entry.update({'two_parameter_Jacobi_ideal':'(v*(t-1))','rank4_divisor':'Q + V(x0) + V(x0^2+x3^2+x4^2)','pfaffian_f':str(s.factor(f))})
 ans.append(entry)
res={'status':'PASS','sympy':s.__version__,'scope':'Exact rank-one H cases c=0 and c=1; proof of all allowed rank-one H via coordinate normal forms and strict stabilizer argument is in TeX.','cases':ans}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
