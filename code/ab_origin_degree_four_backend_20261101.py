"""Exact incidence coordinates and Levi/Bott operations for AB normalization."""
from collections import Counter
from itertools import combinations,combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
x,y,z,w,t=s.symbols('x y z w t');Y=s.Matrix([x,y,z,w]);q=s.Matrix([x,z]);J2=s.Matrix([[0,1],[-1,0]])
ss=s.symbols('s0:3');hh=s.symbols('h0:5');pars=[*ss,*hh]
S=ss[0]*x*x+ss[1]*x*z+ss[2]*z*z;H=sum(hh[i]*x**(4-i)*z**i for i in range(5))
def Nof(S,H):
 a=[s.Poly(s.expand(S),x,z).coeff_monomial(m) for m in (x*x,x*z,z*z)]
 h=[s.Poly(s.expand(H),x,z).coeff_monomial(x**(4-i)*z**i) for i in range(5)]
 return s.Matrix([[h[2]/6,-h[1]/4,h[0]],[h[3]/2,-h[2]/3,h[1]/2],[h[4],-h[3]/4,h[2]/6]])+s.Matrix([[2*a[1],-2*a[0],0],[4*a[2],0,-4*a[0]],[0,2*a[2],-2*a[1]]])/96
N=Nof(S,H);E=N.reshape(9,1).jacobian(pars);inverse=(E.T*E).inv()*E.T
F=q*S/4+J2*s.Matrix([s.diff(H,x),s.diff(H,z)]);mixed=F.jacobian(q)
def Qof(N,vertical=0):
 out=s.zeros(4);A=mixed.subs(dict(zip(pars,inverse*N.reshape(9,1))),simultaneous=True)
 for i,ii in enumerate((0,2)):
  for j,jj in enumerate((1,3)):out[ii,jj],out[jj,ii]=A[i,j],-A[i,j]
 out[1,3],out[3,1]=vertical,-vertical
 return out

def upper(B):
 A=s.eye(4)
 for i,ii in enumerate((0,2)):
  for j,jj in enumerate((1,3)):A[ii,jj]=B[i,j]
 return A

def push(A,Q):return (A*Q.subs(dict(zip(Y,A.inv()*Y)),simultaneous=True)*A.T).applyfunc(s.expand)
def extract(Q):
 A=Q[[0,2],[1,3]];F=A*q/3;S=s.diff(F[0],x)+s.diff(F[1],z);R=F-q*S/4
 H=-s.Poly(R[1],x,z).coeff_monomial(x**3)*x**4/4
 for i in range(1,5):H+=s.Poly(R[0],x,z).coeff_monomial(x**(4-i)*z**(i-1))*x**(4-i)*z**i/i
 return Nof(S,H),s.Matrix([s.Poly(Q[1,3],x,z).coeff_monomial(m) for m in (x*x,x*z,z*z)])


core_weights=[(2,0),(1,1),(0,2),(3,-1),(2,0),(1,1),(0,2),(-1,3)]
vertical_weights=[(3,1),(2,2),(1,3)]
rels=[s.expand(N.extract(rr,cc).det()) for rr in combinations(range(3),2) for cc in combinations(range(3),2)]
gb=s.groebner(rels,*pars,domain=s.QQ,order='grevlex');leading=[p.LM(order='grevlex').exponents for p in gb.polys]
def character(n,weights,exclude=()):
    result=Counter()
    for inds in combinations_with_replacement(range(len(weights)),n):
        powers=tuple(inds.count(i) for i in range(len(weights)))
        if any(all(a>=b for a,b in zip(powers,mm)) for mm in exclude):continue
        result[tuple(sum(weights[i][r] for i in inds) for r in range(2))]+=1
    return result

def decompose(character):
    c=Counter({a:n for a,n in character.items() if n});out=[]
    while c:
        a,b=max(c);m=c[(a,b)];assert a>=b and m>0;out.append((a,b,m))
        for i in range(a-b+1):
            key=a-i,b+i;c[key]-=m;assert c[key]>=0
            if c[key]==0:del c[key]
    return out

def bott(a,b):
    u,v=a+2,b+1
    if 0 in (u,v) or abs(u)==abs(v):return None
    aa,bb=sorted([abs(u),abs(v)],reverse=True);aa-=2;bb-=1
    degree=(u<0)+(v<0)+(u+v<0);A,B=aa-bb,bb
    return degree,(aa,bb),(A+1)*(B+1)*(A+B+2)*(A+2*B+3)//6
