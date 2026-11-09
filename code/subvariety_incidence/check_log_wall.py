"""Exact wall-crossing certificate for a resonant Log boundary on Q^3.

Uses only tracked source and SymPy 1.14; no prior producer, hidden data or CAS.
Run: python code/subvariety_incidence/check_log_wall.py --output work/log_wall.json
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import sympy as s

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=Path('work/log_wall.json'))
args = parser.parse_args()
assert s.__version__ == '1.14.0'

x = s.symbols('x0:5'); x0,x1,x2,x3,x4=x
c,t,v=s.symbols('c t v'); xx=s.Matrix(x)
pairs=list(itertools.combinations(range(5),2))
q=sum(z*z for z in x)+2*c*x0*x2
qU=q-x3*x3-x4*x4
xi=s.Matrix([x0,x1,-2*x2,0,0])
R=s.Matrix([0,0,0,x4,-x3]); S=s.Matrix([x1,-x0,0,0,0])
J=s.zeros(5);J[3,4]=1;J[4,3]=-1
JS=s.zeros(5);JS[0,1]=1;JS[1,0]=-1

def wedge(a,b): return a*b.T-b*a.T

def jacobi(P):
 return {ijk:s.expand(sum(P[i,l]*s.diff(P[j,k],x[l])+
   P[j,l]*s.diff(P[k,i],x[l])+P[k,l]*s.diff(P[i,j],x[l])
   for l in range(5))) for ijk in itertools.combinations(range(5),3)
   for i,j,k in [ijk]}

P0=2*wedge(xi,R)+s.expand(sum(xi[i]*s.diff(qU,x[i]) for i in range(3)))*J
KR=q*J+wedge(xx,R)
P=P0+t*KR
qgrad=s.Matrix([s.diff(q,z) for z in x])
assert all(s.expand(w)==0 for w in P*qgrad)
assert all(z==0 for z in jacobi(P).values())
assert all(s.expand(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])==0
  for i,j,k,l in itertools.combinations(range(5),4))
assert all(s.expand(sum(s.diff(P[i,j],x[i]) for i in range(5))-3*t*R[j])==0 for j in range(5))

# so(q_c): q_c=x^T G_c x, A=G_c^(-1)Z, Z^T=-Z.
G=s.hessian(q,xx)/2
z=s.symbols('z0:10'); Z=s.zeros(5)
for (i,j),value in zip(pairs,z): Z[i,j]=value; Z[j,i]=-value
A=G.inv()*Z; V=A*xx
assert (A.T*G+G*A).applyfunc(s.cancel)==s.zeros(5)
rows=[]
for i,j in pairs:
 expr=s.expand(sum(V[k]*s.diff(P0[i,j],x[k]) for k in range(5))-
      sum(A[i,k]*P0[k,j]+P0[i,k]*A[j,k] for k in range(5)))
 rows.extend(s.Poly(expr,*x).coeffs())
M=s.linear_eq_to_matrix(rows,z)[0]
assert M.shape==(73,10)
assert M*s.eye(10)[:,9]==s.zeros(M.rows,1)
row_indices=[0,1,2,3,4,5,13,33,49]
col_indices=list(range(9))
minor=s.factor(M.extract(row_indices,col_indices).det(method='domain-ge'))
target=13824*c/(c*c-1)**4
assert s.cancel(minor-target)==0
assert M.subs(c,0).rank()==8
assert M.subs(c,0)*s.eye(10)[:,0]==s.zeros(M.rows,1)
assert M.subs(c,0)*s.eye(10)[:,9]==s.zeros(M.rows,1)
assert M.subs(c,s.Rational(1,2)).rank()==9
assert s.factor(minor.subs(c,s.Rational(1,2)))==s.Rational(65536,3)

# At c=0, the extra direction is K_S; Jacobi is v(t+2), exactly.
q0=q.subs(c,0);Pzero=P0.subs(c,0);K0=KR.subs(c,0)
KS=q0*JS+wedge(xx,S)
Pb=Pzero+t*K0+v*KS
assert all(s.expand(w)==0 for w in Pb*s.Matrix([s.diff(q0,zz) for zz in x]))
jb=jacobi(Pb)
nonzero={''.join(map(str,k)):str(s.factor(pol)) for k,pol in jb.items() if pol!=0}
assert list(nonzero)==['013','014','034','134']
assert nonzero=={
 '013':'2*v*x4*(t + 2)*(x0**2 + x1**2 + x2**2 + x3**2 + x4**2)',
 '014':'-2*v*x3*(t + 2)*(x0**2 + x1**2 + x2**2 + x3**2 + x4**2)',
 '034':'2*v*x1*(t + 2)*(x0**2 + x1**2 + x2**2 + x3**2 + x4**2)',
 '134':'-2*v*x0*(t + 2)*(x0**2 + x1**2 + x2**2 + x3**2 + x4**2)'
}
for pol in jb.values():
 assert s.rem(s.Poly(pol,t,v),s.Poly(v*(t+2),t,v))==0
pf=[]
for missing in range(5):
 i,j,k,l=[p for p in range(5) if p!=missing]
 pf.append((-1)**missing*s.expand(Pb[i,j]*Pb[k,l]-Pb[i,k]*Pb[j,l]+Pb[i,l]*Pb[j,k]))
f=s.factor(pf[0]/(2*x0))
assert s.expand(f.subs(t,-2)+3*v*x2*(x2*x2+x3*x3+x4*x4))==0
assert all(s.expand(pf[i]-2*x[i]*f)==0 for i in range(5))

result={
 'status':'PASS','sympy':s.__version__,
 'quadric':'q_c=sum(x_i^2)+2*c*x0*x2; smooth iff c^2 != 1',
 'intrinsic_source':'xi=x0*d0+x1*d1-2*x2*d2; fixed logarithmic boundary with x2 factor',
 'universal_rank_two_pencil':'P_c,t=P0_c+t*(q_c*J+E wedge R)',
 'universal_rank_two_Jacobi_and_Pfaffian':'symbolically zero over Q[c,t]',
 'strict_stabilizer_matrix_shape':list(M.shape),
 'strict_stabilizer_generic_rank':9,
 'strict_stabilizer_minor_rows':row_indices,
 'strict_stabilizer_minor_cols':col_indices,
 'strict_stabilizer_minor':str(minor),
 'stabilizer_off_wall':'C*R for c != 0, c^2 != 1',
 'wall':'c=0',
 'stabilizer_on_wall':'C*R+C*S, rank 8',
 'wall_extension_ideal':'(v*(t+2))',
 'wall_scheme':'Spec C[t,v]/(v(t+2))',
 'wall_nonzero_Jacobi':nonzero,
 'rank4_branch':'t=-2, v != 0; f=-3*v*x2*(x2^2+x3^2+x4^2)',
 'rank4_divisor':'q_0*x2*(x2^2+x3^2+x4^2)',
 'scope':'complete fibre-by-fibre Poisson extension for the one-parameter resonant Log boundary q_c with c^2!=1; not a classification of all exceptional Log points'
}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'PASS','minor':str(minor),'wall_fibre':result['wall_scheme']},ensure_ascii=False))
