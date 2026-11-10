"""Exact finite checks for the invariant resonant Log orthogonality criterion.

This does not replace the general geometric proof in quadric_log_extension.tex.
Run: python code/subvariety_incidence/check_log_orthogonality.py --output work/log_orthogonality.json
"""
from __future__ import annotations
import argparse
from itertools import combinations
import json
from pathlib import Path
import sympy as s

arg=argparse.ArgumentParser(description=__doc__)
arg.add_argument('--output',type=Path,default=Path('work/log_orthogonality.json'))
args=arg.parse_args()
assert s.__version__ == '1.14.0'
x=s.symbols('x0:5'); xx=s.Matrix(x); u0,u1,z,w0,w1=x
pairs=list(combinations(range(5),2))
xi=s.Matrix([u0,u1,-2*z,0,0]); R=s.Matrix([0,0,0,w1,-w0]); E=xx
J=s.zeros(5);J[3,4]=1;J[4,3]=-1

def wedge(v,w):return v*w.T-w*v.T

def jacobi(p):
 return [s.expand(sum(p[i,h]*s.diff(p[j,k],x[h])+p[j,h]*s.diff(p[k,i],x[h])+p[k,h]*s.diff(p[i,j],x[h]) for h in range(5))) for i,j,k in combinations(range(5),3)]

def run_case(b0,b1,c,expected_stab):
 qU=u0**2+u1**2+2*(b0*u0+b1*u1)*z+c*z*z
 q=qU+w0*w0+w1*w1
 G=s.hessian(q,xx)/2
 assert G.det()!=0
 g=s.expand(sum(xi[i]*s.diff(qU,x[i]) for i in range(3)))
 P=2*wedge(xi,R)+g*J
 KR=q*J+wedge(E,R)
 pt=P+s.Symbol('t')*KR
 assert all(s.expand(f)==0 for f in pt*s.Matrix([s.diff(q,a) for a in x]))
 assert all(f==0 for f in jacobi(pt))
 assert all(s.expand(sum(s.diff(P[i,j],x[i]) for i in range(5)))==0 for j in range(5))
 assert all(s.expand(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])==0 for i,j,k,l in combinations(range(5),4))
 # Full strict so(q) stabilizer of the affine homogeneous quadratic P.
 vvar=s.symbols('a0:10'); Z=s.zeros(5)
 for (i,j),a in zip(pairs,vvar):Z[i,j]=a;Z[j,i]=-a
 A=G.inv()*Z; V=A*xx
 expr=[]
 for i,j in pairs:
  pol=s.expand(sum(V[h]*s.diff(P[i,j],x[h]) for h in range(5))-
               sum(A[i,h]*P[h,j]+P[i,h]*A[j,h] for h in range(5)))
  expr += s.Poly(pol,*x).coeffs()
 M=s.linear_eq_to_matrix(expr,vvar)[0]
 rank=M.rank()
 assert rank==10-expected_stab,(b0,b1,c,rank,expected_stab)
 # Display natural Lie generators for qU.
 B=s.Matrix([[0,-1],[1,0]])
 if b0==0 and b1==0:
  block=s.zeros(5)
  block[:2,:2]=B
  assert ((block.T*G+G*block).applyfunc(s.simplify)==s.zeros(5))
  assert s.simplify(block*xx-xi) != s.zeros(5,1)
  assert s.simplify(M*s.Matrix([(-1 if (i,j)==(0,1) else 0) for i,j in pairs]))==s.zeros(M.rows,1)
 # for b!=0, its q-orthogonal Lie candidate has no nonzero B in so2 when H and qU are nondegenerate.
 return {'b0':str(b0),'b1':str(b1),'c':str(c),'strict_so5_rank':rank,'strict_stabilizer_dimension':10-rank,
         'qU_determinant':str(s.factor(G[:3,:3].det())),'g_determinant':str(s.factor((s.hessian(g,s.Matrix(x[:3]))/2).det())),
         'jacobi_with_arbitrary_t':True,'q_casimir_with_arbitrary_t':True,
         'universal_rank2_pfaffian':True}

cases=[
 (s.Integer(0),s.Integer(0),s.Integer(1),2),
 (s.Integer(1),s.Integer(0),s.Integer(2),1),
 (s.Integer(1),s.I,s.Integer(2),1),
 (s.Integer(1),s.Integer(1),s.Integer(3),1),
 (s.Integer(1),s.Integer(0),s.Integer(0),1),
]
result={'status':'PASS','sympy':s.__version__,
 'scope':'Finite exact stabilizer/Jacobi checks at five normal forms. General result proved in current TeX by zero-scheme geometry; no completeness inference from examples.',
 'proof_idea':'Zero(P0) = P(W) disjoint union plane conic (xi(qU)=0), preserving it forces block diagonal Lie symmetries; centralizer of xi in so(qU) is so(H) iff b=0 (when det H and det qU are nonzero).',
 'samples':[run_case(*v) for v in cases]}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(result,indent=2,ensure_ascii=False))
