#!/usr/bin/env python3
"""Exact torus weight decomposition on the 30-dimensional projective tangent."""
import argparse,itertools as it,json
from collections import Counter
from pathlib import Path
import sympy as s
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=Path('work/nilpotent_tangent_weights.json'))
parser.add_argument('--check',type=Path)
args=parser.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('u0 u1 z w0 w1');u0,u1,z,w0,w1=x;E=s.Matrix(x)
rho=s.symbols('rho',nonzero=True)
w=(0,0,1,0,-1)
S=s.diag(1,1,rho,1,1/rho);sub={z:z/rho,w1:rho*w1}
pairs=list(it.combinations(range(5),2))
mons=list(it.combinations_with_replacement(range(5),2))
def sign(a):return (-1)**sum(a[i]>a[j]for i in range(len(a))for j in range(i+1,len(a)))
def biv(q,xi):
 gam={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]}
 B=s.zeros(5)
 for i,j in pairs:
  rest=tuple(k for k in range(5)if k not in(i,j))
  B[i,j]=s.expand(sign((i,j,*rest))*sum((-1)**h*s.diff(q,x[k])*gam.get(tuple(rest[t]for t in range(3)if t!=h),0)for h,k in enumerate(rest)))
  B[j,i]=-B[i,j]
 return B
q=u0**2+2*u1*w0+2*z*w1
G=s.hessian(q,E)/2
P=biv(q,(u0,u1,-2*z))
def K(A):return(q*A*G.inv()/2+(E*(A*E).T-(A*E)*E.T)/2).applyfunc(s.expand)
def lie(M,B):
 v=M*E
 return s.Matrix([[s.expand(sum(v[k]*s.diff(B[i,j],x[k])-M[i,k]*B[k,j]-M[j,k]*B[i,k]for k in range(5)))for j in range(5)]for i in range(5)])
idx=[0,1,2,3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,22,23]
Ms=[]
for k in idx:
 M=s.zeros(5);M[k//5,k%5]=1;Ms.append(M)
A1=s.zeros(5);A1[0,1]=1;A1[3,0]=-1
A2=s.zeros(5);A2[3,2]=1;A2[4,1]=-1
A3=s.diag(0,1,-1,-1,1)
Bs=[lie(M,P)for M in Ms]
Bs +=[biv(f,(u0,u1,-2*z))for f in (u0*w1,w0**2,w0*w1,w1**2)]
Bs +=[biv(q,zv)for zv in ((u1,0,0),(0,u0,0),(u0,0,-z))]
Bs +=[K(A)for A in (A1,A2,A3)]
weights=[w[k//5]-w[k%5]for k in idx]+[1,0,1,2]+[0,0,0]+[1,0,1]
assert len(Bs)==len(weights)==30
for B,wt in zip(Bs,weights):
 transformed=(S*B*S.T).subs(sub,simultaneous=True)
 for i,j in pairs:
  assert s.cancel(rho*transformed[i,j]-rho**wt*B[i,j])==0
tangent=Counter(weights)
raw=Counter(1+w[i]+w[j]-w[a]-w[b]for i,j in pairs for a,b in mons)
gauge=Counter(1+w[i]-w[j]for i in range(5)for j in range(5))
gauge[1]-=1
scale=Counter({0:1})
kernel=tangent+gauge+scale
normal=raw-kernel
assert tuple(map(sum,(raw.values(),gauge.values(),tangent.values(),kernel.values(),normal.values())))==(150,24,30,55,95)
s0,s1,s2,s3,s4,s5,s6,a,b,c=weights[20:]
initial=[s0+c,s1+c,s2+c,s3+c,s4+c,s5+a,s5+c]
assert initial==[2,1,2,3,1,1,1] and s0+c==s2+a and s4+c==s6+a and s2==c
report={
 'date':'2026-11-10','version':'v16','status':'PASS','sympy':s.__version__,
 'cocharacter':'diag(1,1,rho,1,rho^-1)',
 'normalized_action':'rho times pushforward, so P0 is fixed',
 'weights':[-2,-1,0,1,2,3,4],
 'ambient_quadratic_bivector_weight_dimensions':[raw[k]for k in range(-2,5)],
 'euler_gauge_weight_dimensions':[gauge[k]for k in range(-2,5)],
 'projective_scale_weight_dimensions':[int(k==0)for k in range(-2,5)],
 'projective_tangent_weight_dimensions':[tangent[k]for k in range(-2,5)],
 'linearization_normal_weight_dimensions':[normal[k]for k in range(-2,5)],
 'explicit_30_tangent_basis_weights':weights,
 'seven_initial_quadric_weights':initial,
 'nilpotent_h_weight':a+s2,
 'exact_eigenvector_checks':30,
 'claim':'Weight decomposition at P0 and homogeneous quadratic obstruction equations. Higher initial equations and global intersection dimension remain open.'
}
if args.check is not None:
 assert report==json.loads(args.check.read_text(encoding='utf-8')), 'certificate mismatch'
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('PASS: all 30 exact eigenvectors and 95 normal directions by torus weights')
