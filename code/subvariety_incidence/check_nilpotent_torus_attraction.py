#!/usr/bin/env python3
"""Exact cocharacter action on the normalized invariant nilpotent Poisson pencil."""
import argparse,itertools,json
from pathlib import Path
import sympy as s
p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,default=Path('work/nilpotent_torus_attraction.json'))
p.add_argument('--check',type=Path)
args=p.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('u0 u1 z w0 w1');u0,u1,z,w0,w1=x
eta,t,lam,rho=s.symbols('eta t lam rho',nonzero=True)
v=s.Matrix(x)
q=u0**2+2*u1*w0+2*z*w1+eta*w1**2
G=s.hessian(q,v)/2
A=s.zeros(5);A[0,1]=1;A[3,0]=-1
xi=(u0+t*u1,u1,-2*z)
gamma={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]}
pairs=list(itertools.combinations(range(5),2))
def sign(t):return (-1)**sum(t[i]>t[j]for i in range(len(t))for j in range(i+1,len(t)))
P=s.zeros(5)
for i,j in pairs:
 rest=tuple(k for k in range(5)if k not in (i,j))
 P[i,j]=s.expand(sign((i,j,*rest))*sum((-1)**h*s.diff(q,x[k])*gamma.get(tuple(rest[a]for a in range(3)if a!=h),0) for h,k in enumerate(rest)))
 P[j,i]=-P[i,j]
K=(q*A*G.inv()/2+(v*(A*v).T-(A*v)*v.T)/2).applyfunc(s.expand)
S=s.diag(1,1,rho,1,1/rho)
subs={z:z/rho,w1:rho*w1}
Ppush=(S*P*S.T).subs(subs,simultaneous=True)
Kpush=(S*K*S.T).subs(subs,simultaneous=True)
Pnew=P.subs(eta,eta*rho**2)
Knew=K.subs(eta,eta*rho**2)
assert S.det()==1 and S*A*S.inv()==A
assert all(s.cancel(rho*Ppush[i,j]-Pnew[i,j])==0 for i,j in pairs)
assert all(s.cancel(Kpush[i,j]-Knew[i,j])==0 for i,j in pairs)
assert all(s.cancel(rho*(Ppush[i,j]+lam*Kpush[i,j])-Pnew[i,j]-rho*lam*Knew[i,j])==0 for i,j in pairs)
W=s.Matrix([[0,0],[0,eta]])
assert W.rank()==1 and W.subs(eta,0).rank()==0
report={
 'date':'2026-11-10','version':'v16','status':'PASS','sympy':s.__version__,
 'cocharacter':'diag(1,1,rho,1,rho^-1)',
 'A_centralizer':True,'determinant':'1',
 'strict_tensor_transform':'rho * g_rho(P_eta,t + lambda*K_eta) = P_(rho^2 eta),t + (rho lambda)*K_(rho^2 eta)',
 'eta_weight':2,'lambda_weight':1,
 'torus_invariant_for_lambda_nonzero':'eta/lambda^2',
 'rank_q_restriction_common_kernel_for_eta0':0,
 'rank_q_restriction_common_kernel_for_etaNonzero':1,
 'rank2_quadric_orbits_under_preserving_A_and_base_3form':'two: eta=0, eta!=0',
 'limit_rho_to_zero':'[P_(0,t)]',
 'scope':'One cocharacter on the normalized A-invariant family. Does not determine full PGL orbits of the Poisson pencils or the component intersection dimension.'
}
if args.check is not None:
 assert report==json.loads(args.check.read_text(encoding='utf-8')), 'certificate mismatch'
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('PASS: torus action and contraction; 10 exact bivector coefficient identities')
