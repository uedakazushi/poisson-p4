#!/usr/bin/env python3
"""Exact universal identities used by the complete rank-two reconstruction.

No legacy producer or PASS receipt is imported. All spinor inputs, the full
Aff32 b/t family, positive-zero-divisor sources and distinguishing invariants
are rebuilt from the displayed mathematical formulas.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import sympy as s

checks=[]
def check(name,values):
    for v in values:
        if s.cancel(s.expand(v))!=0:raise AssertionError((name,s.expand(v)))
    checks.append(name);print(name,flush=True)
def skew(vals):
    p=s.zeros(5)
    for (i,j),v in zip(it.combinations(range(5),2),vals):p[i,j]=v;p[j,i]=-v
    return p
def div(p,x):return [s.expand(sum(s.diff(p[i,j],x[i]) for i in range(5))) for j in range(5)]
def jac(p,x):return {ijk:s.expand(sum(p[i,l]*s.diff(p[j,k],x[l])+p[j,l]*s.diff(p[k,i],x[l])+p[k,l]*s.diff(p[i,j],x[l]) for l in range(5))) for ijk in it.combinations(range(5),3) for i,j,k in [ijk]}
def pfive(p):
    return [s.expand((-1)**at*(p[i,j]*p[k,l]-p[i,k]*p[j,l]+p[i,l]*p[j,k])) for at in range(5) for i,j,k,l in [[q for q in range(5) if q!=at]]]
def actual(p,x):return s.expand(sum(a*b for a,b in zip(x,pfive(p))))
def wedge(a,b):return a*b.T-b*a.T

a,b,c,d=s.symbols('a b c d');u=(a,b,c,d)
Pgeo=skew((a*a,c*c,a*b+c*d,a*c,-a*b+c*d,d*d,a*d,b*b,b*c,b*d))
S=s.Matrix([[d,b,0,0],[-c,0,b,0],[0,-a,0,d],[0,0,-a,-c]])
v=s.Matrix([b,-d,c,-a]);T=S.col_join(s.Matrix([[a,c,d,b]]))
check('two_plane_resolution_compositions',list(S*v)+list(T*v)+list(T.T*Pgeo))
check('geometric_model_all_five_Pfaffians',pfive(Pgeo))
M=s.zeros(5);M[0,3]=M[3,0]=1;M[1,2]=M[2,1]=-1;M[4,4]=2
check('geometric_model_orthogonal_metric',list(Pgeo.T*M*Pgeo))
assert M.det()==2
mon=(a*a,a*b,a*c,a*d,b*b,b*c,b*d,c*c,c*d,d*d)
assert abs(s.Matrix([[s.Poly(Pgeo[i,j],*u).coeff_monomial(m) for m in mon] for i,j in it.combinations(range(5),2)]).det())==2
checks.append('complete_Veronsese_quadrics_basis_determinant_two')
Psp=skew((b*d,-b*c,-a*d,a*c,b*b,d*d,-a*b-c*d,a*b-c*d,c*c,a*a))
G=s.zeros(5);G[0,0]=-2;G[1,4]=G[4,1]=-1;G[2,3]=G[3,2]=1
check('spinor_Pfaffians_and_metric',pfive(Psp)+list(Psp.T*G*Psp));assert G.det()==-2
coeff=s.symbols('A:20');A=s.Matrix(4,5,coeff)
for row,rank in [(None,20),(0,16),(1,12)]:
    p=Psp.copy()
    if row is not None:
        for k in range(5):p[row,k]=p[k,row]=0
    dz=[sum(s.diff(p[i,j],u[k])*A[k,i] for i in range(5) for k in range(4)) for j in range(5)]
    equations=[s.expand(z).coeff(t) for z in dz for t in u]
    mat,_=s.linear_eq_to_matrix(equations,coeff);assert mat.rank()==rank
    if row is None:assert abs(mat.det())==16
    if row==0:
        for col in mat.nullspace():
            AA=s.Matrix(4,5,col)
            assert all(AA[i,j]==0 for i in range(4) for j in range(1,5))
    checks.append(f'spinor_divergence_row_{row}_rank_{rank}')
x=s.symbols('x0:5');r=s.symbols('r0:8');t=x[1];z=x[4]
forms=(r[0]*t-r[4]*z,-r[6]*x[0]+r[1]*t-r[4]*x[2]+r[2]*z,r[3]*t+r[6]*z,r[4]*x[0]+r[5]*t+r[6]*x[3]+r[7]*z)
P=Psp.subs(dict(zip(u,forms)));P[1,:]=s.zeros(1,5);P[:,1]=s.zeros(5,1)
Q1=r[0]*x[0]*t-r[3]*t*x[3]-r[4]*x[0]*z-r[5]*t*z-r[6]*x[3]*z-s.Rational(1,2)*r[7]*z*z
Q2=r[0]*t*x[2]+r[1]*t*z+s.Rational(1,2)*r[2]*z*z-r[3]*x[0]*t-r[4]*x[2]*z-r[6]*x[0]*z
check('full_eight_parameter_spinor_divergence',div(P,x))
for i,j in it.combinations(range(5),2):
    triple=[k for k in range(5) if k not in (i,j)]
    grad=s.Matrix([[s.diff(f,x[k]) for k in triple] for f in (t,Q1,Q2)])
    check(f'full_eight_parameter_Jacobian_pair_{i}{j}',[grad.det()+(-1)**(i+j-1)*P[i,j]])
check('full_eight_parameter_spinor_Jacobi_and_rank',list(jac(P,x).values())+pfive(P))

# Entire Aff32 actual family, in b and t, without choosing parameter values.
beta,tau=s.symbols('beta tau');N=s.Matrix([0,x[0],x[1],0,x[3]]);H=s.Matrix([0,x[1],2*x[2],beta*x[3],(beta+1)*x[4]])
eta=s.zeros(5);theta=s.zeros(5)
for ij,val in {(0,1):-(2*beta*beta+beta-1)*x[0]*x[1],(0,2):-(beta*beta+beta-1)*x[0]*x[2]-beta*x[1]**2,(0,3):-(2*beta-1)*x[0]*x[3],(0,4):-(beta-1)*x[1]*x[3],(1,2):(beta*beta-4*beta+2)*x[1]*x[2],(2,4):(beta-2)*x[2]*x[4]}.items():
    i,j=ij;eta[i,j]=val;eta[j,i]=-val
for ij,val in {(0,1):beta*(beta-1)*(beta+1)*x[1]**2,(0,2):beta*(beta-1)**2*x[1]*x[2],(0,3):beta*(beta-1)*x[1]*x[3]}.items():
    i,j=ij;theta[i,j]=val;theta[j,i]=-val
family=wedge(N,H)+tau*eta+tau*tau*theta;J=jac(family,x)
check('full_Aff32_beta_tau_projective_Jacobi',[sum((-1)**at*x[k]*J[ids[:at]+ids[at+1:]] for at,k in enumerate(ids)) for ids in it.combinations(range(5),4)])
F=x[3]*(2*(x[0]-(beta-1)*tau*x[1])*x[2]-x[1]**2)*(x[1]*x[3]-(x[0]-beta*tau*x[1])*x[4])
check('full_Aff32_beta_tau_actual_Pfaffian',[actual(family,x)-(beta-1)*tau*F])
check('Aff32_center_exact',[q for q in family.subs(tau,0)-wedge(N,H)])

# Entire positive divisor source H*N^v_A, including H(N)=0 and all Jordan A.
hc=s.symbols('h0:5');mc=s.symbols('m0:16');HH=sum(a*b for a,b in zip(hc,x));YY=s.Matrix(list(s.Matrix(4,4,mc)*s.Matrix(x[:4]))+[0]);NN=s.Matrix([0,0,0,0,1])
pos=HH*wedge(NN,YY)
check('all_H_A_positive_zero_divisor_Poisson',list(jac(pos,x).values())+pfive(pos))
# Weighted logarithmic residue matrices, with arbitrary three parameters.
p,q,vv=s.symbols('p q v');lam=s.zeros(4)
for i,j,val in [(0,1,4*p),(0,2,4*q),(1,2,4*vv),(0,3,-2*(p+q)),(1,3,2*(p-vv)),(2,3,2*(q+vv))]:lam[i,j]=val;lam[j,i]=-val
check('all_horizontal_Log1112_residues',list(lam*s.Matrix([1,1,1,2]))+[lam[0,1]*lam[2,3]-lam[0,2]*lam[1,3]+lam[0,3]*lam[1,2]])
# No projective tangent linear field at the Rat122 distinguishing point.
mx=s.symbols('X:25');field=s.Matrix(5,5,mx)*s.Matrix(x);L=x[0];Q=sum(z*z for z in x);R=sum((i+1)*z*z for i,z in enumerate(x))
exprs=[s.expand(L*sum(field[i]*s.diff(g,x[i]) for i in range(5))-2*g*field[0]) for g in (Q,R)]
mat,_=s.linear_eq_to_matrix([v for f in exprs for v in s.Poly(f,*x).coeffs()],mx)
assert mat.rank()==24
checks.append('Rat122_tangent_linear_field_only_Euler')
for partition,weights in [('5',(0,1,2,3,4)),('41',(0,1,2,3,7))]:
    NN=s.Matrix([0,x[0],x[1],x[2],x[3] if partition=='5' else 0]);HH=s.Matrix([weights[i]*x[i] for i in range(5)])
    minors=[s.Matrix.hstack(s.Matrix(x),NN,HH).extract(ids,range(3)).det() for ids in it.combinations(range(5),3)]
    assert s.gcd_list(minors)==1
    checks.append('Aff'+partition+'_no_divisorial_zero_witness')

parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
result={'version':'v72-rank-two-reconstruction','status':'PASS','checks':checks,'check_count':len(checks),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy':s.__version__,'sampling_used':False,'legacy_producers_imported':False,'limits':['This checker proves displayed algebraic identities, not the reflexive sheaf geometry or the external foliation classifications.','The six full-ambient component matrices are separately rerun with the original code and recorded in reports/rank_two_component_recheck.json.','No independent verification of all147 components.']}
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,indent=2))
