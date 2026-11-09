#!/usr/bin/env python3
"""Exact full-source checks for the shared conductor/conormal reconstruction.

No original producer or PASS receipt is imported. Checks linear completeness in
all 150 cone coefficients, the 84-parameter source, and full polynomial identities.
It does not verify the old componenthood or finite-case recognition arguments.
"""
from __future__ import annotations
import argparse
import hashlib
import platform
import itertools as it
import json
from pathlib import Path
import sympy as s

X=s.symbols('x0 x1 x2 u v');x=X[:3];u,v=X[3:];E=s.Matrix(X)
pairs=list(it.combinations(range(5),2));mons=[X[i]*X[j] for i in range(5) for j in range(i,5)]
checks=[]

def check(name,values):
    for value in values:
        if s.expand(value)!=0:
            raise AssertionError((name,s.expand(value)))
    checks.append(name)
    print(name,flush=True)

def wedge(a,b):
    a=s.Matrix(a);b=s.Matrix(b)
    return a*b.T-b*a.T

def div(p):
    return s.Matrix([sum(s.diff(p[i,j],X[i]) for i in range(5)) for j in range(5)])

def canonical(p):
    return (p-wedge(E,div(p))/5).applyfunc(s.expand)

def field(prefix,tracefree=False):
    p=s.symbols(prefix+':'+str(8 if tracefree else 9))
    mat=s.Matrix(3,3,list(p)+([-p[0]-p[4]] if tracefree else []))
    return mat*s.Matrix(x),p

def embed(y):return list(y)+[0,0]

def jac(p):
    return s.Matrix([s.expand(sum(p[i,l]*s.diff(p[j,k],X[l])+p[j,l]*s.diff(p[k,i],X[l])+p[k,l]*s.diff(p[i,j],X[l]) for l in range(5))) for i,j,k in it.combinations(range(5),3)])

def pf(p):
    out=[]
    for at in range(5):
        i,j,k,l=[i for i in range(5) if i!=at]
        out.append((-1)**at*(p[i,j]*p[k,l]-p[i,k]*p[j,l]+p[i,l]*p[j,k]))
    return s.expand(sum(X[i]*out[i] for i in range(5)))

def coeffs(p):
    return s.Matrix([s.Poly(p[i,j],X).coeff_monomial(m) for i,j in pairs for m in mons])

def chart(p,i,j):
    return s.expand((p[i,j]-X[i]*p[0,j]-X[j]*p[i,0]).subs(X[0],1))

def conditions(p):
    out=[]
    for i in (1,2):
        for j in (3,4):
            out.append(chart(p,i,j).subs({u:0,v:0}))
    q=chart(p,3,4)
    out += [q.subs({u:0,v:0}),s.diff(q,u).subs({u:0,v:0}),s.diff(q,v).subs({u:0,v:0})]
    return [c for f in out for c in s.Poly(f,x[1:]).coeffs()]

rho_par=s.symbols('r:10');cubics=[x[i]*x[j]*x[k] for i in range(3) for j in range(i,3) for k in range(j,3)]
rho=sum(c*m for c,m in zip(rho_par,cubics))
Q=s.zeros(5)
for a,b,k in [(1,2,0),(2,0,1),(0,1,2)]:
    Q[a,b]=s.diff(rho,x[k])/3;Q[b,a]=-Q[a,b]
T,tp=field('t',True);A,ap=field('a');B,bp=field('b');C,cp=field('c')
eu=[0,0,0,1,0];ev=[0,0,0,0,1]
tr=lambda V:sum(s.diff(V[i],x[i]) for i in range(3))
P0=Q+wedge(embed(s.Matrix(x)),embed(T))+wedge(embed(u*A+v*B),eu)+wedge(embed(u*C+v*(3*T-A)),ev)+(-tr(C)*u*u/2+tr(A)*u*v+tr(B)*v*v/2)*wedge(eu,ev)
P0=P0.applyfunc(s.expand);params0=list(rho_par+tp+ap+bp+cp)
lp=s.symbols('l:18');raw1=s.zeros(5)
for at,(i,j) in enumerate(it.combinations(range(3),2)):
    raw1[i,j]=u*sum(lp[3*at+k]*x[k] for k in range(3))+v*sum(lp[9+3*at+k]*x[k] for k in range(3));raw1[j,i]=-raw1[i,j]
vp=s.symbols('n:12')
raw1+=wedge(embed(s.Matrix(vp[:3])*u*u),ev)+wedge(embed(s.Matrix(vp[3:6])*u*v),eu)+wedge(embed(s.Matrix(vp[6:9])*u*v),ev)+wedge(embed(s.Matrix(vp[9:])*v*v),eu)
P1=canonical(raw1);params1=list(lp+vp)
wp=s.symbols('w:9');P2=s.zeros(5)
for at,(i,j) in enumerate(it.combinations(range(3),2)):
    P2[i,j]=sum(wp[3*at+k]*m for k,m in enumerate((u*u,u*v,v*v)));P2[j,i]=-P2[i,j]
P=P0+P1+P2;params=params0+params1+list(wp)
check('all_84_source_columns_divergence_zero',div(P))
check('all_84_source_columns_plane_invariant_abelian',conditions(P))
check('restriction_is_rho',[chart(P,1,2).subs({u:0,v:0})-rho.subs(x[0],1)])
source=coeffs(P).jacobian(params)
assert source.rank()==84;checks.append('source_rank_45_30_9_total_84')
raw_par=s.symbols('g:150');generic=s.zeros(5)
for at,(i,j) in enumerate(pairs):
    generic[i,j]=sum(raw_par[15*at+k]*m for k,m in enumerate(mons));generic[j,i]=-generic[i,j]
constraints=s.Matrix(conditions(generic)).jacobian(raw_par)
assert constraints.rank()==42
check('full150_constraint_matrix_annihilates_source',constraints*source)
ypar=s.symbols('y:25');Y=s.Matrix(5,5,ypar)*E
gauge=coeffs(wedge(E,Y)).jacobian(ypar)
assert gauge.rank()==24
assert source.row_join(gauge).rank()==108
check('full150_constraints_annihilate_all_Euler_gauge',constraints*gauge)
checks.append('full150_kernel_108_equals_source84_plus_Euler24')
for row,norm in [(0,3),(1,4)]:
    for col,z in enumerate((u,v)):
        got=s.diff(chart(P0-Q,1,norm),z).subs({u:0,v:0})
        V=[[A+T,B],[C,4*T-A]][row][col]
        want=(V[1]-x[1]*V[0]).subs(x[0],1)
        got2=s.diff(chart(P0-Q,2,norm),z).subs({u:0,v:0})
        want2=(V[2]-x[2]*V[0]).subs(x[0],1)
        check(f'conormal_matrix_entry_{row}_{col}',[got-want,got2-want2])
alpha,beta,gamma=s.symbols('alpha beta gamma')
hidden=P0.subs({**dict(zip(tp,[0]*8)),**dict(zip(ap,list(alpha*s.eye(3)))),**dict(zip(bp,list(beta*s.eye(3)))),**dict(zip(cp,list(gamma*s.eye(3)))),**dict(zip(rho_par,[0]*10))})
K=-s.Rational(5,2)*gamma*u*u+5*alpha*u*v+s.Rational(5,2)*beta*v*v
check('hidden_three_parameter_binary_quadratic',list(hidden-canonical(K*wedge(eu,ev))))
for k,Pk in enumerate((P0,P1,P2)):
    check(f'normal_weight_{k}',[u*s.diff(Pk[i,j],u)+v*s.diff(Pk[i,j],v)-(k+int(i>=3)+int(j>=3))*Pk[i,j] for i,j in pairs])
# Full polynomial source, not selected coefficient samples.
D=pf(P)
assert all(sum(m[3:]) in (2,3,4,5) for m,c in s.Poly(D,X).terms())
assert all(sum(m[:3])==0 for m,c in s.Poly(D,X).terms() if sum(m[3:])==5)
checks.append('full84_actual_Pfaffian_only_normal_degrees_2_to_5')
check('weight_two_self_Jacobi_zero',jac(P2))
J=jac(P)
for at,(i,j,k) in enumerate(it.combinations(range(5),3)):
    assert all(0<=m[3]+m[4]-sum(a>=3 for a in (i,j,k))<=3 for m,c in s.Poly(J[at],X).terms())
checks.append('full84_Jacobi_only_weights_0_to_3')
# det-twisted Sym^2 has invariant discriminant, including its zero locus.
aa,bb,cc=s.symbols('aa bb cc');m00,m01,m10,m11=s.symbols('m00 m01 m10 m11')
q=aa*u*u+bb*u*v+cc*v*v
act=s.expand(s.diff(q,u)*(m00*u+m01*v)+s.diff(q,v)*(m10*u+m11*v)-(m00+m11)*q)
da=s.expand(act).coeff(u,2).coeff(v,0);db=s.expand(act).coeff(u,1).coeff(v,1);dc=s.expand(act).coeff(v,2).coeff(u,0)
check('twisted_quadratic_discriminant_GL2_invariant',[2*bb*db-4*cc*da-4*aa*dc])
local=s.zeros(5);local[0,1]=1;local[1,0]=-1;local[3,4]=v*v-u**3;local[4,3]=-local[3,4]
check('local_cusp_is_Poisson',jac(local))
a=s.symbols('normalization_parameter')
check('cusp_residue_has_order_two',[s.diff(a*a,a)/(2*a**3)-1/a**2])
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
result={'version':'v70-conductor-conormal-reconstruction','status':'PASS','python':platform.python_version(),'sympy':s.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sampling_used':False,'imports_old_producers':False,'checks':checks,'source_dimensions':{'rank_two_or_unrestricted':[45,30,9],'rank_zero':[35,30,9],'full150_constraints_rank':42,'Euler_rank':24},'limits':['Not an independent verification of the 147 components or old finite recognition cases.','The general algebroid and leading-symbol transport proofs are mathematical arguments in the shared TeX, not certified by these polynomial checks.']}
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
