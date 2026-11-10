#!/usr/bin/env python3
"""Exact source-membership certificate for the isotropic Log rank-four curve.

Verifies the original five-variable homogeneous bivector, its affine chart,
and an actual PGL5 (affine linear) change to the rank-one-long,
nonconstant-g source. No nonlinear Poisson coordinate changes are used.
"""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s

ap = argparse.ArgumentParser()
ap.add_argument('--output', type=Path, default=Path('work/rank4_long_nonconstant_membership.json'))
args = ap.parse_args()
assert s.__version__ == '1.14.0'
u0,h,Z,w0,w1 = s.symbols('u0 h Z w0 w1')
coords=(u0,h,Z,w0,w1)
e=s.Matrix(coords)
q=u0**2+2*h*w0+2*Z*w1
xi=(u0,h,-2*Z)
gamma={(1,2):u0,(0,2):-h,(0,1):-2*Z}

def sign(a):
    return (-1)**sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))

P=s.zeros(5)
for i,j in itertools.combinations(range(5),2):
    r=tuple(k for k in range(5) if k not in (i,j))
    P[i,j]=s.expand(sign((i,j,*r))*sum((-1)**k*s.diff(q,coords[t])*gamma.get(tuple(r[v]for v in range(3)if v!=k),0)for k,t in enumerate(r)))
    P[j,i]=-P[i,j]
G=s.hessian(q,e)/2
A=s.zeros(5)
A[0,1]=1
A[3,0]=-1
v=A*e
assert A.T*G+G*A==s.zeros(5)
K=(q*A*G.inv()/2+(e*v.T-v*e.T)/2).applyfunc(s.expand)
a=s.symbols('a',nonzero=True)
B=P+a*K

def jac(B,vs):
    n=len(vs)
    return [s.expand(sum(B[i,l]*s.diff(B[j,k],vs[l])+B[j,l]*s.diff(B[k,i],vs[l])+B[k,l]*s.diff(B[i,j],vs[l]) for l in range(n))) for i,j,k in itertools.combinations(range(n),3)]
assert jac(B,coords)==[0]*10
x,y,z,w=s.symbols('x y z w')
old=(x,y,z,w)
fs=s.Matrix([u0/h,-w0/h,Z/h,w1/h])
chart=(fs.jacobian(coords)*B*fs.jacobian(coords).T).subs({u0:x,h:1,Z:z,w0:-y,w1:w},simultaneous=True).applyfunc(s.expand)
assert jac(chart,old)==[0]*4
expected=s.zeros(4)
for (i,j),f in {(0,1):a*(y-x**2/2-z*w),(1,2):-6*z**2,(1,3):2*w*z-2*x**2+4*y,(2,3):6*z}.items():
    expected[i,j]=f
    expected[j,i]=-f
assert all(s.expand(val)==0 for val in chart-expected)
modular=[s.expand(sum(s.diff(chart[i,j],old[i]) for i in range(4))) for j in range(4)]
assert modular==[-a,-a*x+10*z,0,10]
N=s.zeros(5)
N[1,0]=-a; N[2,1]=-a; N[2,3]=10; N[4,0]=10
assert N.rank()==2 and (N*N).rank()==1 and N**3==s.zeros(5)
N0=N.subs(a,0)
assert N0.rank()==2 and N0*N0==s.zeros(5)

X,Y,T,W=s.symbols('X Y T W')
changed=s.Matrix([-x/a+10*z/a**2,y/a**2,z,w+10*x/a])
Jac=changed.jacobian(old)
assert s.factor(Jac.det())==-1/a**3
inverse={x:-a*X+10*T/a,y:a*a*Y,z:T,w:W+10*X-100*T/a**2}
actual=(Jac*chart*Jac.T).subs(inverse,simultaneous=True).applyfunc(lambda val:s.factor(val))

# The six-bracket source from
# research/classification/RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md.
g=6*T
U_z=s.S.Zero
U_w=60*T/a**2
d=s.diff(U_z,T)+s.diff(U_w,W)
Ksource=-T*W/a**2+50*T**2/a**4
F=s.expand(U_z*s.diff(Ksource,T)+U_w*s.diff(Ksource,W)+(d+1)*Ksource)
r=Y-X**2/2
formula={
 (0,1):-F-(1+d)*r,
 (0,2):U_z,
 (0,3):U_w,
 (1,2):X*U_z+r*s.diff(g,W)+s.diff(g*Ksource,W),
 (1,3):X*U_w-r*s.diff(g,T)-s.diff(g*Ksource,T),
 (2,3):g,
}
assert s.expand(U_z*s.diff(g,T)+U_w*s.diff(g,W)-d*g)==0
for (i,j), f in formula.items():
    assert s.cancel(actual[i,j]-f)==0, (i,j,actual[i,j],f)
assert jac(actual,(X,Y,T,W))==[0]*4
assert [s.cancel(sum(s.diff(actual[i,j],(X,Y,T,W)[i]) for i in range(4))) for j in range(4)]==[1,X,0,0]
qcanon=s.expand((x*x-2*y+2*z*w).subs(inverse,simultaneous=True))
Qcanon=Y-X**2/2+Ksource
assert s.cancel(qcanon+2*a*a*Qcanon)==0
hs=s.symbols('hs')
Qhom=hs*Y-X**2/2-(T*W)/a**2+50*T**2/a**4
QG=s.hessian(Qhom,s.Matrix([hs,X,Y,T,W]))/2
assert s.factor(QG.det())==-1/(32*a**4)
Pf=s.expand(actual[0,1]*actual[2,3]-actual[0,2]*actual[1,3]+actual[0,3]*actual[1,2])
assert s.cancel(Pf+g*(r+Ksource))==0
report={
 'date':'2026-11-10','version':'v11','status':'PASS','sympy':s.__version__,
 'component':'R4_MODULAR_RANK1_LONG_NONCONSTANT_25',
 'original_chart':['x=u0/h','y=-w0/h','z=Z/h','w=w1/h'],
 'original_chart_brackets':{''.join(str(k) for k in ij):str(chart[ij]) for ij in itertools.combinations(range(4),2)},
 'original_modular':[str(t) for t in modular],
 'modular_Jordan_a_nonzero':[3,1,1], 'modular_Jordan_a_zero':[2,2,1],
 'affine_linear_coordinates':['X=-x/a+10*z/a**2','Y=y/a**2','T=z','W=w+10*x/a'],
 'affine_change_det':str(s.factor(Jac.det())),
 'source_g':'6*T','source_U_z':'0','source_U_w':'60*T/a**2','source_divergence':'0',
 'source_K':'-T*W/a**2+50*T**2/a**4','source_F':str(F),
 'source_bracket_match':'all six exact over QQ(a)',
 'Jacobi':'all strict affine Jacobi coefficients zero',
 'homogeneous_P0_plus_aK1_Jacobi':'all ten strict Jacobi coefficients zero',
 'projective_quadric_formula':'q_original_in_new_chart=-2*a**2*(Y-X**2/2+K)',
 'quadric_Hessian_det':str(s.factor(QG.det())),
 'Pfaffian':'-6*T*(Y-X**2/2+K)',
 'conclusion':'For every a != 0 the projective Poisson tensor is PGL5-equivalent to an actual point of the long nonconstant-g source; by closure P0 is in the rank-four 25-dimensional component. No classification of all components through P0.'
}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
print('PASS rank-four curve in long nonconstant-g component, all a != 0; Q determinant',report['quadric_Hessian_det'])
