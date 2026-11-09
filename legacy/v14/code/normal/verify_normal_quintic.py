#!/usr/bin/env python3
"""Exact certificates for the normal-quintic reduction.
Universal first-jet identities are verified in independent symbols; the
Fisher (a,b)=(1,1) example is reconstructed solely from its quintic.
This does not prove exhaustiveness of normal quintics.
"""
from itertools import combinations
import json
from pathlib import Path
import sympy as s

x=s.symbols('x0:5'); pairs=list(combinations(range(5),2))
q=s.symbols('q0:10')
Q=s.zeros(5)
for a,(i,j) in zip(q,pairs): Q[i,j]=a; Q[j,i]=-a
dq={ (l,k):s.Symbol(f'd{l}_{k}') for l in range(5) for k in range(10)}
def d(expr,l): return s.expand(sum(s.diff(expr,q[k])*dq[l,k] for k in range(10)))
def pfaff4(M,inds):
    a,b,c,e=inds
    return s.expand(M[a,b]*M[c,e]-M[a,c]*M[b,e]+M[a,e]*M[b,c])
def pfaff_vector(M):
    return s.Matrix([(-1)**i*pfaff4(M,[j for j in range(5) if j!=i]) for i in range(5)])
def permutation_sign(indices):
    return (-1)**sum(indices[i]>indices[j] for i in range(len(indices)) for j in range(i+1,len(indices)))
def contract_trivector(T):
    ans={I:0 for I in combinations(range(5),2)}
    for I,v in T.items():
        rem=tuple(i for i in range(5) if i not in I)
        ans[rem]+=permutation_sign(I+rem)*v
    return {I:s.expand(v) for I,v in ans.items()}
pf=pfaff_vector(Q)
assert all(s.expand(v)==0 for v in Q*pf)
div=s.Matrix([sum(d(Q[i,j],i) for i in range(5)) for j in range(5)])
jac={I:s.expand(sum(Q[I[0],l]*d(Q[I[1],I[2]],l)+Q[I[1],l]*d(Q[I[2],I[0]],l)+Q[I[2],l]*d(Q[I[0],I[1]],l) for l in range(5))) for I in combinations(range(5),3)}
D_wedge_Q={I:s.expand(div[I[0]]*Q[I[1],I[2]]-div[I[1]]*Q[I[0],I[2]]+div[I[2]]*Q[I[0],I[1]]) for I in jac}
dalpha={I:s.expand(d(pf[I[1]],I[0])-d(pf[I[0]],I[1])) for I in pairs}
bv=contract_trivector({I:jac[I]-D_wedge_Q[I] for I in jac})
assert all(s.expand(dalpha[I]-bv[I])==0 for I in pairs)
# Divergence of Q times a formal gradient; Hessian cancellation is independent.
f=s.symbols('f0:5')
hvars=s.symbols('h0:15'); h=s.zeros(5)
for a,(i,j) in zip(hvars,combinations(range(6),2)):
    # combinations range6 encode symmetric indices (i,j-1)
    h[i,j-1]=a; h[j-1,i]=a
lhs=sum(d(Q[i,j],i)*f[j]+Q[i,j]*h[i,j] for i in range(5) for j in range(5))
assert s.expand(lhs-(div.T*s.Matrix(f))[0])==0

# Fisher, "Invariants for the elliptic normal quintic", (38),(39), a=b=1.
a=[5*x[i]**2-10*x[(i-1)%5]*x[(i+1)%5]-2*x[(i-2)%5]*x[(i+2)%5] for i in range(5)]
b=[5*x[i]**2-4*x[(i-1)%5]*x[(i+1)%5]+10*x[(i-2)%5]*x[(i+2)%5] for i in range(5)]
P=s.Matrix([[0,a[3],b[1],-b[4],-a[2]],[-a[3],0,a[4],b[2],-b[0]],[-b[1],-a[4],0,a[0],b[3]],[b[4],-b[2],-a[0],0,a[1]],[a[2],b[0],-b[3],-a[1],0]])
pp=pfaff_vector(P)
F=s.expand(sum(x[i]*pp[i] for i in range(5))/5)
assert all(s.expand(s.diff(F,x[i])-pp[i])==0 for i in range(5))
assert all(s.expand(sum(s.diff(P[i,j],x[i]) for i in range(5)))==0 for j in range(5))
assert all(s.expand(sum(P[i,l]*s.diff(P[j,k],x[l])+P[j,l]*s.diff(P[k,i],x[l])+P[k,l]*s.diff(P[i,j],x[l]) for l in range(5)))==0 for i,j,k in combinations(range(5),3))

m2=[x[i]*x[j] for i in range(5) for j in range(i,5)]
mons6=sorted(s.Poly(s.expand(sum(x)**6),*x).monoms(), reverse=True)
row_index={(i,e):k for k,(i,e) in enumerate((i,e) for i in range(5) for e in mons6)}
def annihilator_matrix(F):
    grad=[s.diff(F,z) for z in x]
    rows={}
    for k,(i,j) in enumerate(pairs):
        for l,m in enumerate(m2):
            col=15*k+l
            for target,expr in [(i,m*grad[j]),(j,-m*grad[i])]:
                for e,c in s.Poly(expr,*x).terms():
                    if c:
                        ri=row_index[target,e]
                        rows.setdefault(ri,{})[col]=c
    return s.MutableSparseMatrix(1050,150,{(i,j):v for i,row in rows.items() for j,v in row.items()})
L=annihilator_matrix(F)
# Fraction field exact rank and kernel, not sampling.
ker=L.nullspace()
assert len(ker)==1
pvec=s.Matrix([s.Poly(P[i,j],*x).coeff_monomial(m) for i,j in pairs for m in m2])
assert L*pvec==s.zeros(1050,1)
assert s.Matrix.hstack(ker[0],pvec).rank()==1
fermat=sum(z**5 for z in x)
ferm_rank=annihilator_matrix(fermat).rank()
assert ferm_rank==150
# A normal essential3 cone is not excluded by normality alone; its K has dimension15.
cone=x[0]**5+x[1]**5+x[2]**5
cone_rank=annihilator_matrix(cone).rank()
assert cone_rank==135

# Pure resolution numerical consequences.
t=s.symbols('t'); n=s.symbols('n', integer=True)
num=1-5*t**4+5*t**6-t**10
h=s.cancel(num/(1-t)**3)
assert s.expand(h)==t**7+3*t**6+6*t**5+10*t**4+10*t**3+6*t**2+3*t+1
HP=s.expand(sum(c*s.prod(n-a+i for i in range(1,5))/24 for c,a in [(1,0),(-5,4),(5,6),(-1,10)]))
assert HP==40*n-100
res={
  'universal_signed_pfaffian_syzygies':5,
  'universal_BV_first_jet_identities':10,
  'universal_divergence_gradient_identity':True,
  'Fisher_a_b':[1,1],
  'Fisher_F':str(F),
  'Fisher_skew_matrix':[[str(P[i,j]) for j in range(5)] for i in range(5)],
  'Fisher_gradient_equals_signed_Pfaffians':True,
  'Fisher_divergence_zero':True,
  'Fisher_Jacobi_zero':True,
  'annihilator_matrix_shape':[1050,150],
  'Fisher_annihilator_rank':150-len(ker),
  'Fisher_annihilator_kernel_dimension':len(ker),
  'Fermat_quintic_annihilator_rank':ferm_rank,
  'essential3_normal_cone_annihilator_rank':cone_rank,
  'pure_resolution_shifts':[0,4,6,10],
  'pure_resolution_Betti':[1,5,5,1],
  'Hilbert_numerator':str(num),
  'h_vector':[1,3,6,10,10,6,3,1],
  'Hilbert_polynomial':str(HP),
  'degree':40,
  'arithmetic_genus':101,
  'canonical_shift':5,
  'claims_not_verified_by_computation':['classification/exhaustion of normal quintics','Buchsbaum-Eisenbud theorem','normality of arbitrary input quintic','scheme-level identification of Poisson space with a determinantal locus']
}
Path(__file__).with_suffix('.json').write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in res.items() if k not in ['Fisher_F','Fisher_skew_matrix','claims_not_verified_by_computation']},ensure_ascii=False,indent=2))
