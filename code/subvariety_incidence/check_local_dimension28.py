#!/usr/bin/env python3
"""Independent exact local-dimension certificate at the isotropic Log Poisson point.

All arithmetic uses QQ (SymPy 1.14). No network or historical assets.
The full linearization rank 95 is separately checked in
check_ambient_quadratic_obstruction.py.
"""
import argparse
import itertools as it
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=Path('work/local_dimension28.json'))
args=parser.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('u0 u1 z w0 w1')
u0,u1,z,w0,w1=x
E=s.Matrix(x)
pairs=list(it.combinations(range(5),2))
triples=list(it.combinations(range(5),3))
quads=list(it.combinations(range(5),4))
mon2=list(it.combinations_with_replacement(range(5),2))
mon3=list(it.combinations_with_replacement(range(5),3))

def sign(t):
    return (-1)**sum(t[i]>t[j] for i in range(len(t)) for j in range(i+1,len(t)))

def bivector(q,xi):
    gamma={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]}
    B=s.zeros(5)
    for i,j in pairs:
        rest=tuple(k for k in range(5) if k not in (i,j))
        B[i,j]=s.expand(sign((i,j,*rest))*sum((-1)**r*s.diff(q,x[k])*gamma.get(tuple(rest[t] for t in range(3) if t!=r),0) for r,k in enumerate(rest)))
        B[j,i]=-B[i,j]
    return B

def jac(B):
    return {abc:s.expand(sum(B[i,l]*s.diff(B[j,k],x[l])+B[j,l]*s.diff(B[k,i],x[l])+B[k,l]*s.diff(B[i,j],x[l]) for l in range(5))) for abc in triples for i,j,k in [abc]}

def project_jac(B):
    J=jac(B)
    return {I:s.expand(sum((-1)**j*x[k]*J[tuple(t for t in I if t!=k)] for j,k in enumerate(I))) for I in quads}

def coefficient(F,I,mon):
    return s.Poly(F[I],*x).coeff_monomial(mon)

q=u0**2+2*u1*w0+2*z*w1
xi=s.Matrix([u0,u1,-2*z])
P=bivector(q,xi)
G=s.hessian(q,E)/2
assert G.det()==1 and not any(jac(P).values())
A1=s.zeros(5); A1[0,1]=1; A1[3,0]=-1
A2=s.zeros(5); A2[3,2]=1; A2[4,1]=-1
A3=s.diag(0,1,-1,-1,1)

def K(A):
    C=A*E
    assert A.T*G+G*A==s.zeros(5)
    return (q*A*G.inv()/2+(E*C.T-C*E.T)/2).applyfunc(s.expand)
K1,K2,K3=[K(A) for A in (A1,A2,A3)]
a=s.symbols('a')
Pa=P+a*K1
assert not any(jac(Pa).values())
for omitted in range(5):
    i,j,k,l=[v for v in range(5) if v!=omitted]
    pf=Pa[i,j]*Pa[k,l]-Pa[i,k]*Pa[j,l]+Pa[i,l]*Pa[j,k]
    assert s.expand((-1)**omitted*pf-s.Rational(3,2)*a*u1**2*z*s.diff(q,x[omitted]))==0

# Logarithmic two-form with residues in wedge^2 ker(1,1,1,2).
alpha={}
gamma={(1,2):u0,(0,2):-u1,(0,1):-2*z}
dq=[s.diff(q,v) for v in x]
for i,j,k in triples:
    alpha[i,j,k]=s.expand(dq[i]*gamma.get((j,k),0)-dq[j]*gamma.get((i,k),0)+dq[k]*gamma.get((i,j),0))
log1=s.Matrix([1/u0,-1/u1,0,0,0])
log2=3*s.Matrix(dq)/q-s.Matrix([0,4/u1,2/z,0,0])
for i,j in pairs:
    omega=s.expand(sum(sign((k,i,j))*x[k]*alpha[tuple(sorted((k,i,j)))] for k in range(5) if k not in (i,j)))
    logarithmic=s.cancel(q*u0*u1*z*(log1[i]*log2[j]-log1[j]*log2[i]))
    assert s.cancel(omega-logarithmic)==0

# Two small universal left annihilators of the 150-column Jacobi linearization.
rows1=[((0,1,2,4),u0*z**2*w1,1),((0,1,2,4),u0**3*z,-1),((0,1,2,4),u0*u1*z*w0,1)]
rows2=[((0,1,3,4),u0*u1**2*w0,1),((0,1,2,4),u0*u1**2*z,-s.Rational(2,3))]

def obstruction(B,rows):
    F=project_jac(B)
    return s.expand(sum(weight*coefficient(F,I,mon) for I,mon,weight in rows))

for i,j in pairs:
    for aa,bb in mon2:
        R=s.zeros(5);R[i,j]=x[aa]*x[bb];R[j,i]=-R[i,j]
        Fp=project_jac(P+R)
        Fm=project_jac(P-R)
        for rows in (rows1,rows2):
            lin=sum(weight*(coefficient(Fp,I,mon)-coefficient(Fm,I,mon))/2 for I,mon,weight in rows)
            assert s.expand(lin)==0

# A full 30-dimensional tangent frame modulo Euler gauge and projective scale.
Mindices=[0,1,2,3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,22,23]
Ms=[]
for i in range(5):
    for j in range(5):
        M=s.zeros(5);M[i,j]=1;Ms.append(M)

def lie(M,B):
    v=M*E
    return s.Matrix([[s.expand(sum(v[k]*s.diff(B[i,j],x[k]) for k in range(5))-sum(M[i,k]*B[k,j]+M[j,k]*B[i,k] for k in range(5))) for j in range(5)]for i in range(5)])

def gauge(M):
    v=M*E
    return E*v.T-v*E.T

Bdirections=[lie(Ms[k],P)for k in Mindices]
qmon=[u0*w1,w0**2,w0*w1,w1**2]
Bdirections.extend(bivector(f,xi) for f in qmon)
xi3=[]
for i,j in it.product(range(3),repeat=2):
    if i!=j:
        M=s.zeros(3);M[i,j]=1;xi3.append(M)
for i in (0,1):
    M=s.zeros(3);M[i,i]=1;M[2,2]=-1;xi3.append(M)
Bdirections.extend(bivector(q,xi3[k]*s.Matrix(x[:3])) for k in (0,2,6))
Bdirections.extend((K1,K2,K3))
assert len(Bdirections)==30

def bivector_coordinates(B):
    return [s.Poly(B[i,j],*x).coeff_monomial(x[ii]*x[jj]) for i,j in pairs for ii,jj in mon2]

coords=s.Matrix.hstack(*[s.Matrix(bivector_coordinates(B)) for B in [gauge(M) for M in Ms]+[P]+Bdirections])
assert DomainMatrix.from_Matrix(coords).convert_to(s.QQ).rank()==55
assert all(not any(jac(P+v)[k]-jac(v)[k] for k in triples) for v in Bdirections)

# These two quadrics on the 30 frame variables are coprime.
obsrows=(rows1,rows2)
qdiagonal=[[obstruction(v,row)for v in Bdirections]for row in obsrows]
assert all(all(v==0 for v in row) for row in qdiagonal)
nonzero=[]
for i in range(30):
    for j in range(i+1,30):
        vals=tuple(obstruction(Bdirections[i]+Bdirections[j],rr) for rr in obsrows)
        if any(vals):nonzero.append((i,j,vals))
assert nonzero==[(21,29,(s.Integer(9),s.Integer(0))),(25,27,(s.Integer(0),s.Integer(-4)))]

# Independent dimension-28 logarithmic source: (P^4)^3 x P^14 x P^2.
qg=(u0**2+u1**2+z**2+w0**2+w1**2+u0*u1+u0*w0+u0*w1
    +2*u1*w0+u1*w1+z*w0+2*z*w1+w0*w1)
assert (s.hessian(qg,E)/2).det()==-s.Rational(1,8)
fs=[u0,u1,z,qg]
r1=s.Matrix([1,-1,0,0]);r2=s.Matrix([1,0,-1,0]);r3=s.Matrix([2,0,0,-1])
residues=[{(i,j):v[i]*w[j]-v[j]*w[i]for i,j in it.combinations(range(4),2)}for v,w in [(r1,r2),(r1,r3),(r2,r3)]]
c={ij:residues[0][ij]+2*residues[1][ij]+5*residues[2][ij]for ij in residues[0]}

def omega_source(f,c):
    out={ij:s.S.Zero for ij in pairs}
    for i,j in it.combinations(range(4),2):
        fac=s.prod(f[k]for k in range(4)if k not in (i,j))
        for h,l in pairs:
            out[h,l]+=c[i,j]*fac*(s.diff(f[i],x[h])*s.diff(f[j],x[l])-s.diff(f[i],x[l])*s.diff(f[j],x[h]))
    return {ij:s.expand(v)for ij,v in out.items()}

def source_derivative(f,c,m,d):
    out={ij:s.S.Zero for ij in pairs}
    for i,j in it.combinations(range(4),2):
        oth=[k for k in range(4)if k not in (i,j)]
        fac=s.prod(f[k]for k in oth)
        if m in oth:
            tail=s.prod(f[k]for k in oth if k!=m)
            for h,l in pairs:
                out[h,l]+=c[i,j]*tail*d*(s.diff(f[i],x[h])*s.diff(f[j],x[l])-s.diff(f[i],x[l])*s.diff(f[j],x[h]))
        if m in (i,j):
            for h,l in pairs:
                dih=s.diff(d,x[h])if m==i else s.diff(f[i],x[h])
                dil=s.diff(d,x[l])if m==i else s.diff(f[i],x[l])
                djh=s.diff(d,x[h])if m==j else s.diff(f[j],x[h])
                djl=s.diff(d,x[l])if m==j else s.diff(f[j],x[l])
                out[h,l]+=c[i,j]*fac*(dih*djl-dil*djh)
    return {ij:s.expand(v)for ij,v in out.items()}

variations=[]
for j in range(3):
    variations.extend(source_derivative(fs,c,j,x[i])for i in range(5)if i!=j)
for i,j in mon2:
    if (i,j)!=(4,4):variations.append(source_derivative(fs,c,3,x[i]*x[j]))
variations.extend(omega_source(fs,cd)for cd in residues[1:])
assert len(variations)==28
source_matrix=s.Matrix([[s.Poly(D[I],*x).coeff_monomial(s.prod(x[k]for k in m))for D in variations]for I in pairs for m in mon3])
minor_rows=[2,3,4,6,7,8,9,10,11,16,19,20,21,25,26,27,28,29,30,38,39,42,43,45,46,52,53,55]
det=source_matrix[minor_rows,:].det()
assert det==-57177414000000
assert DomainMatrix.from_Matrix(source_matrix).convert_to(s.QQ).rank()==28
base_omega=omega_source(fs,c)
aug=source_matrix.row_join(s.Matrix([s.Poly(base_omega[I],*x).coeff_monomial(s.prod(x[k]for k in m)) for I in pairs for m in mon3]))
assert DomainMatrix.from_Matrix(aug).convert_to(s.QQ).rank()==29

result={
 'date':'2026-11-10','version':'v10','status':'PASS','sympy':s.__version__,
 'logarithmic_source_type':[1,1,1,2],'logarithmic_source_dimension':28,
 'logarithmic_sample_hessian_det':'-1/8','logarithmic_derivative_rank':28,
 'logarithmic_projective_augmented_rank':29,
 'logarithmic_minor_rows':minor_rows,'logarithmic_minor_det':str(det),
 'ambient_tangent_frame_rank_with_gauge_and_scale':55,
 'ambient_tangent_dimension_from_prior_linearization':30,
 'ambient_1st_witness':[(list(I),str(mon),str(w))for I,mon,w in rows1],
 'ambient_2nd_witness':[(list(I),str(mon),str(w))for I,mon,w in rows2],
 'tangent_quadrics':['9*t21*t29','-4*t25*t27'],
 'verified_tangent_frame_basis_length':30,'exact_local_dimension':28,
 'geometric_claim':'The point lies on the 28-dimensional rank-two log family and in the closure of a rank-four family; reduced components intersect. Ambient nonreducedness and rank-four component identity remain unresolved.'
}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print('PASS: exact local dimension 28; two coprime obstruction quadrics; rank-28 logarithmic source')
