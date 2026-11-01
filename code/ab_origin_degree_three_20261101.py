#!/usr/bin/env python3
"""Degree-three normalization functions of AB are generated in degree one."""
import argparse,hashlib,json,time
from collections import Counter
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
if args.output.exists():raise FileExistsError(args.output)
started=time.monotonic()
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

h=[Counter(),Counter()];dims=[0,0];table=[]
for j in range(4):
    cy=character(3-j,core_weights,leading);cv=character(j,vertical_weights);prod=Counter()
    assert sum(cy.values())==(4-j)**3
    for a,ca in cy.items():
        for b,cb in cv.items():prod[(a[0]+b[0],a[1]+b[1])]+=ca*cb
    for a,b,multiplicity in decompose(prod):
        cohomology=bott(a,b);table.append(dict(vertical_degree=j,Levi_weight=[a,b],multiplicity=multiplicity,cohomology=cohomology))
        if cohomology:
            degree,wt,dim=cohomology;assert degree in(0,1);h[degree][wt]+=multiplicity;dims[degree]+=multiplicity*dim
assert dims==[7285,594] and h[1]==Counter({(9,1):1,(8,0):1})
assert h[0][9,1]==1 and h[0][8,0]==2

# Evaluation of P-cyclic lowest-weight vectors bounds all possible H1 cancellation.
g0,g1,g2,un=s.symbols('g0 g1 g2 un');fiber_variables=(*pars,g0,g1,g2)
weights=[tuple(-v for v in wt) for wt in core_weights+vertical_weights]
weight_monomials={}
for wt in [(-8,0),(-9,-1)]:
    weight_monomials[wt]=[s.prod(fiber_variables[i] for i in inds) for inds in combinations_with_replacement(range(11),3)
        if tuple(sum(weights[i][r] for i in inds) for r in range(2))==wt]
assert set(weight_monomials[-8,0])=={hh[0]**2*g1,hh[0]*hh[1]*g0,hh[0]*ss[0]*g0}
assert weight_monomials[-9,-1]==[hh[0]*g0**2]
generalQ=Qof(N,g0*x*x+g1*x*z+g2*z*z)
radical=s.eye(4);radical[3,2]=un
newN,newg=extract(push(radical,generalQ));assert newN==N
assert s.expand(newg[0]-g0)==12*un*hh[0]
assert s.expand(newg[1]-g1-un*(6*hh[1]-ss[0]/2))==0
T=s.Matrix([[1,un],[0,1]]);levi=s.eye(4)
for i,ii in enumerate((0,2)):
    for j,jj in enumerate((0,2)):levi[ii,jj]=T[i,j]
for i,ii in enumerate((1,3)):
    for j,jj in enumerate((1,3)):levi[ii,jj]=T.inv().T[i,j]
LN,Lg=extract(push(levi,generalQ));lp=inverse*LN.reshape(9,1)
assert lp[0]==ss[0] and lp[3]==hh[0] and s.expand(lp[4]-hh[1])==-4*un*hh[0]
assert Lg[0]==g0 and s.expand(Lg[1]-g1)==-2*un*g0
ca,cb,cc=s.symbols('ca cb cc');candidate=ca*hh[0]**2*g1+cb*hh[0]*hh[1]*g0+cc*hh[0]*ss[0]*g0
radical_derivative=s.diff(candidate.subs({g0:newg[0],g1:newg[1]},simultaneous=True),un).subs(un,0).expand()
levi_derivative=s.diff(candidate.subs(dict(zip(fiber_variables,(*lp,*Lg))),simultaneous=True),un).subs(un,0).expand()
assert s.expand(levi_derivative+(2*ca+4*cb)*hh[0]**2*g0)==0
assert s.expand(radical_derivative-(6*ca+12*cb)*hh[0]**2*hh[1]-(-ca/2+12*cc)*hh[0]**2*ss[0])==0
assert s.diff((hh[0]*newg[0]**2),un,2)==288*hh[0]**3
point={v:0 for v in pars};point[hh[0]]=1;assert N.subs(point).rank()==1

# Check independence on the dense rank-one trace-zero chart, not the free ambient variables.
a,b,d,e,b0,b1,b2=s.symbols('a b d e b0 b1 b2');chartN=s.Matrix([1,a,b])*s.Matrix([[-a*d-b*e,d,e]])
vals=inverse*chartN.reshape(9,1);sub=dict(zip(pars,vals))
forms=[s.Poly(f.subs(sub,simultaneous=True),a,b,d,e,g0,g1,g2) for f in weight_monomials[-8,0]]
mons=sorted(set().union(*(set(f.monoms()) for f in forms)));matrix=s.Matrix([[f.coeff_monomial(mm) for f in forms] for mm in mons]);assert matrix.rank()==3
constraints=[levi_derivative.subs(sub,simultaneous=True),radical_derivative.subs(sub,simultaneous=True)]
coeff=[]
for f in constraints:
    coeff.extend(s.Poly(f,a,b,d,e,g0,g1,g2).coeffs())
C=s.Matrix(coeff).jacobian([ca,cb,cc]);assert C.rank()==2 and C*s.Matrix([1,-s.Rational(1,2),s.Rational(1,24)])==s.zeros(C.rows,1)

# Exact coefficient-polynomial minor: only a lower bound is extracted modulo32003.
chartQ=push(upper(s.Matrix([[b0,b1],[b1,b2]])),Qof(chartN,g0*x*x+g1*x*z+g2*z*z))
variables=(a,b,d,e,b0,b1,b2,g0,g1,g2);R,*unused=ring(','.join(map(str,variables)),QQ)
spatial_monomials=[Y[i]*Y[j] for i,j in combinations_with_replacement(range(4),2)]
polynomials=[R.from_expr(s.Poly(chartQ[i,j],*Y).coeff_monomial(mon)) for i,j in combinations(range(4),2) for mon in spatial_monomials]
prime=32003

def reduce(poly,cache):
    assert all(int(c.denominator)%prime for c in poly.values())
    row={m:int(c.numerator)%prime*pow(int(c.denominator),-1,prime)%prime for m,c in poly.items() if c}
    while row:
        pivot=min(row)
        if pivot not in cache:
            lead=row[pivot];iv=pow(lead,-1,prime);cache[pivot]={m:v*iv%prime for m,v in row.items()};return pivot,lead
        factor=row[pivot]
        for mon,c in cache[pivot].items():
            value=(row.get(mon,0)-factor*c)%prime
            if value:row[mon]=value
            else:row.pop(mon,None)
    return None

basis={};selected=[]
for i,f in enumerate(polynomials):
    if reduce(f,basis):selected.append(i)
assert len(selected)==45
cubic_basis={};pivots=[];products=[];pivot_product=1
for i,j,k in combinations_with_replacement(range(45),3):
    result=reduce(polynomials[selected[i]]*polynomials[selected[j]]*polynomials[selected[k]],cubic_basis)
    if result:
        mon,c=result;pivots.append(mon);products.append((i,j,k));pivot_product=pivot_product*c%prime
assert len(cubic_basis)==6691 and pivot_product!=0
paths=[Path(__file__),ROOT/'research/notes/AB_INCIDENCE_20261101.md',ROOT/'research/notes/AB_ZERO_TRACE_LOCAL_20261101.md',ROOT/'research/notes/AB_ORIGIN_DEGREE_TWO_20261101.md']
report=dict(id='AB_ORIGIN_DEGREE_THREE_20261101',all_passed=True,Borel_Weil_Bott=table,associated_graded_h0=7285,associated_graded_h1=594,H0_multiplicities={str(k):v for k,v in h[0].items()},H1_multiplicities={str(k):v for k,v in h[1].items()},Sym8_evaluation_upper_bound=1,partition91_evaluation_upper_bound=0,normalization_degree_three_dimension=6691,normalization_degree_three_H1=0,prime=prime,coordinate_indices=selected,product_indices=products,pivot_monomials=pivots,triangular_pivot_product_mod_prime=pivot_product,conclusion='All degree-three normalization functions are cubic products of degree-one coefficient functions.',limitations=['Not an all-degree normality theorem.','P-cyclicity and evaluation argument are essential geometric parts of the proof; modular rank is used only as a rational lower bound.'],proof_note='research/notes/AB_ORIGIN_DEGREE_THREE_20261101.md',seconds=time.monotonic()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps({k:report[k] for k in ['all_passed','normalization_degree_three_dimension','triangular_pivot_product_mod_prime','seconds']}))
