#!/usr/bin/env python3
"""Degree-two graded normalization check for AB; not all-degree normality."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
if args.output.exists():raise FileExistsError(args.output)
import sympy as s
from itertools import combinations,combinations_with_replacement
from toric_zero_jet_dimension_20261101 import encode
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

# The three Levi pieces of the degree-two fiber coordinate algebra. Weights
# below are Schur-functor weights in the tautological-dual convention.
core_weights=[(2,0),(1,1),(0,2),(3,-1),(2,0),(1,1),(0,2),(-1,3)]
vertical_weights=[(3,1),(2,2),(1,3)]
pairs=list(combinations_with_replacement(range(8),2))
monomials=[pars[i]*pars[j] for i,j in pairs]
pair_weights=[tuple(core_weights[i][r]+core_weights[j][r] for r in range(2)) for i,j in pairs]
core2=Counter(pair_weights)
relations=[s.Poly(s.expand(N.extract(rr,cc).det()),*pars)
           for rr in combinations(range(3),2) for cc in combinations(range(3),2)]
relation_ranks={}
for wt in sorted(core2):
    indices=[i for i,v in enumerate(pair_weights) if v==wt]
    mat=s.Matrix([[f.coeff_monomial(monomials[i]) for i in indices] for f in relations])
    relation_ranks[wt]=mat.rank();core2[wt]-=mat.rank()
assert sum(relation_ranks.values())==9 and sum(core2.values())==27
cross=Counter(tuple(a[r]+b[r] for r in range(2)) for a in core_weights for b in vertical_weights)
vertical2=Counter(tuple(vertical_weights[i][r]+vertical_weights[j][r] for r in range(2))
                  for i,j in combinations_with_replacement(range(3),2))
def decompose(character):
    c=Counter({a:n for a,n in character.items() if n});answer=[]
    while c:
        a,b=max(c);multiplicity=c[(a,b)]
        assert a>=b and multiplicity>0
        answer.append((a,b,multiplicity))
        for i in range(a-b+1):
            key=(a-i,b+i);c[key]-=multiplicity
            assert c[key]>=0
            if c[key]==0:del c[key]
    return answer
pieces=[decompose(c) for c in (core2,cross,vertical2)]
assert pieces==[[(6,-2,1),(5,-1,1),(4,0,2),(2,2,1)],
               [(6,0,1),(5,1,2),(4,2,2),(3,3,1)],[(6,2,1),(4,4,1)]]
def bott(a,b):
    shifted=(a+2,b+1)
    if 0 in shifted or abs(shifted[0])==abs(shifted[1]):return None
    ordered=sorted(map(abs,shifted),reverse=True)
    partition=(ordered[0]-2,ordered[1]-1)
    length=sum(v<0 for v in shifted)+(sum(shifted)<0)
    A,B=partition[0]-partition[1],partition[1]
    dimension=(A+1)*(B+1)*(A+B+2)*(A+2*B+3)//6
    return dict(degree=length,Sp4_partition=list(partition),dimension=dimension)
bott_table=[];h0=0;h1=0
for j,piece in enumerate(pieces):
    for a,b,multiplicity in piece:
        result=bott(a,b)
        bott_table.append(dict(piece=j,weight=[a,b],multiplicity=multiplicity,cohomology=result))
        if result:
            assert result['degree'] in (0,1)
            if result['degree']==0:h0+=multiplicity*result['dimension']
            else:h1+=multiplicity*result['dimension']
assert (h0,h1)==(845,84)

# Nonzero connecting map: a P-map Sym6(V)->E2 would send y^6 to a scalar
# multiple of h0*g0. The indicated unipotent kills y^6 but not h0*g0.
g0,g1,g2,un=s.symbols('g0 g1 g2 un')
generalQ=Qof(N,g0*x*x+g1*x*z+g2*z*z)
shear=s.eye(4);shear[3,2]=un
newN,newg=extract(push(shear,generalQ))
assert newN==N
assert s.expand(newg[0]-g0)==12*un*hh[0]
assert s.diff(y**6,w)==0
candidate_weights=Counter(tuple(-a[r]-b[r] for r in range(2)) for a in core_weights for b in vertical_weights)
assert candidate_weights[(-6,0)]==1
point={v:0 for v in pars};point[hh[0]]=1
assert N.subs(point).rank()==1 and s.expand(12*hh[0]**2).subs(point)==12

# Dominant incidence chart and an exact rational-rank lower bound through
# a nonzero minor modulo32003. A finite-field upper bound is not used.
a,b,d,e,b0,b1,b2=s.symbols('a b d e b0 b1 b2')
chartN=s.Matrix([1,a,b])*s.Matrix([[-a*d-b*e,d,e]])
chartQ=push(upper(s.Matrix([[b0,b1],[b1,b2]])),Qof(chartN,g0*x*x+g1*x*z+g2*z*z))
variables=(a,b,d,e,b0,b1,b2,g0,g1,g2)
R,*unused=ring(','.join(map(str,variables)),QQ)
spatial_monomials=[Y[i]*Y[j] for i,j in combinations_with_replacement(range(4),2)]
polynomials=[R.from_expr(s.Poly(chartQ[i,j],*Y).coeff_monomial(mon))
             for i,j in combinations(range(4),2) for mon in spatial_monomials]
prime=32003
def row_reduce(poly,cache):
    assert all(int(c.denominator)%prime for c in poly.values())
    row={m:int(c.numerator)%prime*pow(int(c.denominator),-1,prime)%prime for m,c in poly.items() if c}
    while row:
        pivot=min(row)
        if pivot not in cache:
            leading=row[pivot];inverse=pow(leading,-1,prime)
            cache[pivot]={m:v*inverse%prime for m,v in row.items()}
            return pivot,leading
        factor=row[pivot]
        for mon,coefficient in cache[pivot].items():
            value=(row.get(mon,0)-factor*coefficient)%prime
            if value:row[mon]=value
            else:row.pop(mon,None)
    return None
basis={};selected=[]
for i,f in enumerate(polynomials):
    if row_reduce(f,basis):selected.append(i)
assert len(selected)==45
quadratic_basis={};pivots=[];selected_products=[];pivot_product=1
for i,j in combinations_with_replacement(range(45),2):
    result=row_reduce(polynomials[selected[i]]*polynomials[selected[j]],quadratic_basis)
    if result:
        mon,coefficient=result;pivots.append(mon);selected_products.append((i,j))
        pivot_product=pivot_product*coefficient%prime
assert len(quadratic_basis)==761 and pivot_product!=0
paths=[Path(__file__),ROOT/'research/notes/AB_INCIDENCE_20261101.md',
       ROOT/'research/notes/AB_ZERO_TRACE_LOCAL_20261101.md']
report=dict(id='AB_ORIGIN_DEGREE_TWO_20261101',all_passed=True,
    Levi_degree_two_decomposition=pieces,Borel_Weil_Bott=bott_table,
    associated_graded_h0=845,associated_graded_h1=84,
    unipotent_identity='d/du (h0*g0)=12*h0^2; d/du(y^6)=0 for w->w+u*z',
    normalization_degree_two_dimension=761,normalization_degree_two_H1=0,
    prime=prime,coordinate_indices=selected,product_indices=selected_products,
    pivot_monomials=pivots,triangular_pivot_product_mod_prime=pivot_product,
    rational_rank_argument='The nonzero modular minor gives QQ rank>=761; the geometric normalization H0 upper bound761 proves equality.',
    conclusion='All degree-two normalization functions are products of degree-one coefficient functions; no extra normalization generator occurs in degree2.',
    limitations=['Does not prove all-degree generation or normality at the origin.',
                 'Uses the parabolic extension and Frobenius reciprocity argument in the proof note, not just the associated-graded cohomology dimensions.'],
    proof_note='research/notes/AB_ORIGIN_DEGREE_TWO_20261101.md',
    input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
print(json.dumps(dict(all_passed=True,degree_two_dimension=761,modular_minor=pivot_product)))
