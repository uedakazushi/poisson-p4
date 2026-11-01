#!/usr/bin/env python3
"""Normal cone and complete coordinate sections at the zero-trace AB orbit."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s
from itertools import combinations,combinations_with_replacement
from toric_zero_jet_dimension_20261101 import witness,encode
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
if args.output.exists():raise FileExistsError(args.output)
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
b0,b1=s.symbols('b0 b1');B=s.Matrix([[b0,b1],[b1,t]]);Bp=s.Matrix([[b0-b1*b1/t,b1/t],[b1/t,-1/t]])
rotation=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,-1,0]])
h=((rotation*upper(Bp)).inv()*upper(B)).applyfunc(s.cancel);h0=h.subs({b0:0,b1:0});K=h.diff(b1).subs({b0:0,b1:0})*h0.inv()
baseQ=Qof(s.zeros(3),x*x);SB=upper(s.Matrix([[0,0],[0,t]]));base0=push(SB.inv(),baseQ);baseinf=push(h0,base0)
N0,g0=extract(base0);Ni,gi=extract(baseinf)
expected=s.zeros(3);expected[0,2]=-t/12
assert N0==expected and g0==s.Matrix([1,0,0])
expected[0,2]=-s.Rational(1,12)
assert Ni==expected and gi==s.zeros(3,1)
# finite-chart normal free coordinates b0,b1,n00,n01,n02,n12,g0-1,g1,g2.
core=[]
for i,j in ((0,0),(0,1),(0,2),(1,2)):
 NN=s.zeros(3);NN[i,j]=1
 if (i,j)==(0,0):NN[2,2]=-1
 core.append(Qof(NN))
core +=[Qof(s.zeros(3),g) for g in (x*x,x*z,z*z)]
Ni,gi=extract(baseinf);nc=Ni[0,2]
# infinity smooth rankone chart: N=(1,p,q)^t(-p*r-q*s,r,s), coordinates p,q,r,s.
def ci(Q):
 NN,gg=extract(Q)
 return s.Matrix([NN[1,2]/nc,NN[2,2]/nc,NN[0,1],NN[0,2],*gg])
columns=[]
for k in range(9):
 out=s.zeros(9,1)
 if k==0:out[0]=1
 elif k==1:
  out[1]=1/t
  D=K*baseinf+baseinf*K.T-sum(((K*Y)[i]*baseinf.diff(Y[i]) for i in range(4)),s.zeros(4))
  out[2:,0]=ci(D)
 else:out[2:,0]=ci(push(h0,core[k-2]))
 columns.append(out)
T=s.Matrix.hstack(*columns).applyfunc(s.factor)
assert s.factor(T.det())==144*t**5
# global coordinate derivatives including all13 normal coordinates at vertex.
allcore=[]
for i,j in ((0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1)):
 NN=s.zeros(3);NN[i,j]=1
 if i==j:NN[2,2]=-1
 allcore.append(Qof(NN))
allcore +=[Qof(s.zeros(3),g) for g in (x*x,x*z,z*z)]
pairs=list(combinations(range(4),2));mon=[Y[i]*Y[j] for i,j in combinations_with_replacement(range(4),2)]
def vec(Q):return s.Matrix([s.Poly(s.expand(Q[i,j]),*Y).coeff_monomial(m) for i,j in pairs for m in mon])
Ds=[]
for i,j in ((0,0),(0,1)):
 A=s.zeros(4);A[(0,2)[i],(1,3)[j]]=A[(0,2)[j],(1,3)[i]]=1
 Ds.append(vec(A*baseQ+baseQ*A.T-sum(((A*Y)[k]*baseQ.diff(Y[k]) for k in range(4)),s.zeros(4))))
Ds +=[vec(push(SB,C)) for C in allcore]
D=s.Matrix.hstack(*Ds)
# torsion relations finite: t*n10=t*n11=t*n20=t*n21=0.
torsion=[5,6,8,9]
for k in torsion:D[:,k]=D[:,k].subs(t,0)
cs=[];labels=[]
for j in range(D.cols):
 for n in range(5):
  col=D[:,j].applyfunc(lambda f:s.Poly(f,t).coeff_monomial(t**n))
  if any(col):cs.append(col);labels.append([j,n])
W=s.Matrix.hstack(*cs)
coordinate_witness=witness(W)
assert coordinate_witness['rank']==18

# The free quotient of the conormal sheaf has rank9 and degree5. Its sections
# have polynomial degree at most3 on this chart, since entries of T do.
free_columns=[0,1,2,3,4,7,10,11,12]
Df=D[:,free_columns]
Di=(Df*T.inv()).applyfunc(s.cancel)
assert all(s.limit(v,t,s.oo).is_finite for v in Di)
infinity_D=Di.applyfunc(lambda v:s.limit(v,t,s.oo))
free_zero_witness=witness(Df.subs(t,0))
free_infinity_witness=witness(infinity_D)
full_zero_witness=witness(D.subs(t,0))
assert free_zero_witness['rank']==free_infinity_witness['rank']==9
assert full_zero_witness['rank']==13
unknown=s.symbols('v0:36')
section=s.Matrix([unknown[4*j:4*j+4] for j in range(9)])*s.Matrix([1,t,t*t,t**3])
at_infinity=(section.T*T.inv()).applyfunc(s.expand)
equations=[]
for v in at_infinity:
    assert all(-4<=term.as_coeff_exponent(t)[1]<=7 for term in s.Add.make_args(v))
    equations.extend(v.coeff(t,k) for k in range(1,8))
section_matrix=s.linear_eq_to_matrix(equations,unknown)[0]
assert 36-section_matrix.rank()==14

# A strict symplectic shear fixes the tensor and is transitive on the
# exceptional P1 away from its unique singular-source point t=0.
u=s.Symbol('u')
shear=s.eye(4);shear[3,2]=u
J4=s.diag(J2,J2)
assert shear*J4*shear.T==J4 and push(shear,baseQ)==baseQ
L=s.Matrix([[0,0],[1,0],[0,t],[0,1]])
assert shear*L==L.subs(t,t/(1+u*t))*s.diag(1,1+u*t)

# Exact Rees calculation for the graph line N=-(t/12)E13 in the minimal
# nilpotent sl3 cone. Rescaling t by -1/12 gives the following universal chart.
aa,bb,cc,dd,ee,ff,gg,hh,eps=s.symbols('A B C D E F G H eps')
variables=(t,aa,bb,cc,dd,ee,ff,gg,hh,eps)
K=[ee*gg-dd*hh,bb*gg-aa*hh,aa*ee+ee**2+ff*hh,
   bb*dd+ee**2+ff*hh,aa*dd+dd*ee+ff*gg]
lifted=K+[t*dd-eps*(aa*ff-cc*dd),t*ee-eps*(bb*ff-cc*ee),
          t*gg+eps*(aa**2+aa*ee+cc*gg),t*hh+eps*(aa*bb+bb*ee+cc*hh)]
gb=s.groebner(lifted,*variables,order='lex',domain=s.QQ)
assert len(gb.polys)==10
assert all(p.LM(order='lex').exponents[-1]==0 for p in gb.polys)
# Consequently multiplication by eps is injective: this is already the
# eps-saturated Rees ideal. Localized at eps it is the original rank-one ideal.
NM=s.Matrix([[eps*aa,eps*bb,t+eps*cc],[eps*dd,eps*ee,eps*ff],
             [eps*gg,eps*hh,-eps*(aa+ee)]])
minors=[s.expand(NM.extract(rows,cols).det())
        for rows in combinations(range(3),2) for cols in combinations(range(3),2)]
assert all(gb.reduce(v)[1]==0 for v in minors)
divided=[s.cancel(v/eps**2) if not s.expand(v).coeff(eps,1) else s.cancel(v/eps)
         for v in minors]
gb_divided=s.groebner(divided,*variables,order='lex',domain=s.QQ)
assert all(gb_divided.reduce(v)[1]==0 for v in lifted)
initial=K+[t*dd,t*ee,t*gg,t*hh]
gb_initial=s.groebner(initial,*variables[:-1],order='lex',domain=s.QQ)
assert all(gb_initial.reduce(p.as_expr().subs(eps,0))[1]==0 for p in gb.polys)
assert all(s.groebner([p.as_expr().subs(eps,0) for p in gb.polys],
                      *variables[:-1],order='lex',domain=s.QQ).reduce(v)[1]==0 for v in initial)
# The torsion ideal (D,E,G,H) is killed by t and is generated in degree1.
assert all(gb_initial.reduce(t*v)[1]==0 for v in (dd,ee,gg,hh))
assert all(v.subs({dd:0,ee:0,gg:0,hh:0})==0 for v in K)

# Hilbert series of the vertical normal-cone component. Inclusion-exclusion
# on the six leading monomials avoids a finite-degree Hilbert-function guess.
gb_K=s.groebner(K,aa,bb,cc,dd,ee,ff,gg,hh,order='lex',domain=s.QQ)
leading=[p.LM(order='lex').exponents for p in gb_K.polys]
assert len(leading)==6
Z=s.Symbol('Z')
numerator=s.S.Zero
for mask in range(1<<len(leading)):
    chosen=[leading[i] for i in range(len(leading)) if mask&(1<<i)]
    degree=sum(max(m[j] for m in chosen) for j in range(8)) if chosen else 0
    numerator+=(-1)**len(chosen)*Z**degree
assert s.expand(numerator-(1-Z)**3*(1+3*Z+Z*Z))==0
free_hilbert=(1+4*Z)/(1-Z)**10
torsion_hilbert=(1+3*Z+Z*Z)/(1-Z)**10-1/(1-Z)**9
assert s.cancel(free_hilbert+torsion_hilbert-(1+8*Z+Z*Z)/(1-Z)**10)==0

paths=[Path(__file__),ROOT/'research/notes/AB_INCIDENCE_20261101.md',
       ROOT/'research/notes/AB_ONE_VARIABLE_LOCAL_20261101.md']
report=dict(id='AB_ZERO_TRACE_LOCAL_20261101',all_passed=True,field='QQ',
    base_Q='x^2 partial_y wedge partial_w',
    finite_chart_core_N=encode(N0),infinity_chart_core_N=encode(Ni),
    free_normal_transition=encode(T),normal_transition_determinant=str(s.factor(T.det())),
    free_conormal_rank=9,free_conormal_degree=5,free_conormal_h0=14,
    torsion_conormal_length=4,full_conormal_h0=18,
    full_coordinate_differential=encode(D),coefficient_span_labels=labels,
    coordinate_span_witness=coordinate_witness,
    full_zero_witness=full_zero_witness,free_zero_witness=free_zero_witness,
    free_infinity_witness=free_infinity_witness,
    free_global_section_constraint_matrix=encode(section_matrix),
    Rees_groebner_basis=[str(p.as_expr()) for p in gb.polys],
    Rees_leading_exponents=[list(p.LM(order='lex').exponents) for p in gb.polys],
    normal_cone_core_ideal=list(map(str,initial)),
    vertical_component_groebner_basis=[str(p.as_expr()) for p in gb_K.polys],
    vertical_hilbert_numerator=str(s.expand(numerator)),
    tangent_cone_hilbert_series='(1+8Z+Z^2)/(1-Z)^10',embedding_dimension=18,multiplicity=10,
    conclusion='AB is normal at every nonzero point of its zero-trace one-variable locus; combined with prior charts, AB is normal away from the origin.',
    limitations=['The origin is not treated.','No assertion that the completed local ring equals its tangent cone.',
                 'Uses the generalized completion-generation lemma and the all-degree sheaf argument in the proof note.'],
    proof_note='research/notes/AB_ZERO_TRACE_LOCAL_20261101.md',
    input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
print(json.dumps(dict(all_passed=True,embedding_dimension=18,multiplicity=10,output=str(args.output))))
