#!/usr/bin/env python3
"""Generic toric first-jet cover equals nonzero Jacobian two-torsion.

Exact calculations take place on a certified dominant six-dimensional slice.
Equivariance and uniqueness of the rational C9 recovery extend the identities
universally. This theorem concerns the toric component, not full S61 exhaustion.
"""
import argparse,hashlib,json,math,time
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from toric_algebra_fiber import TRIPLES,multiplication_from_jet,associativity_quadrics
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
start=time.time();slice_path=ROOT/'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json';u_path=ROOT/'reports/TORIC_DEGREE5_COVARIANT_20261101.json';d_path=ROOT/'reports/TORIC_DEGREE3_SCALAR_20261101.json'
data=json.loads(slice_path.read_text());model=json.loads(u_path.read_text());ddata=json.loads(d_path.read_text())
assert data['dominance_tangent_determinant']=='8428695240'
v=s.symbols('s0:6');c=s.symbols('c0:20');R,*sv=ring('s0,s1,s2,s3,s4,s5',QQ)
def poly(items):return R.from_dict({tuple(mon):QQ.convert(s.Rational(co)) for mon,co in items})
U=[poly(p) for p in data['C5_sparse']];c4=poly(data['c4_sparse']);c6=poly(data['c6_sparse'])
bs=[sum((QQ.convert(s.Rational(co))*sv[j] for j,co in enumerate(row)),R.zero) for row in data['slice_matrix_b']]
D=[]
for terms in ddata['D3_sparse']:
    total=R.zero
    for mon,co in terms:
        value=R(QQ.convert(s.Rational(co)))
        for i,n in enumerate(mon):
            if n:value*=bs[i]**n
        total+=value
    D.append(total)
def pairing(left,right):
    total=R.zero
    for i,triple in enumerate(TRIPLES):
        complement=tuple(sorted(j^1 for j in triple));sign=(-1)**sum(j%2 for j in triple)
        multiplicity=math.factorial(3)//math.prod(math.factorial(triple.count(j)) for j in range(4))
        total+=sign*multiplicity*left[i]*right[TRIPLES.index(complement)]
    return total
assert pairing(D,U)==54*c4
beta=s.Matrix(model['weight_basis_24_by16'])*s.Matrix(data['slice_matrix_b'])*s.Matrix(v)
m=multiplication_from_jet(beta,dict(zip(TRIPLES,c)));_,E=associativity_quadrics(m);C=s.Matrix([f.subs(dict.fromkeys(c,0)) for f in E]);B=s.Matrix(E).jacobian(c).subs(dict.fromkeys(c,0));H=[s.Poly(s.expand(f-sum(B[i,j]*c[j] for j in range(20))-C[i]),*c) for i,f in enumerate(E)]
Hterms=[]
for f in H:
 row=[]
 for mon,co in f.terms():
  inds=[i for i,n in enumerate(mon) for _ in range(n)];assert len(inds)==2;row.append((inds[0],inds[1],QQ.convert(co)))
 Hterms.append(row)
HH=s.Matrix([[f.coeff_monomial(c[i]*c[j]) for i in range(20) for j in range(i,20)] for f in H]);hrows=list(HH.T.rref()[1]);assert len(hrows)==45
Cpoly=[R.from_expr(f) for f in C];Bpoly=[[R.from_expr(f) for f in row] for row in B.tolist()]
@lru_cache(maxsize=210)
def UU(i,j):return U[i]*U[j]
F=[]
for i,row in enumerate(Hterms):
 h=sum((co*UU(a,b) for a,b,co in row),R.zero);F.append(-h/QQ(3)-QQ(27)*c4*Cpoly[i])
print('F built',time.time()-start,flush=True);UU.cache_clear()
Bcoeff=[B.diff(z).applyfunc(s.Rational) for z in v];B0=Bcoeff[0];rows=list(B0.T.rref()[1]);assert len(rows)==20;inv=B0[rows,:].inv();As=[inv*mat[rows,:] for mat in Bcoeff[1:]]
def frac(x):return Fraction(int(s.numer(x)),int(s.denom(x)))
invrows=[[(j,frac(co)) for j,co in enumerate(row) if co] for row in inv.tolist()];arows=[[[(j,frac(co)) for j,co in enumerate(row) if co] for row in a.tolist()] for a in As]
Fco=[]
for i in rows:
 dic={}
 for mon,co in F[i].items():
  assert sum(mon)==10;dic[mon[1:]]=Fraction(int(co.numerator),int(co.denominator))
 Fco.append(dic)
from itertools import combinations_with_replacement
monomials=[]
for degree in range(10):
 for indices in combinations_with_replacement(range(5),degree):monomials.append(tuple(indices.count(i) for i in range(5)))
values={};degree=-1
for mon in monomials:
 if sum(mon)!=degree:degree=sum(mon);print('recursion',degree,time.time()-start,flush=True)
 f=[d.get(mon,Fraction(0)) for d in Fco];val=[sum((co*f[j] for j,co in row),Fraction(0)) for row in invrows]
 for k in range(5):
  if not mon[k]:continue
  old=list(mon);old[k]-=1;previous=values[tuple(old)]
  for i,row in enumerate(arows[k]):val[i]-=sum((co*previous[j] for j,co in row),Fraction(0))
 values[mon]=val
V=[]
for i in range(20):V.append(R.from_dict({(9-sum(mon),*mon):QQ(co.numerator,co.denominator) for mon,vec in values.items() if (co:=vec[i])}))
print('V built terms',sum(len(f) for f in V),time.time()-start,flush=True)
for i in range(96):assert sum((Bpoly[i][j]*V[j] for j in range(20)),R.zero)==F[i],i
print('FIRST identity all96 PASS',time.time()-start,flush=True)
def enc(f):return [[list(mon),str(co)] for mon,co in sorted(f.items())]
assert pairing(D,V)==162*c6

BU=[sum((Bpoly[i][j]*U[j] for j in range(20)),R.zero) for i in range(96)];BV=F
# The51 complementary constant relations annihilate C,BU,BV. Therefore both
# remaining residue vectors lie in the45-dimensional columnspace of H.
hcols=list(HH[hrows,:].rref()[1]);lift=HH[:,hcols]*HH[hrows,hcols].inv();assert lift*HH[hrows,:]==HH
lift_rows=[[(j,QQ.convert(co)) for j,co in enumerate(row) if co] for row in lift.tolist()]
for vector in (Cpoly,BU,BV):
    for i,row in enumerate(lift_rows):
        assert vector[i]==sum((co*vector[hrows[j]] for j,co in row),R.zero),(i,'columnspace')
print('ALL51 columnspace relations forC,BU,BV PASS',time.time()-start,flush=True)
for count,i in enumerate(hrows):
 h=sum((co*(U[a]*V[b]+V[a]*U[b]) for a,b,co in Hterms[i]),R.zero)
 assert h+QQ(54)*c4*BU[i]+QQ(486)*c6*Cpoly[i]==R.zero,('second',i)
 if count%10==0:print('second',count,time.time()-start,flush=True)
print('SECOND45 PASS',time.time()-start,flush=True)
for count,i in enumerate(hrows):
 h=sum((co*V[a]*V[b] for a,b,co in Hterms[i]),R.zero)
 assert h+QQ(162)*c6*BU[i]-QQ(27)*c4*BV[i]+QQ(729)*c4*c4*Cpoly[i]==R.zero,('third',i)
 if count%5==0:print('third',count,time.time()-start,flush=True)
print('THIRD45 PASS',time.time()-start,flush=True)

# Elementary inverse and trace identities for the monic Weierstrass cubic.
a,b,z=s.symbols('a b z');f=z**3-27*a*z-54*b;derivative=3*z*z-27*a
inverse=(a*z*z-3*b*z-18*a*a)/(486*(a**3-b*b))
assert s.cancel(s.rem(s.together(derivative*inverse-1),f,z))==0
assert s.rem(z*derivative-54*a*z-162*b,f,z)==0
assert s.expand(a*(54*a)-18*a*a*3)==0
assert s.expand(a*(162*b)-3*b*(54*a))==0
report={'id':'TORIC_JACOBIAN_TWO_TORSION_20261101','all_passed':True,'field':'QQ[s0,...,s5] with certified dominantSp4-saturation',
 'variables':list(map(str,v)),'slice_matrix_b':data['slice_matrix_b'],'C9_sparse':list(map(enc,V)),'C9_total_monomials':sum(len(f) for f in V),
 'B0_pivot_rows':rows,'B0_minor_determinant':str(B0[rows,:].det()),'C9_recovery':'Unique rational solution of B(V)=-H(U,U)/3-27*c4*C; polynomialdegree9onthisslice.',
 'homogeneous_H_rank':45,'H_independent_rows':hrows,'H_independent_columns':hcols,'H_columnspace_lift':[[str(co) for co in row] for row in lift.tolist()],
 'all51_complementary_relations_annihilate_C_BU_BV':True,
 'first_residue_identity_equations_checked':96,'second_residue_identity_independent_rows_checked':45,'third_residue_identity_independent_rows_checked':45,
 'second_and_third_full96_follow_from_columnspace_relations':True,
 'scalar_pairings':{'D_U_equals_54_c4':True,'D_V_equals_162_c6':True},
 'generic_cubic':'X^3-27*c4(beta)*X-54*c6(beta)',
 'generic_isomorphism':'c_i=(C5_i(beta)*X+C9_i(beta))/(3*X^2-27*c4(beta))',
 'generic_inverse':'X=<D3(beta),c>',
 'inverse_derivative_mod_cubic':str(inverse),'universal_trace_c_zero':True,
 'universal_toric_S61_trace':'sum_of_three_toric_q = 2*Q51+Q52',
 'theorem':'The generic toric first-jet three-cover is canonically isomorphic to Jac(C_beta)[2] minus the zero section, with the displayed inverse maps.',
 'proof_dependencies':['Dominant six-dimensional slice and Fisher invariants.','UniversalSp4covariance ofC5andD3.','Genericrank20 ofB, ensuringuniqueequivariant rationalC9.','Establishedentire96equationgenericquotientdimension3.'],
 'limitations':['No assertion that the full S61 genericfiber has onlyfivepoints.','No assertion that rational C9 extends polynomially over allW16.','No classification of special fibers alongDelta=0.'],
 'proof_note':'research/notes/TORIC_JACOBIAN_TWO_TORSION_20261101.md'}
sources=[Path(__file__),ROOT/'code/toric_algebra_fiber.py',slice_path,u_path,d_path,ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json',ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json']
report['source_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
print(json.dumps({'all_passed':True,'generic_toric_Jacobian_2torsion_isomorphism':True,'full_S61_exhaustion_proved':False,'seconds':round(time.time()-start,2),'output':str(args.output)}))
