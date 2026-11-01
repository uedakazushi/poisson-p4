#!/usr/bin/env python3
"""Exact one-parameter toric inverse algebra and geometric S3 certificate.

The accompanying proof uses geometric irreducibility of GL4/S5 and the
non-square discriminant on this line; it does not determine the full S61 fiber.
"""
import sys,json,time,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from toric_algebra_fiber import multiplication_from_jet, TRIPLES
from weight1_structure import constraints
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
started=time.monotonic();t=s.Symbol('t');K=QQ.frac_field(t)
data=json.loads((ROOT/'reports/TORIC_ALGEBRA_FIBER_20261031.json').read_text())
cert=json.loads((ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json').read_text())
diridx=0
v=s.Matrix([s.Rational(x) for x in data['fixed_first_jet']])+t*s.Matrix.vstack(*constraints()).nullspace()[diridx]
base=multiplication_from_jet(v,{a:s.S.Zero for a in TRIPLES})
zero=s.zeros(24,1)
kernels=[multiplication_from_jet(zero,{a:s.S.One if a==b else s.S.Zero for a in TRIPLES}) for b in TRIPLES]
R,*cs=ring(','.join('c'+str(i) for i in range(20)),K)
m=[[[R(K.from_sympy(base[i][j,k]))+sum((cs[z]*int(kernels[z][i][j,k]) for z in range(20)),R.zero) for k in range(4)] for j in range(4)] for i in range(4)]
def assoc(m,ring):
 q=[[sum((m[i][a][b]*m[j][b][a] for a in range(4) for b in range(4)),ring.zero)/3 for j in range(4)] for i in range(4)]
 out=[]
 for i in range(4):
  for j in range(i+1,4):
   for a in range(4):
    for b in range(4):
     out.append(sum((m[i][a][k]*m[j][k][b]-m[j][a][k]*m[i][k][b] for k in range(4)),ring.zero)+(q[j][b] if a==i else ring.zero)-(q[i][b] if a==j else ring.zero))
 return q,out
q,eq=assoc(m,R);print('equations',time.monotonic()-started,flush=True)
lin=[sum((eq[i]*K.from_sympy(s.Rational(c)) for i,c in a['constant_combination_of_96_equations']),R.zero) for a in cert['universal_linear_consequences']]
assert all(max(map(sum,e.keys()),default=0)<=1 for e in lin)
L=[[e.coeff(cs[i]) for i in range(18)] for e in lin]
B=[[-e.coeff(R.one),-e.coeff(cs[18]),-e.coeff(cs[19])] for e in lin]
DM=DomainMatrix(L,(18,18),K)
D=DM.det();print('linear determinant',K.to_sympy(D).factor(),flush=True)
sol=DM.inv().matmul(DomainMatrix(B,(18,3),K)).to_list();print('solved',time.monotonic()-started,flush=True)
S,x,y=ring('x,y',K)
ss=[sum((S(c)*p for c,p in zip(row,[S.one,x,y])),S.zero) for row in sol]+[x,y]
mm=[[[S(K.from_sympy(base[i][j,k]))+sum((ss[z]*int(kernels[z][i][j,k]) for z in range(20)),S.zero) for k in range(4)] for j in range(4)] for i in range(4)]
qq,rr=assoc(mm,S)
rows=cert['residual_leading_pivot_rows_zero_based']; mons=[x*x,x*y,y*y]
C=[[rr[i].coeff(z) for z in mons] for i in rows]
CMD=DomainMatrix(C,(3,3),K);eps=CMD.det(); inv=CMD.inv().to_list()
low=[rr[i]-sum((S(C[a][b])*mons[b] for b in range(3)),S.zero) for a,i in enumerate(rows)]
nf=[-sum((S(inv[i][j])*low[j] for j in range(3)),S.zero) for i in range(3)]
assert all(max(map(sum,f.keys()),default=0)<=1 for f in nf)
for f in rr: assert f-sum((S(f.coeff(mons[j]))*(mons[j]-nf[j]) for j in range(3)),S.zero)==0
vec=lambda f:[f.coeff(S.one),f.coeff(x),f.coeff(y)]
nx,nxy,ny=map(vec,nf)
MX=DomainMatrix([[K.zero,nx[0],nxy[0]],[K.one,nx[1],nxy[1]],[K.zero,nx[2],nxy[2]]],(3,3),K)
MY=DomainMatrix([[K.zero,nxy[0],ny[0]],[K.zero,nxy[1],ny[1]],[K.one,nxy[2],ny[2]]],(3,3),K)
assert MX.matmul(MY)==MY.matmul(MX)
T=s.Symbol('z');char=MY.charpoly(); cp=sum(K.to_sympy(a)*T**(3-i) for i,a in enumerate(char))
disc=s.factor(s.discriminant(cp,T));print('exact three-dimensional algebra computed',flush=True)
result={'direction':diridx,'beta':[str(a) for a in v],'linear_determinant':str(s.factor(K.to_sympy(D))),'leading_determinant':str(s.factor(K.to_sympy(eps))),'cubic':str(s.factor(cp)),'discriminant':str(disc),'nf':[str(f.as_expr()) for f in nf],'all96_residuals_in_3quadratic_span':True,'Mx':[[str(K.to_sympy(a)) for a in row] for row in MX.to_list()],'My':[[str(K.to_sympy(a)) for a in row] for row in MY.to_list()]}
identity=DomainMatrix.eye(3,K).to_dense()
for lhs,coeffs in [(MX.matmul(MX),nx),(MX.matmul(MY),nxy),(MY.matmul(MY),ny)]:
 assert lhs==identity.scalarmul(coeffs[0])+MX.scalarmul(coeffs[1])+MY.scalarmul(coeffs[2])
unit=DomainMatrix([[K.one],[K.zero],[K.zero]],(3,1),K)
cyclic=unit.hstack(MY.matmul(unit)).hstack(MY.matmul(MY).matmul(unit))
cyclic_det=cyclic.det();assert cyclic_det
num,den=s.fraction(s.cancel(disc))
num_factor=s.factor_list(num,t);den_factor=s.factor_list(den,t)
assert len(num_factor[1])==len(den_factor[1])==1
Q,nexp=num_factor[1][0];D,dexp=den_factor[1][0]
assert (s.degree(Q,t),nexp,s.degree(D,t),dexp)==(12,2,10,3)
u,w,g=s.gcdex(D,s.diff(D,t),t);assert g==1
u2,w2,g2=s.gcdex(D,Q,t);assert g2==1
assert s.expand(u*D+w*s.diff(D,t))==1
assert s.expand(u2*D+w2*Q)==1
assert s.factor(disc.subs(t,0))==s.Rational(cert['cubic_discriminant'])
result.update({'id':'TORIC_MONODROMY_20261101','all_passed':True,'field':'QQ(t)',
 'linear_solutions_c0_through_c19':[str(f.as_expr()) for f in ss],
 'cyclic_determinant':str(s.factor(K.to_sympy(cyclic_det))),
 'discriminant_numerator_factor':str(Q),'discriminant_denominator_factor':str(D),
 'discriminant_factor_exponents':[2,3],
 'squarefree_denominator_bezout':[str(u),str(w)],
 'coprime_num_den_bezout':[str(u2),str(w2)],
 'geometric_monodromy':'S3, using the geometric proof note and prior generic degree three theorem',
 'limitations':['No full S61 finiteness or exhaustion.','The line avoids finitely many localization zeros; the ten odd discriminant poles are not asserted to lie in the affine toric locus.','No full branch divisor or classification of all exceptional fibers.'],
 'proof_note':'research/notes/TORIC_MONODROMY_20261101.md',
 'input_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'code/toric_algebra_fiber.py',ROOT/'code/weight1_structure.py',ROOT/'reports/TORIC_ALGEBRA_FIBER_20261031.json',ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json']}})
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps({'all_passed':True,'non_square_over_Ct':True,'generic_geometric_monodromy':'S3','output':str(args.output)}))
