#!/usr/bin/env python3
"""Independent finite5x5 pushforward audit of the full CA self-flag3jet."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
PARENT='AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63'
STEM=PARENT+'_REVIEW'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
INDEX={(p,m):i for i,(p,m) in enumerate(itertools.product(PAIRS,MONS))}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def dm(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ)
def product(A,B):return (dm(A)*dm(B)).to_Matrix()
def same(A,B):assert A.shape==B.shape and (A-B).applyfunc(s.expand)==s.zeros(*A.shape)

def physical_operators(M,B=None,C=None,order=3):
 """All coefficients of g P(g^-1 X) g^T by direct four-factor products.

 g=I+aM+a²B+a³C. Expand its matrix inverse directly, without logarithms,
 exponentials, Lie operators, or any author helper.
 """
 B=s.zeros(5) if B is None else B;C=s.zeros(5) if C is None else C
 g=[s.eye(5),M,B,C][:order+1]
 inverse=[s.eye(5),-M,M*M-B,-C+M*B+B*M-M**3][:order+1]
 for degree in range(order+1):
  assert sum((g[i]*inverse[degree-i] for i in range(degree+1)),s.zeros(5))==(s.eye(5) if degree==0 else s.zeros(5))
 grows=[[(d,j,Q(str(A[i,j]))) for d,A in enumerate(g) for j in range(5) if A[i,j]] for i in range(5)]
 # Output tensor wedge uses columns of g, substitution uses rows of inverse.
 gcols=[[(d,j,Q(str(A[j,i]))) for d,A in enumerate(g) for j in range(5) if A[j,i]] for i in range(5)]
 invrows=[[(d,j,Q(str(A[i,j]))) for d,A in enumerate(inverse) for j in range(5) if A[i,j]] for i in range(5)]
 terms=[{} for _ in range(order+1)]
 for col,((i,j),(r,t)) in enumerate(itertools.product(PAIRS,MONS)):
  for da,u,ca in gcols[i]:
   for db,v,cb in gcols[j]:
    if u==v or da+db>order:continue
    pair=(min(u,v),max(u,v));sign=1 if u<v else -1
    for dc,l,cc in invrows[r]:
     for dd,m,cd in invrows[t]:
      degree=da+db+dc+dd
      if degree>order:continue
      row=INDEX[pair,tuple(sorted((l,m)))];key=(row,col)
      terms[degree][key]=terms[degree].get(key,Q(0))+sign*ca*cb*cc*cd
 return [[(r,c,q) for (r,c),q in T.items() if q] for T in terms]

def act(op,v):
 ans=[Q(0)]*150
 for r,c,a in op:
  if v[c]:ans[r]+=a*v[c]
 return s.Matrix([s.Rational(q.numerator,q.denominator) for q in ans])

def main(output):
 started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports'/PARENT;report=load(folder/'report.json');assert report['all_passed']
 for n,h in report['input_sha256'].items():assert sha(ROOT/n)==h,n
 for n,h in report['evidence_sha256'].items():assert sha(folder/n)==h,n
 setup=load(folder/'setup.json');cert=load(folder/'certificate.json.gz')
 source=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';data=load(source)
 accepted=ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json';assert load(accepted)['all_passed']
 h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w])
 def canonical(P):
  d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
  return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
 def raw(f=0,g=0,H=0):
  P=s.zeros(5);P[0,2]=f;P[3,4]=g*x*x;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
  return canonical(P-P.T)
 def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b]) for (i,j), (a,b) in itertools.product(PAIRS,MONS)])
 cubics=[s.prod((x,h,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
 C=s.Matrix.hstack(*[vec(raw(f=f)) for f in [x*x,h*x,h*h]],vec(raw(g=1)),*[vec(raw(H=H)) for H in cubics])
 same(C,s.Matrix(data['fixed_fiber_linear_basis150_by20']))
 center=s.Matrix(data['center_fixed20']);P0=C*center
 Pphysical=s.zeros(5)
 for j,p in [(0,x*x),(3,h*x),(4,h*h+z*z)]:Pphysical[2,j]=p;Pphysical[j,2]=-p
 same(P0,vec(Pphysical))
 free=data['affine_slice_free_indices'];assert free==[1,2,3,9,10,14,15,17,18,19]
 F=[C[:,j] for j in free];fullflag=[tuple(v) for v in data['flag_positions']]
 assert len(fullflag)==9 and len(set(fullflag))==9 and (1,0) in fullflag
 flag=[v for v in fullflag if v!=(1,0)];assert [list(v) for v in flag]==setup['second_and_third_order_flag8']
 E=[]
 for i,j in flag:
  B=s.zeros(5);B[i,j]=1;E.append(B)
 M=s.zeros(5);M[1,0]=1;M[0,3]=-2;same(M,s.Matrix(setup['tangent_M']))
 aa=s.symbols('f0:9');generic=s.zeros(5)
 for a,(i,j) in zip(aa,fullflag):generic[i,j]=a
 assert generic**4==s.zeros(5) and (s.eye(5)+generic).det()==1
 rows=list(dm(C.T).rref()[1]);left=C[rows,:].inv()
 project=lambda V:V-product(C,product(left,V[rows,:]))
 base=physical_operators(M);fflag=[]
 for a,b in fullflag:
  L=s.zeros(5);L[a,b]=1;fflag.append(project(act(physical_operators(L,order=1)[1],P0)))
 first=s.Matrix.hstack(*fflag);assert dm(first).rank()==8
 normalized=s.zeros(9,1);normalized[fullflag.index((1,0))]=1;normalized[fullflag.index((0,3))]=-2
 same(product(first,normalized),s.zeros(150,1))
 Bops=[physical_operators(M,B) for B in E]
 unitops=[physical_operators(B,order=1)[1] for B in E]
 D=project(s.Matrix.hstack(*[act(base[1],v) for v in F],*[act(op,P0) for op in unitops]))
 constant=project(act(base[2],P0));same(D,s.Matrix(setup['all150_outside_fiber_order2_linear_matrix']))
 same(constant,s.Matrix(setup['all150_order2_constant']))
 pc=setup['pivot_columns'];pr=setup['pivot_rows'];inverse=D[pr,pc].inv();assert len(pc)==len(pr)==15
 same(inverse,s.Matrix(setup['inverse_pivot']));same(product(D[:,pc],product(inverse,D[pr,:])),D)
 origin=s.zeros(18,1);kernel=s.zeros(18,3);fc=[i for i in range(18) if i not in pc];assert fc==[12,14,17]
 for j,i in enumerate(fc):kernel[i,j]=1
 sd=-product(inverse,constant[pr,:]);sk=-product(inverse,D[pr,fc])
 for j,i in enumerate(pc):origin[i]=sd[j];kernel[i,:]=sk[j,:]
 same(origin,s.Matrix(cert['order2_affine_origin']));same(kernel,s.Matrix(cert['order2_affine_embedding']))
 same(product(D,origin)+constant,s.zeros(150,1));same(product(D,kernel),s.zeros(150,3))
 labels=[(),*[(i,) for i in range(18)],*[(i,10+j) for i in range(10) for j in range(8)]]
 assert [list(v) for v in labels]==cert['order3_low_monomials']
 # Direct finite-matrix a³ coefficient, with B varied one basis unit at a time.
 cols=[act(base[3],P0),*[act(base[2],v) for v in F],
       *[act(ops[3],P0)-act(base[3],P0) for ops in Bops],
       *[act(op,F[i]) for i in range(10) for op in unitops]]
 R3=project(s.Matrix.hstack(*cols));same(R3,s.Matrix(cert['order3_raw_quotient_coefficient_matrix']))
 quotient=R3-product(D[:,pc],product(inverse,R3[pr,:]));same(quotient,s.Matrix(cert['order3_coker_coefficient_matrix_before_substitution']))
 v0,v1,v2=s.symbols('v0 v1 v2');vv=(v0,v1,v2);affine=origin+kernel*s.Matrix(vv)
 words=[m for degree in range(3) for m in itertools.combinations_with_replacement(range(3),degree)]
 assert [list(m) for m in words]==cert['free_monomials']
 expansion=s.zeros(len(labels),len(words))
 for i,label in enumerate(labels):
  poly=s.Poly(s.prod(affine[j] for j in label),*vv)
  for j,word in enumerate(words):expansion[i,j]=poly.coeff_monomial(s.prod(vv[k] for k in word))
 residual=product(quotient,expansion);same(residual,s.Matrix(cert['all150_order3_coker_after_full_higher_variable_elimination']))
 er=cert['independent_coker_rows'];assert len(er)==3
 q=[s.expand(sum(c*s.prod(vv[i] for i in word) for c,word in zip(residual[r,:],words))) for r in er]
 conic=4*v0*v1+4*v0+7*v1*v1+18*v1+8
 expected=[2*v2*(2*v0+3*v1+5)/5,-conic/2,-v2*(v0+4*v1+5)/5]
 assert all(s.expand(a-b)==0 for a,b in zip(q,expected))
 basis=residual[er,:];bc=list(dm(basis).rref()[1]);same(product(residual[:,bc],basis[:,bc].inv())*basis,residual)
 H=4*v0+7*v1+11;cof=[-H/6,2*v2/3,-2*H/3]
 assert s.expand(sum(a*b for a,b in zip(cof,q))-v2)==0 and s.expand(-2*q[1]-conic)==0
 # Reverse containment: q1,q3 are multiples of v2, q2=-conic/2.
 assert s.cancel(q[0]/v2).is_polynomial(*vv) and s.cancel(q[2]/v2).is_polynomial(*vv)
 t=s.Symbol('t');param={v0:-7*t/4-1+3/(4*t),v1:t-1,v2:0}
 assert s.cancel(conic.subs(param))==0
 assert s.expand(conic-(v1+1)*(4*v0+7*v1+11)+3)==0
 # Thus v1+1 is a unit (inverse (4v0+7v1+11)/3), proving exact Gm.
 witness=s.Matrix([-2,0,0]);low=origin+kernel*witness
 values=s.Matrix([s.prod(low[j] for j in label) for label in labels]);intrinsic=product(R3,values)
 high=s.zeros(18,1);hs=-product(inverse,intrinsic[pr,:])
 for j,i in enumerate(pc):high[i]=hs[j]
 same(low,s.Matrix(cert['finite3jet_first_source_rho10']+cert['finite3jet_second_order_flag_B8']))
 same(high,s.Matrix(cert['finite3jet_higher_first_source_P2_free10']+cert['finite3jet_third_order_flag_C8']))
 B=sum((low[10+i]*E[i] for i in range(8)),s.zeros(5));C3=sum((high[10+i]*E[i] for i in range(8)),s.zeros(5))
 full=physical_operators(M,B,C3);P1=product(C[:,free],low[:10,:]);P2=product(C[:,free],high[:10,:])
 witness_outputs=[P0,act(full[1],P0)+P1,act(full[2],P0)+act(full[1],P1)+P2,
                  act(full[3],P0)+act(full[2],P1)+act(full[1],P2)]
 for value in witness_outputs:same(project(value),s.zeros(150,1));same(product(C,product(left,value[rows,:])),value)
 assert high[1]==-1 and free[1]==2 and kernel[1,:]==s.zeros(1,3)
 fun=-product(inverse[pc.index(1),:],R3[pr,:])
 extra=s.factor(sum(c*s.prod(affine[j] for j in label) for c,label in zip(fun,labels)))
 assert s.expand(extra-(v1+2)*(3*v0+7*v1+5)/2)==0
 files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',folder/'report.json',folder/'setup.json',folder/'certificate.json.gz',source,accepted]
 files += [ROOT/n for n in report['input_sha256']]
 result=dict(id=STEM,all_passed=True,independent_method='finite5x5 physical pushforward and inverse; no Lie/log/author helper imports',
  physical_fixed_CA150_by20_reconstructed=True,all9_relative_flags_and_free10=True,first_flag_quotient_rank=8,
  nilpotent_opposite_flag_N4_and_det1=True,all150_order2_matrix_rank=15,order2_columns=18,order2_kernel_rank=3,
  all150_order3_raw_and_full_higher_cokernel_reconstructed=True,order3_constant_cokernel_rank=3,
  exact_ideal_generators=['v2',str(conic)],both_ideal_containments_by_explicit_cofactors=True,
  necessary_curve_is_Gm=True,explicit_finite3jet_witness=['-2','0','0'],direct_all150_membership_through3=True,
  extra_CA_coefficient_index2_at_a2='-1',unique_higher_extra_CA_polynomial=str(extra),
  formal_or_algebraic_pair_branch_proved=False,order4_or_allorder_lift_claimed=False,
  nonordinary_pair_or_conductor_or_B18_full_odd_ideal_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-started)
 (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',STEM,result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
