#!/usr/bin/env python3
"""Full tangent obstructions and an actual MPS degeneration for rank-three Q.

All 150 quadratic bivector variables and 350 projective Jacobi equations are
retained.  The two exact obstruction functionals and the actual curve together
give a reduced-support argument, not a claim that the linear tangent plane is
integrable.  See RANK_THREE_QUADRIC_SUPPORT.md.
"""
import argparse,hashlib,itertools,json,resource,sys,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
sys.path.insert(0,str(ROOT/'legacy/v03/code'))
import poisson_p4_components_check as p
resource.setrlimit(resource.RLIMIT_CPU,(90,90))
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic()
h,u,v,z,w,a,b=s.symbols('h u v z w a b');xx=(h,u,v,z,w)
P=s.zeros(5)
source=json.loads((ROOT/'reports/CLASSIFICATION_RANK3_TANGENT_EXTRA_SOURCES_20261107.json').read_text())['result']['square']
for ij,f in source['witness_brackets'].items():
 i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=dict(zip(map(str,xx),xx)));P[j,i]=-P[i,j]
def sparse(M):
 return {(i,j):{e:int(c) for e,c in s.Poly(M[i,j],*xx).terms() if c}
         for i,j in p.PAIRS if M[i,j]!=0}
def vector(M):
 return s.Matrix([s.Poly(M[i,j],*xx).coeff_monomial(e) for (i,j),e in p.BIVECTOR_BASIS])
def matrix(col):
 R=s.zeros(5)
 for ((i,j),e),c in zip(p.BIVECTOR_BASIS,col):
  R[i,j]+=c*s.prod(x**k for x,k in zip(xx,e))
 for i,j in p.PAIRS:R[j,i]=-R[i,j]
 return R
def dm(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ)
D=s.Matrix(p.linearized_jacobi_matrix(sparse(P)).tolist())
K=dm(D).nullspace().to_Matrix().T
orbit=[]
for i,j in itertools.product(range(5),repeat=2):
 A=s.zeros(5);A[i,j]=1;X=A*s.Matrix(xx)
 orbit.append(vector(sum([X[k]*P.diff(xx[k]) for k in range(5)],s.zeros(5))-A*P-P*A.T))
M=s.Matrix(p.vertical_bivector_matrix().tolist()).row_join(s.Matrix.hstack(*orbit,vector(P)))
piv=dm(M.row_join(K)).rref()[1]
ks=[j-M.cols for j in piv if j>=M.cols]
assert len(ks)==2 and K.cols==49
directions=[]
for k in ks:
 col=K[:,k];den=s.ilcm(*[q.q for q in col]);col*=den
 gcd=s.igcd(*[int(c) for c in col]);col/=gcd
 directions.append(matrix(col))
R,S=directions
assert D*vector(R)==s.zeros(350,1) and D*vector(S)==s.zeros(350,1)
T=P+a*R+b*S
J={}
for i,j,k in p.TRIPLES:
 J[i,j,k]=s.expand(sum(T[i,l]*s.diff(T[j,k],xx[l])+T[j,l]*s.diff(T[k,i],xx[l])+T[k,l]*s.diff(T[i,j],xx[l]) for l in range(5)))
four={}
for ids in p.QUADRUPLES:
 four[ids]=s.expand(sum((-1)**q*xx[ids[q]]*J[tuple(t for t in ids if t!=ids[q])] for q in range(4)))
eq=[]
for f in four.values():
 for c in s.Poly(f,*xx).coeffs():
  c=s.factor(c)
  if c!=0 and c not in eq and -c not in eq:eq.append(c)
def pf(M,ids):
 i,j,k,l=ids
 return s.expand(M[i,j]*M[k,l]-M[i,k]*M[j,l]+M[i,l]*M[j,k])
divisor=s.factor(sum((-1)**q*xx[q]*pf(T,tuple(j for j in range(5) if j!=q)) for q in range(5)))
result={'source_sha256':hashlib.sha256((ROOT/'reports/CLASSIFICATION_RANK3_TANGENT_EXTRA_SOURCES_20261107.json').read_bytes()).hexdigest(),
 'kernel_dimension':K.cols,'orbit_vertical_scale_rank':len(piv)-2,'quotient_dimension':len(ks),
 'directions':[{f'{i},{j}':str(Q[i,j]) for i,j in p.PAIRS if Q[i,j]!=0} for Q in directions],
 'full_projective_Jacobi_equations':[str(q) for q in eq], 'projective_divisor':str(divisor),
 'elapsed_seconds':time.monotonic()-start}
O=s.Matrix([[s.expand(s.Poly(four[ids],*xx).coeff_monomial(e)).coeff(a,aa).coeff(b,bb)
             for aa,bb in [(2,0),(1,1),(0,2)]] for ids,e in p.FOUR_VECTOR_BASIS])
dpiv=dm(D).rref()[1];J0=D[:,list(dpiv)];rows=list(dm(J0.T).rref()[1]);minor=J0[rows,:]
minor_inverse=dm(minor).inv()
correction=minor_inverse.matmul(dm(O[rows,:])).to_Matrix()
obstructions=O-J0*correction
assert obstructions[rows,:]==s.zeros(len(rows),3)
canonical,opiv=dm(obstructions.T).rref()
genrows=list(opiv)
obstruction_generators=[s.factor(sum(obstructions[i,j]*m for j,m in enumerate((a*a,a*b,b*b)))) for i in genrows]
result['quadratic_obstruction']={'derivative_rank':len(dpiv),'augmented_rank':len(dm(D.row_join(O)).rref()[1]),
 'increment':len(opiv),'selected_derivative_columns':list(dpiv),'selected_derivative_rows':rows,
 'obstruction_rows':genrows,'generators':list(map(str,obstruction_generators)),
 'all_original_rows_retained':True}
functionals=[]
for i in genrows:
 left=s.zeros(1,350);left[0,i]=1
 coeff=dm(J0[i,:]).matmul(minor_inverse).to_Matrix()
 for j,k in enumerate(rows):left[0,k]-=coeff[0,j]
 assert left*D==s.zeros(1,150)
 assert left*O==obstructions[i,:]
 functionals.append({'row':i,'nonzero_entries':{str(j):str(left[0,j]) for j in range(350) if left[0,j]},
                     'quadratic_pairing':list(map(str,list(left*O)))})
result['quadratic_obstruction']['left_null_functionals']=functionals
if not opiv or obstructions[:,0]==s.zeros(350,1):
 Q=[]
 active=list(range(3)) if not opiv else [0]
 for j in active:
  col=s.zeros(150,1)
  for i,k in enumerate(dpiv):col[k]=-correction[i,j]
  assert D*col+O[:,j]==s.zeros(350,1)
  Q.append(matrix(col))
 Tc=T+a*a*Q[0]+a*b*Q[1]+b*b*Q[2] if not opiv else P+a*R+a*a*Q[0]
 Jc={}
 for i,j,k in p.TRIPLES:
  Jc[i,j,k]=s.expand(sum(Tc[i,l]*s.diff(Tc[j,k],xx[l])+Tc[j,l]*s.diff(Tc[k,i],xx[l])+Tc[k,l]*s.diff(Tc[i,j],xx[l]) for l in range(5)))
 residuals=set()
 for ids in p.QUADRUPLES:
  f=s.expand(sum((-1)**q*xx[ids[q]]*Jc[tuple(t for t in ids if t!=ids[q])] for q in range(4)))
  residuals.update(s.factor(q) for q in s.Poly(f,*xx).coeffs() if q!=0)
 result['second_order_corrections']=[{f'{i},{j}':str(Q0[i,j]) for i,j in p.PAIRS if Q0[i,j]!=0} for Q0 in Q]
 result['corrected_plane_exact_Poisson']=not residuals
 result['corrected_plane_residuals']=sorted(map(str,residuals))
 if not residuals:
  result['corrected_plane_projective_divisor']=str(s.factor(sum((-1)**q*xx[q]*pf(Tc,tuple(j for j in range(5) if j!=q)) for q in range(5))))
  # The exact curve, not its non-Poisson tangent plane, is used for this limit.
  assert active==[0]
  old_weights=(1,0,0,-1,-1)
  for N,target in [(P,1),(R,0),(Q[0],-1)]:
   assert all(sum(e[k]*old_weights[k] for k in range(5))-old_weights[i]-old_weights[j]==target
              for i,j in p.PAIRS for e,c in s.Poly(N[i,j],*xx).terms() if c)
  Pstar=Tc.subs(a,1)
  div=s.Matrix([sum(s.diff(Pstar[i,j],xx[i]) for i in range(5)) for j in range(5)])
  E=s.Matrix(xx);Pi=Pstar-(E*div.T-div*E.T)/5
  assert all(s.expand(sum(s.diff(Pi[i,j],xx[i]) for i in range(5)))==0 for j in range(5))
  forms=(4*h+4*u-w,-4*h+12*u+9*w,4*h-4*u+3*w,v,z)
  A=s.Matrix([[s.diff(f,x) for x in xx] for f in forms]);yy=s.symbols('x y l v z')
  back=dict(zip(xx,A.inv()*s.Matrix(yy)));weights=(0,0,1,0,-2)
  changed=(A*Pi*A.T).applyfunc(lambda f:s.expand(f.subs(back,simultaneous=True)))
  limit=s.zeros(5);hist={}
  for i,j in p.PAIRS:
   for e,c in s.Poly(changed[i,j],*yy).terms():
    if not c:continue
    wt=sum(e[k]*weights[k] for k in range(5))-weights[i]-weights[j]
    assert wt>=0
    hist[wt]=hist.get(wt,0)+1
    if wt==0:limit[i,j]+=c*s.prod(yy[k]**e[k] for k in range(5))
   limit[j,i]=-limit[i,j]
  limit_divisor=s.factor(sum((-1)**q*yy[q]*pf(limit,tuple(j for j in range(5) if j!=q)) for q in range(5)))
  assert s.expand(limit_divisor+16*s.prod(yy))==0
  result['orbit_limits']={'original_weights':list(old_weights),'P_R_Q_weights':[1,0,-1],
   'coordinate_forms':list(map(str,forms)),'coordinate_determinant':str(A.det()),
   'divergence_free_gauge_vector':list(map(str,div/5)),
   'SNC_weights':list(weights),'transformed_weight_histogram':hist,'negative_terms':0,
   'transformed_brackets':{f'{i},{j}':str(changed[i,j]) for i,j in p.PAIRS if changed[i,j]!=0},
   'limit_brackets':{f'{i},{j}':str(limit[i,j]) for i,j in p.PAIRS if limit[i,j]!=0},
   'limit_divisor':str(limit_divisor)}
result['elapsed_seconds']=time.monotonic()-start
result['producer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
result['status']='exact support inputs verified; formal uniqueness proof is in the note'
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'kernel_dimension':K.cols,'quadratic_obstruction_generators':result['quadratic_obstruction']['generators'],
                  'actual_curve':result.get('corrected_plane_exact_Poisson'),
                  'SNC_limit':result.get('orbit_limits',{}).get('limit_divisor'),
                  'elapsed_seconds':result['elapsed_seconds']},indent=2))
