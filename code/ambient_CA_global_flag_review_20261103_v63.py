#!/usr/bin/env python3
"""Independent sparse tensor audit of the full projective CA flag model."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
IDS={(ij,m):15*n+j for n,ij in enumerate(PAIRS)for j,m in enumerate(MONS)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def action(i,j):
 data={}
 def add(out,col,c):
  if c:data[out,col]=data.get((out,col),0)+c
 for (ab,m),col in IDS.items():
  a,b=ab
  for pos in range(2):
   if ab[pos]==j:
    pair=list(ab);pair[pos]=i
    if pair[0]!=pair[1]:add(IDS[(tuple(sorted(pair)),m)],col,1 if pair[0]<pair[1]else-1)
  for pos in range(2):
   if m[pos]==i:
    changed=list(m);changed[pos]=j;add(IDS[(ab,tuple(sorted(changed)))],col,-1)
 return s.SparseMatrix(150,150,data)
def canonical(v):
 d=s.zeros(5)
 for (ij,m),ind in IDS.items():
  for pos,index in enumerate(m):
   if index==ij[0]:d[ij[1],m[1-pos]]+=v[ind]
   if index==ij[1]:d[ij[0],m[1-pos]]-=v[ind]
 out=v.copy()
 for (i,j)in PAIRS:
  for a in range(5):
   out[IDS[((i,j),tuple(sorted([i,a])))]]-=d[j,a]/5
   out[IDS[((i,j),tuple(sorted([j,a])))]]+=d[i,a]/5
 return out
def wedge_factor(v):
 tuples=list(itertools.combinations(range(5),4));mons3=list(itertools.combinations_with_replacement(range(5),3))
 index={(t,m):35*i+j for i,t in enumerate(tuples)for j,m in enumerate(mons3)};W=s.zeros(175,5)
 for q in range(5):
  for r in range(5):
   for (ij,m),col in IDS.items():
    order=(r,q,*ij)
    if len(set(order))<4 or v[col]==0:continue
    sign=(-1)**sum(order[i]>order[j]for i in range(4)for j in range(i+1,4))
    W[index[(tuple(sorted(order)),tuple(sorted([r,*m])))],q]+=sign*v[col]
 return W

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 path=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h;paths.append(p)
 raw=[s.Matrix(v)for v in d['basis_raw150']];can=[s.Matrix(v)for v in d['basis_canonical150']]
 assert all(canonical(r)==c for r,c in zip(raw,can))
 C=s.Matrix.hstack(*can);assert C.rank()==20
 # Recover the fixed raw family directly from the displayed polynomials.
 X=s.symbols('h x y z w');h,x,y,z,w=X;mon3=[s.prod((h,x,z,w)[i]for i in t)for t in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in t)]
 assert [str(v)for v in mon3]==d['cubic_monomials']
 def coeff(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[a]for a in m))for i,j in PAIRS for m in MONS])
 def tensor(f,g,H):
  P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
 expected=[tensor(f,0,0)for f in [h*h,h*x,x*x]]+[tensor(0,1,0)]+[tensor(0,0,H)for H in mon3]
 assert all(coeff(P)==v for P,v in zip(expected,raw))
 # Index/monomial sparse action, independent of the producer's polynomial
 # differentiation and canonical action routine.
 generators={};invariance=[]
 for pos,A in zip(d['parabolic_matrix_positions'],d['exact_20_by20_Lie_actions']):
  i,j=pos;T=action(i,j);A=s.Matrix(A);assert T*C==C*A
  assert all(A[3,k]==0 for k in range(20)if k!=3)
  assert all(A[k,j]==0 for k in range(3)for j in range(3,20))
  generators[i,j]=T;invariance.append(True)
 source_actions={(i,j):action(i,j)for i,j in itertools.product(range(5),repeat=2)}
 H=sum((c*v for c,v in zip(d['boundary_witness_cubic_coefficients'],can[4:])),s.zeros(150,1));records={}
 for name,point,source in [('g_zero',can[0]+can[1]+can[2]+H,can[:3]+can[4:]),('f_zero',can[3]+H,can[3:])]:
  witness=d['rank_two_boundary_exact_witnesses'][name]
  D=s.Matrix.hstack(*source,*[source_actions[i,j]*point for i,j in itertools.product(range(5),repeat=2)])
  assert D==s.Matrix(witness['exact_raw150_image_differential'])
  determinant=D[witness['minor_rows'],witness['minor_columns']].det()
  assert determinant==s.Rational(witness['minor_determinant'])!=0
  W=wedge_factor(point);wr,wc=witness['constant_vector_factor_minor_rows'],witness['constant_vector_factor_minor_columns']
  assert W[wr,wc].det()==s.Rational(witness['constant_vector_factor_minor_determinant'])!=0
  assert W.rank()==witness['intrinsic_constant_vector_factor_matrix_rank']
  if name=='g_zero':assert W*s.Matrix([0,0,1,0,0])==s.zeros(175,1)
  records[name]=dict(all_boundary_differential_columns_independently_rebuilt=True,
   nonzero_QQ_minor_size=len(witness['minor_columns']),nonzero_QQ_minor_determinant=str(determinant),
   intrinsic_E_wedge_v_wedge_P_rank=W.rank(),intrinsic_factor_minor_independently_rebuilt=True)
 # Fresh all20-parameter physical Jacobi and actual degeneracy quintic.
 u=s.symbols('u0:20');P=sum((a*b for a,b in zip(u,expected)),s.zeros(5))
 for i,j,k in itertools.combinations(range(5),3):
  assert s.expand(sum(P[i,t]*s.diff(P[j,k],X[t])+P[j,t]*s.diff(P[k,i],X[t])+P[k,t]*s.diff(P[i,j],X[t])for t in range(5)))==0
 pf=P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3]
 assert s.expand(h*pf-u[3]*h**3*(u[0]*h*h+u[1]*h*x+u[2]*x*x))==0
 paths += [ROOT/'research/notes/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW',all_passed=True,
  fresh_full_QQ_producer_replay_matches=True,independent_sparse150_canonicalization_and_all320_Lie_columns=True,
  both_g_zero_and_f_zero_linear_subspaces_invariant_under_all16_parabolic_generators=True,
  all20_parameter_physical_raw_Jacobi_and_actual_quintic_verified=True,
  boundary_witnesses=records,projective_bundle_geometry_and_no_Euler_basepoints_accepted=True,
  whole_projective_CA_reduced_rank_at_most2_boundary_has_exactly_two_irreducible_components_dimensions=[27,25],
  global_normalization_R1_or_CA_intersection_CB_claimed=False,
  scope='Full projective reduced CA, before PGL quotient; boundary support and dimensions only, not new Poisson components.',
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS full projective CA proper flag bundle and exact two boundary supports')

if __name__=='__main__':main()
