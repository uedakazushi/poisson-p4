#!/usr/bin/env python3
"""Independent full-rank1 tensors and exact target-coordinate cotangent audit."""
import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def same(a,b):return {k:v for k,v in a.items()if k!='seconds'}=={k:v for k,v in b.items()if k!='seconds'}
def parse(v,local):return s.sympify(re.sub(r'\blambda\b','_lambda',str(v)),locals=local)
def matrix(entries,local):return s.Matrix([[parse(v,local)for v in row]for row in entries])

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--rank2-replay',type=Path,required=True)
 ap.add_argument('--whole-replay',type=Path,required=True)
 ap.add_argument('--tangent-replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 names=['AMBIENT_RANK2_INTERSECTION_UNIQUE','AMBIENT_WHOLE_REDUCED_INTERSECTION','AMBIENT_INTERSECTION_TANGENT_EXCESS']
 paths=[Path(__file__)];reports=[]
 for name,replay in zip(names,[args.rank2_replay,args.whole_replay,args.tangent_replay]):
  path=ROOT/f'reports/{name}_20261103_v63/report.json';d=load(path);assert same(d,load(replay));reports.append(d);paths.append(path)
  for file,h in d['input_sha256'].items():
   pth=ROOT/file;assert sha(pth)==h,file;paths.append(pth)
 rank2,whole,tan=reports
 x,y,z,w=s.symbols('x y z w');Y=s.Matrix([x,y,z,w]);X=s.symbols('x1:5')
 p,sv,r,lam=s.symbols('p s r lambda');k=s.symbols('k0:9');free=s.symbols('c0 c2 c5')
 local={str(v):v for v in [*Y,*X,p,sv,r,lam,*k,*free]}
 local['_lambda']=lam
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 # All nonzero critical-rank2 points reduce to B=C=1. Check the two
 # actual finite covariant tensor transformations without producer helpers.
 A,B,C,a,b=s.symbols('A B C a b',nonzero=True);beta=s.zeros(4)
 beta[0,2]=A*x+2*C*z;beta[0,3]=B*x;beta=beta-beta.T
 shear=s.Matrix([[1,0,0,0],[0,1,0,-A/(2*C)],[A/(2*C),0,1,0],[0,0,0,1]])
 assert (shear.T*J*shear-J).applyfunc(s.cancel)==s.zeros(4)
 inv=shear.inv();after=(inv.T*beta.subs(dict(zip(Y,inv*Y)),simultaneous=True)*inv).applyfunc(s.cancel)
 assert after==beta.subs(A,0)
 torus=s.diag(a,1/a,b,1/b);inv=torus.inv()
 after=(inv.T*after.subs(dict(zip(Y,inv*Y)),simultaneous=True)*inv).applyfunc(s.cancel)
 assert after==(beta.subs({A:0,B:B*b/a**2,C:C/(a*b*b)},simultaneous=True)).applyfunc(s.cancel)
 assert rank2['canonical_full_projective_flag_fiber']==['r','q','w**2']
 fb=load(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json')
 bi,qi=fb['beta_indices_zero_based'],fb['q_indices_zero_based'];adapter=s.Matrix(fb['beta16_to_intrinsic24_adapter'])
 physical={};rank1_records={}
 for name,body in whole['families'].items():
  Q=matrix(body['actual_quadratic_matrix'],local).subs(dict(zip(X,Y)),simultaneous=True)
  assert Q+Q.T==s.zeros(4)
  flag=s.Matrix([parse(v,local)for v in body['flag']])
  div=s.Matrix([s.expand(sum(s.diff(Q[i,j],Y[i])for i in range(4)))for j in range(4)])
  factor=2 if name.endswith('_CA')else 4;linear=x+lam*z if name.startswith('y_chart')else z
  assert (div-factor*r*linear*flag).applyfunc(s.expand)==s.zeros(4,1)
  for i,j,h in itertools.combinations(range(4),3):
   assert s.expand(flag[i]*Q[j,h]-flag[j]*Q[i,h]+flag[h]*Q[i,j]).subs(r,0)==0
  u=s.Matrix([parse(v,local)for v in body['normalized_61_coordinates']])
  assert adapter*u[bi,0]==s.Matrix(whole['canonical_rank1_beta'])
  assert u[qi,0].subs(dict.fromkeys([r,*k],0))==s.zeros(45,1)
  qM=u[qi,0].jacobian([r,*k]);assert qM==matrix(body['quadratic_q_direction_matrix'],local)
  assert (u[qi,0]-qM*s.Matrix([r,*k])).applyfunc(s.expand)==s.zeros(45,1)
  assert body['zero_perspective']
  physical[name]=Q;rank1_records[name]=dict(full_divergence_identity=True,zero_q_origin=True,r_zero_has_common_flag_factor=True)
 for chart in ['y_chart','off_point']:
  assert (physical[chart+'_CA']-physical[chart+'_CB']).subs(r,0).applyfunc(s.expand)==s.zeros(4)
 # Direct inverse in ACTUAL target K, not the producer's truncated source
 # recursion: all selected blocks have constant nonzero determinant.
 comp=load(ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json');matrices={};tangent_records={}
 for name,body in comp['families'].items():
  pars=[local[v]for v in body['source_fiber_parameters']]
  u=s.Matrix([s.sympify(v,locals=local)for v in body['full_61_coordinates']]);qv=u[qi,0]
  M=u.jacobian(pars);base=u.subs(dict.fromkeys(pars,0));selected=body['actual_q_selected_rows'];targetrows=[qi[i]for i in selected]
  L=M[targetrows,:];L0=L.subs({p:0,sv:0});assert s.expand(L.det()-L0.det())==0 and L0.det()!=0
  ks=s.symbols(f'K0:{len(pars)}');K=s.Matrix(ks);offset=base[targetrows,0].subs({p:0,sv:0})
  source=L.inv()*(L0*K+offset-base[targetrows,0])
  exact=(M*source+base).applyfunc(s.cancel)
  assert all(s.denom(v)==1 or s.denom(v).is_Number for v in exact)
  assert exact[targetrows,0]==L0*K+offset
  origin={p:0,sv:0,**dict.fromkeys(ks,0)}
  for j,v in enumerate(pars):
   if str(v)in ['k0','k2','k5']:origin[ks[j]]=free[['k0','k2','k5'].index(str(v))]
  regular=list(ks)+[sv]
  columns=[exact.diff(p,n).subs(origin)/s.factorial(n)for n in [2,3]]
  columns += [exact.diff(v).subs(origin)for v in regular]
  columns += [exact.diff(v).diff(p).subs(origin)for v in regular]
  T=s.Matrix.hstack(*columns).applyfunc(s.expand);rank=T.rank();expected=19 if name=='N'else 21
  saved=matrix(tan['families'][name]['beta_slice_embedded_tangent_matrix'],local)
  assert rank==expected and T.row_join(saved).rank()==expected
  assert all(s.Poly(v,*free).total_degree()<=1 for v in T)
  matrices[name]=T;tangent_records[name]=dict(independent_exact_target_K_inverse=True,
   slice_tangent_dimension=rank,whole_S61_tangent_dimension=rank+7,
   saved_matrix_column_span_exactly_matches=True,matrix_coefficients_linear_in_free_parameters=True,
   independently_rebuilt_matrix_sha256=hashlib.sha256(json.dumps(list(map(str,T))).encode()).hexdigest())
 union=matrices['CA'].row_join(matrices['CB']);assert union.rank()==22
 assert matrices['N'].row_join(matrices['CA']).rank()==21
 assert matrices['N'].row_join(matrices['CB']).rank()==21
 for values in [(0,0,0),(1,0,0),(2,3,5)]:assert union.subs(dict(zip(free,values))).rank()==22
 paths += [ROOT/f'research/notes/{n}_20261103_v63.md'for n in names[:2]]
 paths += [ROOT/'research/notes/AMBIENT_INTERSECTION_NONREDUCED_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_INTERSECTION_SYMPLECTIC_SCOPE_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_INTERSECTION_20261103_v63_REVIEW.md',
  ROOT/'reports/AMBIENT_PINCHED_MODELS_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_COMPONENT_R1_NONCM_20261103_v63_REVIEW/report.json']
 result=dict(id='AMBIENT_INTERSECTION_20261103_v63_REVIEW',all_passed=True,fresh_QQ_three_producer_replays_match=True,
  critical_rank2_finite_tensor_reduction_independently_checked=True,rank1_full_families=rank1_records,
  whole_reduced_intersection_equals_N_accepted='WITHIN ALL S61; after restoration only generic-Poisson-rank4 open',
  tangent_families=tangent_records,symbolic_union_rank=22,all_specializations_union_rank_upper=22,
  actual_intersection_tangent_S61_lower=27,reduced_N_tangent_S61_dimension=26,
  actual_scheme_intersection_nonreduced_along_common_nonnormal_stratum_accepted=True,
  full_scheme_intersection_ideal_or_nilpotence_order_claimed=False,
  global_projective_intersection_components_supported_at_generic_Poisson_rank_at_most2_excluded=False,
  input_sha256={str(path.relative_to(ROOT)):sha(path)for path in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS whole S61 reduced intersection and actual target cotangent excess')

if __name__=='__main__':main()
