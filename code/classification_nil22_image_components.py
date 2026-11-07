#!/usr/bin/env python3
"""Reconstruct the true quadratic-support certificates of two nil22 sources.

No producer import, higher jet, or local-uniqueness inference is made here.
The universal source and first quadratic report are immutable arithmetic inputs.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
from classification_three_hyperplanes_quadric_independent import ambient_matrices
ROOT=Path(__file__).resolve().parents[1]

def main():
 started=time.monotonic()
 source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
 quad_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_QUADRATIC_20261107.json'
 sources=json.loads(source_path.read_text())['cases'];quadratics=json.loads(quad_path.read_text())['cases']
 xx=s.symbols('h x1:5');h,q1,q2,x,y=xx;mm=s.symbols('m0:6');env={str(z):z for z in xx+mm}
 dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
 rank=lambda A:len(dm(A).rref()[1])
 def decode(data):
  P=s.zeros(5)
  for ij,f in data.items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
  return P
 def sparse(P):return {ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())for ij in PAIRS if P[ij]!=0}
 records=[]
 for index,label,d_rank,n_count in [(4,'triple_A',101,1),(5,'triple_B',100,2)]:
  packet=sources[index];saved=quadratics[index];assert label in packet['name']
  family=decode(packet['universal_brackets']);values={env[z]:s.Rational(a)for z,a in packet['full_source_specification']['witness'].items()}
  P=decode(packet['witness_brackets']);assert family.subs(values)==P
  directions=[family.diff(z).subs(values)for z in mm]
  D,V,G,F=ambient_matrices(P,xx,directions)
  assert (rank(D),rank(V),rank(G),rank(F))==(d_rank,24,48,48)
  assert D*F==s.zeros(350,F.cols)
  assert jacobi(sparse(family))==s.zeros(350,1)
  N=list(map(decode,saved['normals']));assert len(N)==n_count
  NN=s.Matrix.hstack(*[s.Matrix([sparse(R).get(ij,{}).get(ex,0)for ij,ex in COLUMNS])for R in N])
  assert D*NN==s.zeros(350,n_count)
  assert rank(F.row_join(NN))==150-d_rank
  Qcols=[jacobi(sparse(N[i])) if i==j else jacobi(sparse(N[i]+N[j]))-jacobi(sparse(N[i]))-jacobi(sparse(N[j]))for i,j in itertools.combinations_with_replacement(range(n_count),2)]
  Q=s.Matrix.hstack(*Qcols);funcs=[]
  for item in saved['quadratic_cokernel_functionals']:
   row=s.zeros(1,350)
   for i,a in item['nonzero_coefficients'].items():row[int(i)]=s.Rational(a)
   assert row*D==s.zeros(1,150)
   assert list(row*Q)==list(map(s.Rational,item['quadratic_pairing']))
   funcs.append(row)
  L=s.Matrix.vstack(*funcs);pair=L*Q
  normalized=pair[:,:n_count].inv()*L
  expected=s.Matrix([[1]]) if label=='triple_A' else s.Matrix([[1,0,-s.Rational(3,124)],[0,1,s.Rational(5,31)]])
  assert normalized*D==s.zeros(n_count,150);assert normalized*Q==expected
  coordinates=list(dm(D).rref()[1]);rows=list(dm(D[:,coordinates].T).rref()[1]);minor=dm(D[rows,coordinates]).det();assert minor!=0
  coordinate=s.zeros(150,d_rank)
  for j,i in enumerate(coordinates):coordinate[i,j]=1
  cut=coordinate.row_join(NN);assert cut.cols==102;assert rank(F.row_join(cut))==150
  a,b=s.symbols('a b')
  if label=='triple_A':
   generators=[a*a];radical_nonzero='1';radical_identity='a^2 is a generator'
  else:
   A=s.Rational(3,124);B=s.Rational(5,31);generators=[a*a-A*b*b,a*b+B*b*b]
   radical_nonzero=s.factor(B*B-A);assert radical_nonzero!=0
   assert s.expand(b*generators[0]-(a-B*b)*generators[1]-radical_nonzero*b**3)==0
   radical_identity='b*q1-(a-B*b)*q2=(B^2-A)*b^3, B^2-A!=0'
  # Full unscaled affine source and the actual quintic, independent of all6 m.
  Q0=decode(packet['full_six_magnetic_source']);aff=Q0.subs(h,1)
  modular=s.Matrix([s.expand(sum(s.diff(aff[i,j],xx[i])for i in range(1,5)))for j in range(1,5)])
  assert modular==s.Matrix([x,y,0,0])
  pf=s.expand(Q0[1,2]*Q0[3,4]-Q0[1,3]*Q0[2,4]+Q0[1,4]*Q0[2,3])
  cubic=x**3-9*h*x*x-3*h*y*y+108*h**3 if label=='triple_A' else x**3+h*x*y-h**3
  assert s.expand(h*pf+h*h*cubic)==0
  assert len(s.factor_list(cubic)[1])==1
  infinity={h:0,x:0,y:1}
  gradient=[s.diff(cubic,z)for z in(h,x,y)];at_infinity=[z.subs(infinity)for z in gradient]
  if label=='triple_A':
   assert at_infinity==[-3,0,0]
   singular={h:1,x:6,y:0};assert cubic.subs(singular)==0;assert all(z.subs(singular)==0 for z in gradient)
   assert s.hessian(cubic.subs(h,1),(x,y)).det().subs(singular)!=0
   incidence='The unique double hyperplane h=0 avoids the cubic node and meets it at one smooth flex with multiplicity3.'
  else:
   assert at_infinity==[0,0,0]
   local=s.Poly(cubic.subs(y,1),h,x);tangent=sum(co*h**ex[0]*x**ex[1]for ex,co in local.terms()if sum(ex)==2)
   assert tangent==h*x
   incidence='The unique double hyperplane h=0 passes through the cubic node and is one of its two tangents.'
  records.append(dict(name=packet['name'],source_case_index=index,all_six_actual_parameters_retained=True,
   full_D_rank=d_rank,vertical_rank=24,actual_group_Euler_scalar_rank=48,actual_F_rank=48,full_kernel_dimension=150-d_rank,true_normal_count=n_count,
   full_D_times_F_and_N_zero=True,F_plus_N_rank=150-d_rank,true_literal_slice_dimension=102,actual_cut_count=48,
   literal_slice_coordinate_columns=coordinates,implicit_D_rows=rows,implicit_D_minor=str(minor),F_plus_literal_slice_rank=150,
   normalized_left_functionals=[{str(i):str(a)for i,a in enumerate(normalized.row(j))if a}for j in range(n_count)],
   normals=saved['normals'],all_D_columns_annihilated=True,all_normal_quadratic_pairings=[[str(a)for a in row]for row in expected.tolist()],
   quadratic_generators=list(map(str,generators)),radical_nonzero_coefficient=str(radical_nonzero),radical_identity=radical_identity,
   reduced_true_slice_support_is_origin=True,whole_source_is_full_component=True,projective_component_dimension=23,full_ambient_tangent_dimension=125-d_rank,
   local_uniqueness_or_smoothness_claimed=False,whole_source_affine_modular=list(map(str,modular)),whole_source_actual_quintic=str(-h*h*cubic),
   cubic_irreducible=True,infinity_gradient=list(map(str,at_infinity)),intrinsic_double_H_incidence=incidence))
 inputs=[source_path,quad_path,Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py',ROOT/'code/classification_three_hyperplanes_quadric_independent.py']
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),cases=records,
  support_logic='The true102-dimensional slice is cut by48 independent linear forms complementary to actualF48. Its independentD rows eliminate101 or100 coordinates formally. Verified fullD-annihilating functionals yield maximal-radical initial equations in the remaining1 or2 variables. Hence the cut has dimension0 and Krull bounds the fullambient localdimension by48, matching the actual source lower48. Its irreducible PGL-saturated closure is a full23-dimensional component. No local uniqueness, smoothness or exact nil-thickness follows.',
  known94_distinction='Dimensions other than23 are excluded after componenthood. All prior23-dimensional components other than R4_ZERO_LOG_5H have a reduced rank4 degeneracy divisor (MPS have an SNC point), whereas each present whole closure satisfies the proper closed condition D=H^2F3. The zero-log23 whole closure has D=H^5, also closed, and the present sources have an irreducible nodal cubic factor. Thus neither equals a prior component. The added half/third25 components do not affect this comparison.',
  pair_distinction='Both whole6-magnetic sources have unique doubleH and unique irreducible cubic. In triple_A H avoids the cubic node; in triple_B H passes through its node. If the two23-dimensional component closures coincided, their PGL-saturated constructible dense images would contain intersecting open subsets, contradicting this invariant incidence.',
  independent_review_required=True,seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_COMPONENTS_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(all_passed=True,new_component_candidates=[r['name']for r in records],seconds=result['seconds'])))
if __name__=='__main__':main()
