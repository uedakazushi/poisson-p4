#!/usr/bin/env python3
"""Reconstruct the full true two-normal certificate for the nil31 third source.

The previous complete source and quadratic vectors remain immutable inputs.
This extracts explicit normalized left functionals and a maximal-radical
identity, without computing a higher jet or claiming local uniqueness.
"""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
from classification_three_hyperplanes_quadric_independent import ambient_matrices
ROOT=Path(__file__).resolve().parents[1]

def main():
 started=time.monotonic();source_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json';quad_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_THIRD_QUADRATIC_SUPPORT_20261107.json'
 packet=json.loads(source_path.read_text())['branches'][1];saved=json.loads(quad_path.read_text());xx=s.symbols('h x1:5');env=dict(zip(map(str,xx),xx))
 def decode(data):
  P=s.zeros(5)
  for ij,f in data.items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
  return P
 def sparse(P):return{ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())for ij in PAIRS if P[ij]!=0}
 P=decode(packet['witness_brackets']);directions=list(map(decode,packet['all_ten_actual_parameter_derivatives']))
 D,V,G,F=ambient_matrices(P,xx,directions);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
 assert D*F==s.zeros(350,F.cols);assert len(dm(D).rref()[1])==98;assert len(dm(F).rref()[1])==50
 N=list(map(decode,saved['normals']));assert len(N)==2
 vectors=[]
 for R in N:vectors.append(s.Matrix([sparse(R).get(ij,{}).get(ex,0)for ij,ex in COLUMNS]))
 NN=s.Matrix.hstack(*vectors);assert D*NN==s.zeros(350,2);assert len(dm(F.row_join(NN)).rref()[1])==52
 Q=s.Matrix.hstack(jacobi(sparse(N[0])),jacobi(sparse(N[0]+N[1]))-jacobi(sparse(N[0]))-jacobi(sparse(N[1])),jacobi(sparse(N[1])))
 funcs=[]
 for item in saved['quadratic_cokernel_functionals']:
  row=s.zeros(1,350)
  for i,a in item['nonzero_coefficients'].items():row[int(i)]=s.Rational(a)
  assert row*D==s.zeros(1,150);assert list(row*Q)==list(map(s.Rational,item['quadratic_pairing']));funcs.append(row)
 L=s.Matrix.vstack(*funcs);pair=L*Q;assert pair.shape==(2,3)
 normalized=pair[:,:2].inv()*L;expected=s.Matrix([[1,0,-s.Rational(113435,3825851)],[0,1,s.Rational(214061,7651702)]])
 assert normalized*D==s.zeros(2,150);assert normalized*Q==expected
 t0,t1=s.symbols('t0 t1');A=-expected[0,2];B=expected[1,2];q1=t0*t0-A*t1*t1;q2=t0*t1+B*t1*t1;nonzero=s.factor(B*B-A);assert nonzero!=0
 assert s.expand(t1*q1-(t0-B*t1)*q2-nonzero*t1**3)==0
 columns=list(dm(D).rref()[1]);rows=list(dm(D[:,columns].T).rref()[1]);minor=dm(D[rows,columns]).det();assert minor!=0
 # Literal complementary slice:98 independent coordinate columns plus2 normals.
 coordinate=s.zeros(150,98)
 for j,i in enumerate(columns):coordinate[i,j]=1
 cut=coordinate.row_join(NN);assert len(dm(F.row_join(cut)).rref()[1])==150
 r=-s.Rational(1,3);h,a,b,c,d=xx;u,v,t,q=s.symbols('u v t q');AA=s.symbols('A0:6');localenv={**env,**{str(z):z for z in(u,v,t,q)+AA}}
 whole=decode(packet['whole_homogeneous_h_Casimir_brackets']).subs(h,1)
 modular=s.Matrix([s.expand(sum(s.diff(whole[i,j],xx[i])for i in range(1,5)))for j in range(1,5)])
 assert modular==s.Matrix([0,-2*c/3,-s.Rational(2,3),0])
 assert len(set(modular.free_symbols)-{c})==0
 here=Path(__file__);inputs=[source_path,quad_path,here,ROOT/'code/classification_rank3_mixed_component_review.py',ROOT/'code/classification_three_hyperplanes_quadric_independent.py']
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=here.read_text(),
  full_D_rank=98,actual_F_rank=50,full_kernel_dimension=52,true_normal_count=2,full_D_times_N_and_F_zero=True,F_plus_N_rank=52,
  true_literal_slice_dimension=100,actual_cut_count=50,literal_slice_coordinate_columns=columns,implicit_D_rows=rows,implicit_D_minor=str(minor),F_plus_literal_slice_rank=150,
  normals=saved['normals'],normalized_left_functionals=[{str(i):str(a)for i,a in enumerate(row)if a}for row in [normalized.row(0),normalized.row(1)]],
  all_D_columns_annihilated=True,all_three_quadratic_pairings=[[str(a)for a in row]for row in expected.tolist()],
  quadratic_generators=list(map(str,[q1,q2])),radical_nonzero_coefficient=str(nonzero),radical_identity='t1*q1-(t0-B*t1)*q2=(B^2-A)*t1^3, with B^2-A!=0',reduced_true_slice_support_is_origin=True,
  whole_source_is_full_component=True,projective_component_dimension=25,full_ambient_tangent_dimension=27,local_uniqueness_or_smoothness_claimed=False,
  intrinsic_triple_hyperplane='h=0 in the normalized projective family, equivalently old_z=0',whole_source_modular_on_triple_H_infinity=list(map(str,modular)),
  entire_modular_type='A rank1, A^2=0, constant component outside imA in the length2 chain; affine modular field has no zero for any source parameter.',
  half_distinction='The half component has unique triple hyperplane old_h=0 and its entire dense source modular field there is old_y*d_old_x+old_z*d_old_y, nil31 with zero set and constant in imA. The third entire dense source triple-H modular field is nil211-long and has no zero. PGL equivalence must preserve the intrinsic triple hyperplane and becomes affine in its chart, preserving these types. Two same-dimensional component closures therefore cannot coincide.',
  known92_distinction='After componenthood: all dimensions other than25 are excluded. The third full rank4 closure has closed divisor H^3L^2; this excludes the SNC-meeting T01/T39/T40, positiveNEW25 with rank5 smooth residual quadric, and CN31 generic independent221 divisor. Generic-rank2 components are also excluded by closed rank.',
  support_logic='The literal100-dimensional slice is cut by50 independent linear functions complementary to actualF50. Its98 independent D rows eliminate98 coordinates formally; the two verified left functionals have initial forms q1,q2 with maximal radical. Thus the cut has dimension0. Krull bounds the full150-coordinate local dimension by50, equal to the actual source lower50. Its irreducible PGL-saturated closure is a full25-dimensional component. No local uniqueness, smoothness, product decomposition or exact nil-thickness is inferred.',seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_THIRD_COMPONENT_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,component25=True,radical_coefficient=str(nonzero),seconds=result['seconds'])))
if __name__=='__main__':main()
