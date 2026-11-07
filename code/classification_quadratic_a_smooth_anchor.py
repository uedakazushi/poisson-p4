#!/usr/bin/env python3
"""One weakest actual dimension squeeze at the SAME old lambda7 anchor.

The joint actual family uses inverse magnetic normalizers of the all-lambda
Poisson curve. No new anchor, local Groebner basis, unit solve or jet is used.
"""
from pathlib import Path
import hashlib,itertools as it,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 begun=time.monotonic();xx=s.symbols('h x1 x2 x3 x4');h,q1,q2,x,y=xx;sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';source=json.loads(sourcepath.read_text())['cases'][7];params=s.symbols('m0:6')+(s.Symbol('lam'),);env={str(v):v for v in xx+params};whole=s.zeros(5)
 for ij,text in source['universal_brackets'].items():i,j=map(int,ij.split(','));whole[i,j]=s.sympify(text,locals=env);whole[j,i]=-whole[i,j]
 values=dict(zip(params,[1,0,0,0,0,0,7]));P=whole.subs(values);frozen=json.loads((ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json').read_text())['frame'];assert frozen['D_rank']==101 and frozen['actual_full_source_F_rank']==48
 for ij,text in frozen['canonical_point'].items():i,j=map(int,ij.split(','));assert s.expand(P[i,j]-s.sympify(text,locals=env))==0
 arcpath=ROOT/'reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SPECTRAL_ARC_20261107.json';arc=json.loads(arcpath.read_text())['curves'][1];assert arc['status']=='ALL_LAMBDA_ACTUAL_POISSON_CURVE_FOUND';ep,lam=s.symbols('epsilon lambda');uu=s.symbols('u1 u2 v1 v2');ss=(s.Symbol('h'),*uu);env2={str(v):v for v in ss+(ep,lam)};env2['lam']=lam;family=s.zeros(5)
 for ij,text in arc['actual_regular_homogeneous_family'].items():i,j=map(int,ij.split(','));family[i,j]=s.sympify(text.replace('lambda','lam'),locals=env2);family[j,i]=-family[i,j]
 C=s.Matrix([h,q1,x,y,y-q2]).jacobian(s.Matrix(xx));back=C.inv();moved=(2*back*family*back.T).subs(dict(zip(ss,C*s.Matrix(xx))),simultaneous=True).subs(lam,7).applyfunc(s.expand);assert(moved.subs(ep,0)-P).applyfunc(s.expand)==s.zeros(5);R=moved.diff(ep).subs(ep,0)
 def vector(Q):return s.Matrix([s.Poly(Q[ij],*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def sparse(Q):return{ij:dict(s.Poly(Q[ij],*xx).terms())for ij in PAIRS if Q[ij]!=0}
 assert jacobi(sparse(P))==s.zeros(350,1);DR=(jacobi(sparse(P+R))-jacobi(sparse(P))-jacobi(sparse(R))).applyfunc(s.expand);assert DR==s.zeros(350,1)
 E=s.Matrix(xx);columns=[]
 for i,j in it.product(range(5),repeat=2):
  A=s.zeros(5);A[i,j]=1;field=A*E;columns.append(vector(E*field.T-field*E.T));columns.append(vector((A*P+P*A.T-xx[j]*P.diff(xx[i])).applyfunc(s.expand)))
 columns.append(vector(P));directions=[whole.diff(v).subs(values)for v in params];columns.extend(map(vector,directions));F=s.Matrix.hstack(*columns);actual=F.row_join(vector(R));dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);rank=len(dm(F).rref()[1]);pc=list(dm(actual).rref()[1]);pr=list(dm(actual[:,pc].T).rref()[1]);assert rank==48 and len(pc)==49;minor=dm(actual[pr,pc]).det();assert minor!=0
 paths=[Path(__file__),sourcepath,ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json',arcpath,ROOT/'reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_REVIEW_20261107.json',ROOT/'code/classification_rank3_mixed_component_review.py'];out=dict(task_id='E_multiple_hyperplane',all_passed=True,same_existing_anchor='lambda7,m=(1,0,0,0,0,0), universal source integer scale2, exactly frozen graded point',accepted_D_rank=101,projective_tangent_upper=24,all7_actual_source_directions_retained=True,old_actual_F_rank=48,new_joint_actual_F_rank=49,new_arc_direction={str(i)+','+str(j):str(R[i,j])for i,j in PAIRS if R[i,j]!=0},complete_actual_matrix=[[str(z)for z in row]for row in actual.tolist()],rank49_minor_rows=pr,rank49_minor_columns=pc,rank49_minor_determinant=str(minor),full350_new_DR_zero=True,actual_family_justification='Take the all-lambda literal curve and act by the inverse all6mag affine normalizer, with inverse tensor scale. Near this SAME anchor alpha*beta*I is a unit; the actual joint rational family is regular and has endpoint the whole7source. Its seven parameter derivatives plus epsilon derivative and all PGL directions are precisely this matrix. Hence rank49 is an actual single-family image lower bound, not a tangent span of unrelated families.',actual_projective_dimension_lower=24,smooth_unique_full_component=True,known_full_component_projective_dimension=24,known_component_union='MPS40 by independent literal source/SNC review',old_zero_support_goal_false=True,old_cut_reduced_support_dimension=1,old_cut_reason='The true48cut complete18fixed-variable model contains this actual transverse direction and total local dimension is24; by the implicit model and Krull cut the remaining truecut has dimension1. Exact parametrization or scheme thickness not claimed.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun);path=ROOT/'reports/CLASSIFICATION_QUADRATIC_A_SMOOTH_ANCHOR_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(smooth_unique24=True,rank_actual=49,old_zero_support_false=True,seconds=out['seconds'])))
if __name__=='__main__':main()
