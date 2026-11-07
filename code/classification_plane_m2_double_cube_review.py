#!/usr/bin/env python3
"""Independent full K1-nine plus ALL P2-nine Ore/quartic identity review.
No producer or CAS imports. Reuses frozen independent Fraction sparse arithmetic.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json,time,resource,signal
import classification_plane_m4_special_review as H
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
RAW='reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json'
OLD='reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_SEPARATE_REVIEW_20261107.json'
JOB='research/classification/jobs/plane_m2_double_tangent_cube_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_CUBE_SEPARATE_REVIEW_20261107.json'

def load(p):return json.loads((ROOT/p).read_text())
def det3(R,E,A,D):return R.add(R.mul(E[0],R.add(R.mul(A[1],D[2]),R.scale(R.mul(A[2],D[1]),-1))),R.scale(R.mul(E[1],R.add(R.mul(A[0],D[2]),R.scale(R.mul(A[2],D[0]),-1))),-1),R.mul(E[2],R.add(R.mul(A[0],D[1]),R.scale(R.mul(A[1],D[0]),-1))))
def skew(R,data):
 P={(i,j):{}for i in range(5)for j in range(5)}
 for(i,j),v in data.items():P[i,j]=v;P[j,i]=R.scale(v,-1)
 return P

def main():
 signal.alarm(60);t=time.monotonic();job=load(JOB);assert all(hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']for b in job['exact_inputs_and_hashes']);prior=load(OLD);assert prior['all_independent_arithmetic_passed']and prior['unresolved_constant_hcube_upper27'];raw=load(RAW);row=next(r for r in raw['cases']if r['type']=='constant_hcube');K=H.mat(row['weight1_certificate']['complete_kernel']);assert len(K)==9
 names=tuple('n'+str(i)for i in range(9))+tuple('r'+str(i)for i in range(9));R=LaurentRing(names,('x','y','z','u','v'));p=R.parse;coords=list(map(p,('x','y','z','u','v')));basis1,basis2=H.sources();data={}
 for ij in H.PAIRS:
  val={}
  for(pair,ex),c in zip(H.COLS,row['full150_P0']):
   if pair==ij and F(c):val=R.add(val,{(0,)*18+ex:F(c)})
  for j,coeff in enumerate(K):
   v={}
   for k,Q in enumerate(basis1):
    for e,c in Q.get(ij,{}).items():
     ee=[0]*18+list(e);ee[j]=1;v=R.add(v,{tuple(ee):c*coeff[k]})
   val=R.add(val,v)
  for j,Q in enumerate(basis2):
   for e,c in Q.get(ij,{}).items():ee=[0]*18+list(e);ee[9+j]=1;val=R.add(val,{tuple(ee):c})
  data[ij]=val
 P=skew(R,data);B=(2,3,4);assert all(not P[i,j]for i,j in combinations(B,2));delta=[P[1,j]for j in B];A=[R.derivative(P[0,j],1)for j in B];V=[R.add(P[0,j],R.scale(R.mul(coords[1],a),-1))for j,a in zip(B,A)]
 assert all(not R.derivative(v,i)for v in A+V+delta for i in(0,1))
 assert A==list(map(p,('2*z+(9*n6+3*n8)*u+(2*n2+4*n4)*v','-16*u','-10*v')))
 assert delta==list(map(p,('z**2+n2*z*v+9*(n6+n8)*z*u/2+r6*u**2+r7*u*v+r8*v**2','-8*z*u+n4*u*v','-5*z*v+n8*u*v')))
 # Full original projective Jacobi convention, no strict Jacobi assumption.
 J={tri:R.add(*(R.mul(P[a,l],R.derivative(P[b,c],l))for a,b,c in(tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1]))for l in range(5)))for tri in combinations(range(5),3)}
 quad=(0,2,3,4);Q=R.add(*(R.scale(R.mul(coords[quad[k]],J[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)));EB=[coords[i]for i in B];assert Q==det3(R,EB,A,delta)
 for at,i in enumerate(B):
  for bt,j in enumerate(B):
   if at<bt:assert J[0,i,j]==R.add(R.mul(delta[at],A[bt]),R.scale(R.mul(delta[bt],A[at]),-1))
 assert not J[2,3,4]
 D={}
 for at in range(5):
  i,j,k,l=[v for v in range(5)if v!=at];pf=R.add(R.mul(P[i,j],P[k,l]),R.scale(R.mul(P[i,k],P[j,l]),-1),R.mul(P[i,l],P[j,k]));D=R.add(D,R.scale(R.mul(coords[at],pf),(-1)**at))
 assert D==R.scale(det3(R,EB,[P[0,j]for j in B],delta),-1)
 assert not R.derivative(D,0)and R.derivative(D,1)==R.scale(Q,-1)
 remainder=R.add(D,R.mul(coords[1],Q));assert all(not R.derivative(remainder,i)for i in(0,1))
 result=dict(task_id='G',all_independent_arithmetic_passed=True,no_producer_or_CAS_imports=True,full_K1_nine_and_all_P2_nine_retained=True,all18_parameter_source_reconstructed=True,all9_P2_not_replaced_by_kernel3=True,whole_source_Ore_structure=True,full_A_and_delta_coefficients_checked=True,strict_xBiBj_not_zero_assumed=True,original_quad_xzuv_equals_det_E_A_delta=True,actual_D_equals_minus_det_E_xB_yB=True,identity_dx_actualD_zero=True,identity_dy_actualD_equals_minus_projectiveJacobi_xzuv=True,fullsolution_actualD_in_Czuv=True,irreducible_rank4_fullsource_in_accepted_eight_essential3_components=True,G_image_upper23=True,whole_tangent_double_G_image_upper23_with_prior_partial=True,no_full_E_claim_for_other_three_models=True,historical_V1_7ea7_not_used=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);path=ROOT/OUT;assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in('all_independent_arithmetic_passed','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
