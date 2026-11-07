#!/usr/bin/env python3
"""Two whole9-parameter actual limits into the accepted nil31 half/third25.

Verify the exact punctured GL/source identities and all polynomial endpoints;
the Laurent parameter discovery is not needed to reproduce this certificate.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_rank3_mixed_component_review import PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]

def main():
 begun=time.monotonic();h,q1,q2,x,y,e=s.symbols('h q1 q2 x y e');xx=(h,q1,q2,x,y);aff=xx[1:]
 mm=s.symbols('m0:6');n,c,d=s.symbols('n c d');params=mm+(n,c,d);env={str(z):z for z in xx+params};env.update(dict(zip(['x1','x2','x3','x4'],aff)))
 nil22_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';nil31_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json';original_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
 nil22=json.loads(nil22_path.read_text())['cases'];nil31=json.loads(nil31_path.read_text())['branches'];original=json.loads(original_path.read_text())['cases'][0]
 def decode(data,localenv=env):
  P=s.zeros(5)
  for ij,f in data.items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=localenv);P[j,i]=-P[i,j]
  return P
 def encode(P):return {str(i)+','+str(j):str(s.factor(P[i,j]))for i,j in itertools.combinations(range(P.rows),2)if P[i,j]!=0}
 def sparse(P):return {ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())for ij in PAIRS if P[ij]!=0}
 def pf(P):return s.factor(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
 M0,M1,M2,M3,M4,M5=mm;records=[]
 for label,index in [('quadratic_B_to_half',8),('triple_C_to_third',6)]:
  target_h=decode(nil22[index]['full_six_magnetic_source']);oldtarget=target_h
  if label.startswith('triple'):
   perm=(3,1,2,0,4);target_h=target_h.extract(perm,perm).subs({h:x,x:h},simultaneous=True)
   R=s.Matrix([0]+[s.cancel(target_h[0,j]/h)for j in range(1,5)])
   assert all(not s.denom(z).free_symbols.intersection(xx)for z in R)
   target_h=(target_h-s.Matrix(xx)*R.T+R*s.Matrix(xx).T).applyfunc(s.expand)
   infinity=dict(permutation=[3,1,2,0,4],Euler_gauge_vector=list(map(str,R)),all_identity_checked=True)
  else:infinity=None
  assert target_h.row(0)==s.zeros(1,5);target=target_h[1:,1:].subs(h,1)
  if label.startswith('quadratic'):
   X,Y,Z,W=s.symbols('X Y Z W');old=(X,Y,Z,W);ns=s.symbols('n0:15');local={str(z):z for z in ns};local.update(dict(zip(['x1','x2','x3','x4'],old)))
   source=s.zeros(4);bind={ns[6]:0,ns[7]:0,ns[9]:0,ns[10]:-s.Rational(1,2),ns[11]:0}
   for (i,j),f in zip(itertools.combinations(range(4),2),original['full_six_brackets']):source[i,j]=s.sympify(f,locals=local).subs(bind);source[j,i]=-source[i,j]
   periods={ns[0]:-(e*M0+1)/e**2,ns[1]:-M3/e,ns[2]:(c-e*M4)/e**2,ns[3]:(6*d-4*e*M5+1)/(4*e**2),ns[4]:-M1/e,ns[5]:-1/e,ns[8]:-(e*M2-n)/e**2,ns[12]:n,ns[13]:c,ns[14]:d}
   B=s.Matrix([[0,1,0,-1/e],[e,0,0,0],[0,0,1,0],[0,0,0,1]]);sub=dict(zip(old,[q2/e,q1+y/e,x,y]));expected_pf=-x*x/2;expected_mod=s.Matrix([x,y+e*q1,0,0]);component='accepted nil31-image half25 component'
   source_recognition='Original complete15 source with n6=n7=n9=n11=0,n10=-1/2 and n5=-1/e!=0: exactly the accepted whole10-parameter half branch, with original h as its intrinsic triplehyperplane.'
  else:
   a,b,cc,dd=s.symbols('a b cc dd');old=(a,b,cc,dd);AA=s.symbols('A0:6');u,v,t,q=s.symbols('u v t q');local={str(z):z for z in AA+(u,v,t,q)};local.update(dict(zip(['h','x1','x2','x3','x4'],[h,a,b,cc,dd])))
   source=decode(nil31[1]['whole_homogeneous_h_Casimir_brackets'],local)[1:,1:].subs(h,1)
   periods={AA[0]:-M3/e,AA[1]:-(e*M1+3)/e**2,AA[2]:(c-e*M4)/e**2,AA[3]:-M0/e,AA[4]:-(e*M2-n)/e**2,AA[5]:(9*d-6*e*M5+2)/(6*e**2),u:n,v:c,t:d,q:-3/e}
   B=s.Matrix([[0,0,1,-1/e],[0,e,0,0],[1,0,0,0],[0,0,0,1]]);sub=dict(zip(old,[x,q2/e,q1+y/e,y]));expected_pf=x*x;expected_mod=-s.Rational(2,3)*s.Matrix([1,y+e*q1,0,0]);component='accepted nil31-image third25 component'
   source_recognition='Normalized whole10-parameter third branch of NIL31_IMAGE_ZERO_COMPONENTS, with q=-3/e!=0.'
  assert B.det()==-e
  family=(B*source*B.T).subs(sub,simultaneous=True).subs(periods,simultaneous=True).applyfunc(s.cancel)
  assert (family.subs(e,0)-target).applyfunc(s.expand)==s.zeros(4)
  for z in family:assert s.Poly(z,*aff,e,*params).total_degree()>=0;assert s.Poly(z,*aff).total_degree()<=2
  strict=[s.expand(sum(family[i,l]*s.diff(family[j,k],aff[l])+family[j,l]*s.diff(family[k,i],aff[l])+family[k,l]*s.diff(family[i,j],aff[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)]
  assert strict==[0]*4;assert s.expand(pf(family)-expected_pf)==0
  modular=s.Matrix([s.expand(sum(s.diff(family[i,j],aff[i])for i in range(4)))for j in range(4)])
  assert (modular-expected_mod).applyfunc(s.expand)==s.zeros(4,1)
  homogeneous=s.zeros(5)
  for i,j in itertools.combinations(range(4),2):homogeneous[i+1,j+1]=s.expand(h*h*family[i,j].subs({z:z/h for z in aff},simultaneous=True));homogeneous[j+1,i+1]=-homogeneous[i+1,j+1]
  assert jacobi(sparse(homogeneous)).applyfunc(s.expand)==s.zeros(350,1)
  records.append(dict(name=label,all_nine_original_parameters=list(map(str,params)),original_source_index=index,source_component=component,source_recognition=source_recognition,
   changed_infinity=infinity,actual_GL_matrix=[[str(z)for z in row]for row in B.tolist()],GL_determinant=str(B.det()),source_Laurent_parameters={str(z):str(s.factor(f))for z,f in periods.items()},
   whole_regular_affine_family=encode(family),whole_regular_homogeneous_family=encode(homogeneous),whole_original_endpoint_identity=True,
   complete_strict_Jacobi_and_projective350_zero=True,whole_family_degree_at_most_two=True,actual_affine_Pfaffian=str(expected_pf),whole_modular=list(map(str,modular)),
   whole_source_covered_by_accepted_component=True,logic='For each fixed choice of all9 original parameters, every e!=0 tensor is an actual GL transform of the accepted half/third fullsource and has rank4. The exact tensor is polynomial in e with original source endpoint at0. Closedness gives membership at0, and the projective infinity permutation/Euler identity returns the entire original source. No uniqueness or matching-divisor inference is used.'))
 inputs=[nil22_path,nil31_path,original_path,Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py']
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),cases=records,new_component_count=0,independent_review_required=True,seconds=time.monotonic()-begun)
 dest=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_COMPONENT_LIMITS_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,whole_sources_covered=[r['name']for r in records],seconds=result['seconds'])))
if __name__=='__main__':main()
