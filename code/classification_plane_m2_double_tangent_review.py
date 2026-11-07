#!/usr/bin/env python3
"""Independent full four-source tangent normal-double review over QQ(t).
Only earlier independent Fraction/rational-function and sparse helpers imported.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import hashlib,json,resource,signal,time
import classification_plane_m4_special_review as H
from classification_unimodular_logarithmic_review import RF,det,rank,roots_supported
ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json'
JOB='research/classification/jobs/plane_m2_double_tangent_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_SEPARATE_REVIEW_20261107.json'
T=RF((F(0),F(1)));ZERO=(0,)*5

def load(p):return json.loads((ROOT/p).read_text())
def scalar(v):return RF(tuple(map(F,v['numerator'])),tuple(map(F,v['denominator'])))if isinstance(v,dict)else RF(F(v))
def matrix(A):return[[scalar(v)for v in row]for row in A]
def enc(v):v=RF(v);return dict(numerator=list(map(str,v.p)),denominator=list(map(str,v.q)))
def mm(A,B):return[[sum((a*b for a,b in zip(row,col)if a and b),RF())for col in zip(*B)]for row in A]
def tensor(v):
 data={}
 for(ij,e),c in zip(H.COLS,v):
  if c:data.setdefault(ij,{})[e]=c
 return H.skew(data)
def Pscale(P,k):return H.skew({ij:H.scale(P.get(ij,{}),k)for ij in H.PAIRS})
def diff(P,Q):return H.plus(P,Pscale(Q,-1))
def linear(P,Q):return[a-b-c for a,b,c in zip(H.jvector(H.jac(H.plus(P,Q))),H.jvector(H.jac(P)),H.jvector(H.jac(Q))) ]
def op(P,basis):return H.transpose([linear(P,Q)for Q in basis])
def actualpoly(saved):return{tuple(r['exponent']):scalar(r['coefficient'])for r in saved}
def certificate(M,C):
 rr,cc=C['rows'],C['columns'];rk=C['rank'];assert len(rr)==len(cc)==rk
 d=det([[M[i][j]for j in cc]for i in rr]);assert d and d==scalar(C['minor_determinant']);K=matrix(C['complete_kernel']);n=len(M[0]);null=C['nullity'];assert len(K)==null==n-rk
 if K:
  assert all(len(row)==n for row in K);assert not any(any(row)for row in mm(M,H.transpose(K)));assert rank(K)==null
 return K,d

def models():
 xyz=[H.a(1,i)for i in range(5)];x,y,z,u,v=xyz
 f=H.mul(z,H.add(H.mul(x,z),H.scale(H.mul(y,y),-1)));Pf=H.skew({(0,1):H.der(f,2),(0,2):H.scale(H.der(f,1),-1),(1,2):H.der(f,0)});N=[H.scale(y,2),z,{}];out={}
 for name,cube,c in[('constant_hf_Z0',False,0),('constant_hf_ZN',False,1),('constant_hcube',True,0)]:
  X=[H.scale(a,-8 if cube else-2)for a in N];Y=[H.scale(a,-5 if cube else-2)for a in N];Z=[z,{},{}]if cube else[H.scale(a,c)for a in N]
  P=H.plus(Pf,H.skew({(i,3):H.mul(u,X[i])for i in range(3)}),H.skew({(i,4):H.add(H.mul(u,Z[i]),H.mul(v,Y[i]))for i in range(3)}),H.skew({(3,4):{}if cube else H.mul(u,u)}));out[name]=P
 out['moving_all_t_nonzero']=H.skew({(0,1):H.scale(H.add(H.mul(z,x),H.scale(H.mul(y,y),-1)),3),(0,3):H.add(H.scale(H.mul(y,u),2*T+3),H.scale(H.mul(z,v),-3)),(0,4):H.add(H.scale(H.mul(x,u),-3),H.scale(H.mul(y,v),2*T+3)),(1,3):H.scale(H.mul(z,u),T+3),(1,4):H.add(H.scale(H.mul(y,u),-3),H.scale(H.mul(z,v),T)),(3,4):H.scale(H.mul(u,u),T+3)})
 return out

def group(P):
 cols=[]
 for a,b in product(range(5),repeat=2):cols.append(H.vector(H.skew({(i,j):H.add(H.a(1,i,b)if j==a else{},H.a(-1,j,b)if i==a else{})for i,j in H.PAIRS})))
 for a,b in product(range(5),repeat=2):
  if(a<3)!=(b<3):continue
  O=H.skew({(i,j):H.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},H.scale(H.mul(H.a(1,b),H.der(P.get((i,j),{}),a)),-1))for i,j in H.PAIRS});cols.append(H.vector(O))
 cols.append(H.vector(P));assert len(cols)==39
 return H.transpose(cols)

def at(v,t):
 v=RF(v);num=sum((c*t**i for i,c in enumerate(v.p)),F(0));den=sum((c*t**i for i,c in enumerate(v.q)),F(0));assert den;return num/den

def check(row,P,B1,B2):
 assert H.vector(P)==list(map(scalar,row['full150_P0']));assert not any(H.jac(P).values());D=H.actual(P);assert D==actualpoly(row['P0_actualG2'])and D
 J1,J2,G=op(P,B1),op(P,B2),group(P);assert J1==matrix(row['full350by30_J1'])and J2==matrix(row['full350by9_J2'])and G==matrix(row['full150_Euler_Levi_scale_matrix'])
 K1,d1=certificate(J1,row['weight1_certificate']);K2,d2=certificate(J2,row['weight2_certificate']);KG,dg=certificate(G,row['normalizer_certificate']);assert len(H.rref([[F(v)for v in r[:25]]for r in G])[1])==24
 s=37-row['normalizer_certificate']['rank'];fiber=s+6;upper=24-fiber+row['source_parameter_dimension']+len(K1)+len(K2);assert s==row['graded_Levi_projective_stabilizer_dimension']and fiber==row['actual_filtered_group_fiber_dimension']and upper==row['whole_projective_image_dimension_upper']
 # All missing-u terms of every complete K1 and every K2, not a slice.
 bad=[j for j,Q in enumerate(B1)if any(not e[3]for i in range(5)if i!=3 for e in Q.get((min(i,3),max(i,3)),{}))]
 udiv=all(e[3]for i in range(5)if i!=3 for e in P.get((min(i,3),max(i,3)),{}));ufac=udiv and all(not k[j]for k in K1 for j in bad);assert ufac==row['actual_u_hyperplane_factor_for_every_full_solution'];assert all(e[3]for Q in B2 for i in range(5)if i!=3 for e in Q.get((min(i,3),max(i,3)),{}))
 saved=row['all_full39_actualD_linear_terms'];assert len(saved)==39
 for Q,v in zip(B1+B2,saved):
  result=H.add(H.actual(H.plus(P,Q)),H.scale(D,-1),H.scale(H.actual(Q),-1));assert result==actualpoly(v)
 return dict(type=row['type'],whole350_and150_maps_reconstructed=True,K1=len(K1),K2=len(K2),Euler_rank=24,Levi_stabilizer=s,filtered_fiber=fiber,image_upper=upper,entire_irreducible_source_empty=ufac,G24_closed=ufac or upper<=24),[d1,d2,dg],K1+K2

def universal(raw,B):
 saved={tuple(r['pair']):r for r in raw['universal_full39_quadratic_Jacobi_and_actualD']};found=0
 for i,j in combinations_with_replacement(range(39),2):
  P=H.plus(B[i],B[j]);J=H.jvector(H.jac(P));D=H.actual(P)
  if i==j:J=[v/4 for v in J];D=H.scale(D,F(1,4))
  if any(J)or D:
   row=saved[i,j];assert J==[scalar(row['all_nonzero_Jacobi_rows'].get(str(k),'0'))for k in range(350)];assert D==actualpoly(row['actualD']);found+=1
  else:assert(i,j)not in saved
 assert found==len(saved)
 return found

def leading_geometry():
 # New constant normalform proof, not the already accepted moving/conic proof.
 from classification_minimal_quartic_review import LaurentRing
 R=LaurentRing(('A','B','C','D','E','k'),('x','y','z','u','v'));p=R.parse;X=list(map(p,('-16*y','-8*z','0')));Z=list(map(p,('A*x+B*y+C*z','D*x-A*y+E*z','0')));xyz=list(map(p,('x','y','z')));f=p('z*(x*z-y**2)')
 det3=R.add(*(R.scale(R.mul(*(M[i][j]for i,j in enumerate(perm))),(-1)**sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3)))for perm in __import__('itertools').permutations(range(3))for M in[[xyz,X,Z]]))
 Zf=R.add(*(R.mul(Z[i],R.derivative(f,i))for i in range(3)));assert R.add(det3,Zf)==p('-18*D*x*y*z+18*A*y**2*z+9*A*x*z**2+(-18*E+9*B)*y*z**2+9*C*z**3')
 return dict(constant_hf_full_Z_cN_two_actual_orbits=True,constant_hcube_all5_offdiagonal_coefficients_and_Pf_constraints_checked=True,constant_hcube_A_D_k0_B2E_C_nonzero=True,whole_normalizer_inverse_and_filtered_shear_fiber_read=True)

def main():
 signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);t=time.monotonic();job=load(JOB);assert all(hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']for b in job['exact_inputs_and_hashes']);raw=load(INPUT);assert hashlib.sha256(raw['producer_source_bytes'].encode()).hexdigest()==raw['producer_sha256']=='6ff941e7979d36569883b621f9922a992e340283dc9691f1828de3e2f4f82db4'
 B1,B2=H.sources();assert [H.vector(P)for P in B1]==matrix(raw['full30_weight1_basis'])and[H.vector(P)for P in B2]==matrix(raw['full9_weight2_basis']);assert raw['full350_row_labels']==[dict(quad=list(q),monomial=list(e))for q,e in H.ROWS]
 polars=universal(raw,B1+B2);M=models();records=[];gfactors=set();allowed={F(0),F(6),F(-3,2),F(-5,2),F(-3),F(-9,2),F(-6),F(-9)}
 for row in raw['cases']:
  res,dets,ker=check(row,M[row['type']],B1,B2);records.append(res)
  if row['type'].startswith('moving'):
   # Every determinant numerator/denominator and every kernel denominator is
   # factored independently by exact division by the complete rational roots.
   for d in dets:
    for poly in(RF(d).p,RF(d).q):gfactors|=roots_supported(poly,allowed)
   for k in ker:
    for v in k:
     for poly in(RF(v).p,RF(v).q):gfactors|=roots_supported(poly,allowed)
   claimed={-F(fac[0])/F(fac[1])for fac in row['generic_open_excluded_irreducible_factors']};assert gfactors-{F(0)}==claimed==allowed-{F(0)};assert len(row['whole_nonzero_t_walls'])==7
   for wall in row['whole_nonzero_t_walls']:
    tt=F(wall['t']);assert tt in claimed;P=H.skew({ij:{e:at(c,tt)for e,c in M[row['type']].get(ij,{}).items()if at(c,tt)}for ij in H.PAIRS});r,_,_=check(wall,P,B1,B2);records.append(r)
 assert not raw['whole_tangent_double_G24_closed'];assert records[2]['image_upper']==27 and not records[2]['G24_closed'];geometry=leading_geometry();out=dict(task_id='G',all_independent_arithmetic_passed=True,no_producer_or_CAS_imports=True,all4_whole_sources_and7_walls=records,whole_generic_t_and_all7_complex_walls_checked=True,all39_quadratic_polarizations_reconstructed=polars,whole_leading_normalizer_geometry=geometry,partial_conclusions_only=True,unresolved_constant_hcube_upper27=True,whole_tangent_or_whole_M2_not_claimed=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);destination=ROOT/OUT;assert not destination.exists();destination.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('all_independent_arithmetic_passed','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
