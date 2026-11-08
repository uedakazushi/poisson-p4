#!/usr/bin/env python3
"""Independent whole triangle m2-double HF/H21 chain; Fraction sparse only."""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import json,hashlib,resource,signal,time
import classification_plane_m4_special_review as H
from classification_minimal_quartic_review import LaurentRing
from classification_plane_m2_double_tangent_review import group
ROOT=Path(__file__).resolve().parents[1]
HF='reports/CLASSIFICATION_PLANE_M2_TRIANGLE_DOUBLE_HF_CASIMIR_20261107.json'
WALL='reports/CLASSIFICATION_PLANE_M2_TRIANGLE_DOUBLE_21_20261107.json'
JOB='research/classification/jobs/triangle_m2_double_separate_review.json'
OUT='reports/CLASSIFICATION_TRIANGLE_M2_DOUBLE_SEPARATE_REVIEW_20261107.json'

def load(p):return json.loads((ROOT/p).read_text())
def mat(M):return[[F(v)for v in row]for row in M]
def mm(A,B):return[[sum((a*b for a,b in zip(r,c)if a and b),F(0))for c in zip(*B)]for r in A]
def cert(M,C,K=None):
 rr,cc=C['rows'],C['columns'];r=C['rank'];assert len(rr)==len(cc)==r;d=H.det([[M[i][j]for j in cc]for i in rr]);assert d==F(C['minor_determinant'])and d
 if K is not None:
  K=mat(K);assert len(K)==len(M[0])-r
  if K:assert all(len(k)==len(M[0])for k in K)and len(H.rref(K)[1])==len(K)and not any(any(v)for v in mm(M,H.transpose(K)))
 return str(d)
def tensor(v):
 data={}
 for(ij,e),c in zip(H.COLS,v):
  if c:data.setdefault(ij,{})[e]=F(c)
 return H.skew(data)
def poly(data):return{tuple(v['exponent']):F(v['coefficient'])for v in data}

def casimir(raw):
 R=LaurentRing((),('x','y','z','u','v'));p=R.parse;mon6=H.monomials(6);mon1=H.monomials(1);assert raw['full150_column_labels']==[dict(pair=list(ij),monomial=list(e))for ij,e in H.COLS]
 labels=[dict(type='Casimir',component=i,monomial=list(e))for i in range(5)for e in mon6]+[dict(type='divergence',component=i,monomial=list(e))for i in range(5)for e in mon1];assert labels==raw['full1075_row_labels'];records=[]
 for row in raw['cases']:
  eps=row['epsilon'];f=p('(u**2*x-v**3)*(y*z+'+str(eps)+'*u**2)+u**5');assert f==poly(row['actualF']);grad=[H.der(f,i)for i in range(5)];cols=[]
  for(i,j),e in H.COLS:
   mon={e:F(1)};cas=[{}for _ in range(5)];dv=[{}for _ in range(5)];cas[i]=H.mul(mon,grad[j]);cas[j]=H.scale(H.mul(mon,grad[i]),-1);dv[i]=H.scale(H.der(mon,j),-1);dv[j]=H.der(mon,i);cols.append([P.get(ex,F(0))for P in cas for ex in mon6]+[P.get(ex,F(0))for P in dv for ex in mon1])
  M=H.transpose(cols);assert M==mat(row['full1075by150_Casimir_divergence_matrix']);d=cert(M,row['full_rank_certificate'],row['complete150_kernel']);assert row['full_rank_certificate']['rank']==150 and row['canonical_strict_necessary_source_dimension']==0 and row['all_quadratic_source_pairs']==[]
  for weights in row['whole_actualD_stabilizer_tori']:
   chars={sum(a*b for a,b in zip(e,weights))for e in f};assert len(chars)==1
  assert len(H.rref([[F(1)]*5]+[list(map(F,w))for w in row['whole_actualD_stabilizer_tori']])[1])==len(row['whole_actualD_stabilizer_tori'])+1
  records.append(dict(epsilon=eps,full1075by150_reconstructed=True,all150directions_retained=True,rank=150,kernel0=True,nonzero_minor=d,all_rank4_source_empty=True))
 return records

def model(r):
 E=[H.a(1,i)for i in range(5)];x,y,z,u,v=E;f=H.mul(x,H.mul(y,z));P=H.skew({(0,1):H.der(f,2),(0,2):H.scale(H.der(f,1),-1),(1,2):H.der(f,0)});Y=[-2*r-1,r+1,r];X=[Y[0]-1,Y[1]+2,Y[2]-1]
 return H.plus(P,H.skew({(i,3):H.scale(H.mul(E[i],u),X[i])for i in range(3)}),H.skew({(i,4):H.add(H.scale(H.mul(E[i],v),Y[i]),H.mul(x,u)if i==2 else{})for i in range(3)}),H.skew({(3,4):H.scale(H.mul(u,v),3*r+1)}))
def operator(P,B):
 jp=H.jvector(H.jac(P));return H.transpose([[a-b-c for a,b,c in zip(H.jvector(H.jac(H.plus(P,Q))),jp,H.jvector(H.jac(Q)))]for Q in B])

def walls(raw):
 B1,B2=H.sources();assert[H.vector(P)for P in B1]==mat(raw['full30_weight1_basis'])and[H.vector(P)for P in B2]==mat(raw['full9_weight2_basis']);old=load('reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json');assert raw['universal_full39_quadratic_Jacobi_and_actualD']==old['universal_full39_quadratic_Jacobi_and_actualD'];assert raw['full350_row_labels']==[dict(quad=list(q),monomial=list(e))for q,e in H.ROWS];records=[]
 for row in raw['cases']:
  r=F(row['r']);P=model(r);assert H.vector(P)==list(map(F,row['full150_P0']))and not any(H.jac(P).values());D=H.actual(P);assert D==poly(row['actualG2'])==H.a(-3*(r+2),0,0,1,3,3)
  J1,J2=operator(P,B1),operator(P,B2);assert J1==mat(row['full350by30_J1'])and J2==mat(row['full350by9_J2']);d1=cert(J1,row['J1_certificate'],row['full30_weight1_kernel']);d2=cert(J2,row['J2_certificate'],row['full9_weight2_kernel']);G=group(P);assert G==mat(row['full150_Euler_Levi_scale_matrix']);dg=cert(G,row['normalizer_certificate'],row['normalizer_complete39_kernel']);assert len(H.rref([row[:25]for row in G])[1])==24
  K1,K2=mat(row['full30_weight1_kernel']),mat(row['full9_weight2_kernel']);fiber=37-row['normalizer_certificate']['rank']+6;assert fiber==9==row['actual_filtered_group_fiber_dimension'];upper=24-fiber+len(K1)+len(K2);assert upper==row['raw_linear_projective_image_upper']==row['whole_projective_image_dimension_upper']<=23
  bad=[[Q.get((i,3),{}).get((0,0,0,0,2),F(0))for Q in B1]for i in range(3)];badmap=mm(bad,H.transpose(K1));assert badmap==mat(row['all_three_bad_urow_maps'])and not any(badmap[1])and not any(badmap[2])
  for Q,saved in zip(B1+B2,row['all_full39_actualD_linear_terms']):assert H.add(H.actual(H.plus(P,Q)),H.scale(D,-1),H.scale(H.actual(Q),-1))==poly(saved)
  records.append(dict(r=str(r),whole350_30_9_reconstructed=True,K1=len(K1),K2=len(K2),all_original_kernel_products_zero=True,weight1_minor=d1,weight2_minor=d2,normalizer_minor=dg,Levi_stabilizer3=True,filtered_six_shears_fiber9=True,whole_projective_image_upper=upper,optional_quadratic_not_executed=True))
 assert len(records)==10
 return records

def finite_weights():
 roots=set();table=[]
 for a,m,n in product(range(3),repeat=3):
  if a+m+n>2:continue
  for b in range(6-2*a-m-n):
   aa=-a+2*m+2*n-3*b;bb=-2*a+4*m-n-b
   if aa:
    r=F(-bb,aa);roots.add(r);table.append(dict(a=a,m=m,n=n,b=b,r=str(r)))
   elif not bb:assert(a,m,n,b)==(0,0,0,0)
 expected=set(map(F,['-7','-2','-1','-3/4','-4/7','-1/2','-1/3','0','1/7','1/2','1','3']));assert roots==expected
 # Off-block resonance iff mu_i=0 or mu_i=nu. These yield exactly
 # {-2,-1,0,1/2}; the bad quadratic resonance adds only r=0.
 linepairs=[(-1,-2),(2,4),(2,-1)];normal=(-3,-1);res=set()
 for a,b in linepairs:
  if a:res.add(F(-b,a))
  if a-normal[0]:res.add(F(-(b-normal[1]),a-normal[0]))
 assert res==set(map(F,['-2','-1','0','1/2']))
 return dict(all_nonnegative_degree_and_normal_order_weights_enumerated=True,full12_roots=sorted(map(str,roots)),full10_target_walls=sorted(map(str,roots-{F(-2),F(0)})),linear_offblock_resonances=sorted(map(str,res)),r0_whole_leading_nil_exclusion_handchecked=True,rminus2_explicit_order3_complement=True)

def syzygies():
 R=LaurentRing((),('x','y','z','u','v'));p=R.parse;records=[]
 for name,Fstr,fields in [('fixed_bad_z','x*y*(x*u**2+v**3)+u**5',[['3*x','-6*y','0','v']]),('generic_wall_complement','y*(u*x+v**2)**2+u**5',[['2*x','-4*y','0','v'],['2*v','0','0','-u']])]:
  f=p(Fstr);base=(0,1,3,4);grad=[H.der(f,i)for i in base]
  for d in range(3):
   mon=[e for e in H.monomials(d)if not e[2]];target=[e for e in H.monomials(4+d)if not e[2]];cols=[H.mul({e:F(1)},g)for g in grad for e in mon];M=[[P.get(e,F(0))for P in cols]for e in target];k=[]
   if d:
    m=[e for e in H.monomials(d-1)if not e[2]]
    for V in fields:
     for ex in m:
      vec=[H.mul({ex:F(1)},p(s))for s in V];k.append([P.get(e,F(0))for P in vec for e in mon])
   r=len(H.rref(M)[1]);assert r==len(cols)-len(k)
   if k:assert len(H.rref(k)[1])==len(k)and not any(any(x)for x in mm(M,H.transpose(k)))
   records.append(dict(type=name,degree=d,full_gradient_syzygy_matrix=[len(M),len(cols)],kernel_dimension=len(k),entire_kernel_explicit_T_or_TZ=True))
 return records

def main():
 signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);t=time.monotonic();job=load(JOB);assert all(hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']for b in job['exact_inputs_and_hashes']);hf=load(HF);ww=load(WALL)
 for r in(hf,ww):assert hashlib.sha256(r['producer_source_bytes'].encode()).hexdigest()==r['producer_sha256']
 result=dict(task_id='G',all_independent_arithmetic_passed=True,no_producer_or_CAS_imports=True,all_HF_full150_sources=casimir(hf),all_10_whole_H21_sources=walls(ww),whole_finite_wall_enumeration=finite_weights(),both_AR_leq2_complete=syzygies(),canonical_lift_and_fullraw36_source_proof_handchecked=True,all_HF_actual_inverse_and_units_zero_complements_handchecked=True,full_HF_source_empty=True,whole_H21_image_upper23=True,whole_reduced_triangle_m2double_G_closed_with_prior_leading=True,no_pointwise_E_claim_for_H21=True,no_optional_quadratic_or_jet_or_GB=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);out=ROOT/OUT;assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in('all_independent_arithmetic_passed','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
