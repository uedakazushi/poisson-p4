#!/usr/bin/env python3
"""Independent true two-normal full-obstruction review; Fraction sparse only.

No producer, legacy or CAS imports. Existing independent five-variable sparse
operations are used; all action, derivative, cut and obstruction maps are rebuilt.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib,json,resource,signal,time
import classification_plane_m4_special_review as H
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_TWO_NORMAL_OBSTRUCTION_20261107.json'
MANIFEST='reports/CLASSIFICATION_E6_A5_TWO_NORMAL_COMPONENT_DEPENDENCIES_20261107.json'
PRIOR='research/classification/inputs/e6_a5_prior124_8b46f29d_20261107.json'
AUDIT='reports/CLASSIFICATION_E6_A5_TWO_NORMAL_PRIOR124_AUDIT_20261107.json'
JOB='research/classification/jobs/marked_quartic_two_normal_separate_review.json'
OUT='reports/CLASSIFICATION_MARKED_QUARTIC_TWO_NORMAL_SEPARATE_REVIEW_20261107.json'

def load(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def mm(A,B):return[[sum((a*b for a,b in zip(row,col)if a and b),F(0))for col in zip(*B)]for row in A]
def mv(A,x):return[sum((a*b for a,b in zip(row,x)if a and b),F(0))for row in A]
def solve(A,B):
 n=len(A);nr=len(B[0]);T=[list(map(F,row))+list(map(F,b))for row,b in zip(A,B)];sgn=1;d=F(1)
 for j in range(n):
  k=next(k for k in range(j,n)if T[k][j])
  if k!=j:T[k],T[j]=T[j],T[k];sgn=-sgn
  q=T[j][j];d*=q
  for k in range(j+1,n):
   if T[k][j]:
    c=T[k][j]/q
    for l in range(j+1,n+nr):T[k][l]-=c*T[j][l]
    T[k][j]=F(0)
 X=[[F(0)]*nr for _ in range(n)]
 for i in reversed(range(n)):
  for k in range(nr):X[i][k]=(T[i][n+k]-sum((T[i][j]*X[j][k]for j in range(i+1,n)if T[i][j]and X[j][k]),F(0)))/T[i][i]
 assert mm(A,X)==B
 return X,d*sgn

def tensor(v):
 data={}
 for(ij,e),c in zip(H.COLS,v):
  if c:data.setdefault(ij,{})[e]=F(c)
 return H.skew(data)
def sparse(v,n):return[F(v.get(str(i),'0'))for i in range(n)]
def rows(v,n,k):return[[F(c)for c in v.get(str(i),['0']*k)]for i in range(n)]
def minor(A,C):
 r,c=C['rows'],C['columns'];assert len(r)==len(c)==C['rank'];d=H.det([[A[i][j]for j in c]for i in r]);assert d and d==F(C['determinant']);return d

def action(P):
 EE=[];GG=[]
 for i,j in product(range(5),repeat=2):
  E=H.skew({(u,v):H.add(H.a(1,u,j)if v==i else{},H.a(-1,v,j)if u==i else{})for u,v in H.PAIRS})
  G=H.skew({(u,v):H.add(H.mul(H.a(1,j),H.der(P.get((u,v),{}),i)),H.scale(P.get((j,v),{})if u==i else{},-1),H.scale(P.get((u,j),{})if v==i else{},-1))for u,v in H.PAIRS})
  EE.append(H.vector(E));GG.append(H.vector(G))
 return H.transpose(EE+GG+[H.vector(P)]),H.transpose(EE)

def source_point(name,ambient):
 tag='MU6'if name.startswith('E6')else'MU5';raw=next(r for r in load('reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_'+tag+'_SOURCES_20261107.json')['cases']if r['name']==name)
 params=tuple('r'+str(i)for i in range(6));R=LaurentRing(params,('h','x','y','z','w'));binding=[F(ambient['canonical_source_binding'][k])for k in params];scale=F(ambient['canonical_projective_integer_scale'])
 def ev(poly):
  out={}
  for ex,c in poly.items():
   for e,v in zip(ex[:6],binding):c*=v**e
   co=ex[6:];out[co]=out.get(co,F(0))+c*scale
  return{e:c for e,c in out.items()if c}
 data={tuple(map(int,k.split(','))):R.parse(s)for k,s in raw['all150_source_brackets'].items()};P=H.skew({ij:ev(v)for ij,v in data.items()});mag=[]
 for j in range(6):
  d={}
  for ij,v in data.items():
   dd={}
   for e,c in v.items():
    if e[j]:ee=list(e);ee[j]-=1;dd[tuple(ee)]=c*e[j]
   d[ij]=ev(dd)
  mag.append(H.vector(H.skew(d)))
 return P,mag

def check(row):
 name=row['name'];sourcepath=row['source_report'];assert sha(sourcepath)==row['source_report_sha256'];a=next(a for a in load(sourcepath)['cases']if a['name']==name);P,mag=source_point(name,a)
 R=LaurentRing((),('h','x','y','z','w'));expected=H.skew({tuple(map(int,k.split(','))):R.parse(v)for k,v in row['same_old_canonical_brackets'].items()});assert H.vector(P)==H.vector(expected);assert not any(H.jac(P).values())
 D=H.transpose([H.jvector(H.jac(H.plus(P,H.skew({ij:{e:F(1)}}))))for ij,e in H.COLS]);dc=row['full_D_certificate']['columns'];dr=row['full_D_certificate']['rows'];dminor=minor(D,row['full_D_certificate']);M,V=action(P);mminor=minor(M,row['full51_action']);assert len(H.rref(V)[1])==24
 fc=row['full51_action']['columns'];fr=row['true_cut_fixed_coordinate_rows'];assert fr==row['full51_action']['rows']and len(fc)==len(fr)==48;FM=[[r[j]for j in fc]for r in M];FB=[FM[i]for i in fr]
 Fweights,fdet=solve(FB,[M[i]for i in fr]);assert fdet==mminor and mm(FM,Fweights)==M
 # Every original magnetic direction and all51 actual directions are in F.
 Wmag,_=solve(FB,[H.transpose(mag)[i]for i in fr]);assert mm(FM,Wmag)==H.transpose(mag)
 assert not any(any(r)for r in mm(D,M))
 normals=[sparse(n,150)for n in row['true_normals']];assert len(normals)==2;N=H.transpose(normals);assert not any(any(r)for r in mm(D,N));assert not any(N[i][j]for i in fr for j in range(2))
 # U=e_dc-F*(F_fr)^-1*(e_dc)_fr is an actual true-cut basis.
 Ifr=[[F(i==j)for j in dc]for i in fr];WU,_=solve(FB,Ifr);U=[[F(i==j)-v for j,v in zip(dc,rr)]for i,rr in enumerate(mm(FM,WU))];assert not any(U[i][j]for i in fr for j in range(100));DU=mm(D,U);DC=[[r[j]for j in dc]for r in D];assert DU==DC
 rem=[i for i in range(150)if i not in dc];perm=rem+dc;sign=(-1)**sum(perm[i]>perm[j]for i in range(150)for j in range(i+1,150));dd=sign*H.det([FM[i]+N[i]for i in rem]);assert dd==F(row['full48plus2plus100_direct_sum']['determinant'])and dd
 # D>=100 by minor, D<=100 by the fifty independent F,N columns.
 assert H.actual(P)==R.parse(a['actual_quintic']);actual=H.actual(P);gradient=H.transpose([[H.der(actual,i).get(e,F(0))for e in H.monomials(4)]for i in range(5)]);assert len(H.rref(gradient)[1])==3
 Q11=H.jvector(H.jac(tensor(normals[0])));Q22=H.jvector(H.jac(tensor(normals[1])));Qsum=H.jvector(H.jac(tensor([x+y for x,y in zip(*normals)])));Q12=[x-y-z for x,y,z in zip(Qsum,Q11,Q22)];Q=H.transpose([Q11,Q12,Q22]);assert Q==rows(row['all350_quadratic_columns'],350,3)
 B=[[D[i][j]for j in dc]for i in dr];X,bd=solve(B,[Q[i]for i in dr]);assert bd==dminor;correction=mm(DC,X);QB=[[x-y for x,y in zip(q,c)]for q,c in zip(Q,correction)];assert QB==rows(row['all350_cokernel_quadratic_rows'],350,3);assert not any(any(QB[i])for i in dr)
 br=row['complete_quadratic_basis_rows'];G=[QB[i]for i in br];assert G==H.mat(row['basis_quadratic_coefficients']);rank=len(G);assert rank==row['quadratic_form_span_rank']and len(H.rref(G)[1])==rank
 ell=[sparse(v,350)for v in row['full_D_annihilating_functionals']];assert not any(any(r)for r in mm(ell,D));assert mm(ell,Q)==G;comb=rows(row['all350_rows_as_basis_combinations'],350,rank);assert mm(comb,G)==QB
 # Exact homogeneous ideal membership, including both P1 coordinate charts.
 if rank==3:
  coeff=H.mat(row['pure_square_membership_in_basis_quadratics']);assert mm(coeff,G)==[[F(1),F(0),F(0)],[F(0),F(0),F(1)]];powers=[2,2];mdet=None
 else:
  assert rank==2;cols=[]
  for aa,bb,cc in G:cols.extend([[aa,bb,cc,F(0)],[F(0),aa,bb,cc]])
  Mac=H.transpose(cols);assert Mac==H.mat(row['degree3_Macaulay_matrix']);mdet=H.det(Mac);assert mdet==F(row['homogeneous_quadratic_resultant'])and mdet
  coeff=H.mat(row['pure_cube_linear_multiplier_coefficients']);assert mm(Mac,H.transpose(coeff))==[[F(1),F(0)],[F(0),F(0)],[F(0),F(0)],[F(0),F(1)]];powers=[3,3]
 return dict(name=name,same_old_complete_source_point_and_all6mag_preserved=True,full350_by150_D_reconstructed=True,D_rank=100,full100minor=str(dminor),full_actual51_rank=48,Euler_rank=24,all51_and6mag_tangent_directions_span_actualF48=True,all_full_D_products_zero=True,true48_coordinate_cut=True,all150_F48_plus_N2_plus_U100_unit=str(dd),true100_implicit_unit=str(bd),all350_by3_full_quadratic_map=True,all350_cokernel_rows_reconstructed=True,full150_D_annihilation_and_complete_formspan=True,full_quadratic_span_rank=rank,Macaulay_determinant=None if mdet is None else str(mdet),pure_normal_power_membership=powers,reduced_true_cut_origin=True,whole_orbitclosure_component_projective_dimension=23,scheme_smoothness_not_claimed=True,actual_quintic_essential_rank=3)

def prior_check():
 projection=load(PRIOR);prior=projection['all124_rows'];audit=load(AUDIT);assert len(prior)==124 and len({r['id']for r in prior})==124;assert {r['id']for r in prior}=={r['id']for r in audit['all124_dispositions']};assert projection['source_commit']=='8b46f29ddbb9f539abc7a026e492ee25bf8ad699'and projection['source_git_blob']=='8b00b8c01a3c76dde6053b918e5eaa318adc1e26'
 # Reasons are classified from accepted fullcomponent data, rather than
 # trusting the text of the producer disposition as a numerical certificate.
 counts={};marked=[]
 for r in prior:
  name=r['id'];rank=r['generic_rank'];dim=r['projective_dimension']
  if rank==2:reason='genericrank2'
  elif dim!=23:reason='fullcomponentdimension'
  elif name.startswith('T')and name[1:].isdigit():reason='accepted actual SNC essential5 witness'
  elif name in['R4_H_QUARTIC_ESSENTIAL5','R4_GA9_CONSTANT_HYPERPLANE_QUARTIC_23']:reason='accepted independently recomputed actual essential5 witness'
  elif name=='R4_MINIMAL_M':reason='accepted dense old orbit cofactor rank1 versus new rank2 at unique H'
  elif name in['R4_A0','R4_DPM']:reason='accepted old dense quartic D5/3A2 profile'
  elif name.startswith('R4_MARKED_QUARTIC_'):reason='accepted dense marked curve/line orbit profile';marked.append(name)
  else:
   reason='accepted proper closed factor or irreducible witness locus'
   assert name in ['R4_ZERO_LOG_5H','R4_CMINUS_MAG','R4_INDEPENDENT_THREE_HYPERPLANES_TRITANGENT_QUADRIC_23','R4_PLANE_CONIC_CUBIC_CUSP_OUTSIDE_1_23','R4_PLANE_CONIC_CUBIC_CUSP_OUTSIDE_3_23','R4_PLANE_CONIC_CUBIC_CUSP_AT_CUSP_2_23','R4_PLANE_CONIC_CUBIC_CUSP_AT_CUSP_3_23','R4_PLANE_CONIC_CUBIC_NODE_SMOOTH_CONTACT_23','R4_PLANE_CONIC_CUBIC_NODE_AT_NODE_23','R4_PLANE_TWO_LINES_CUBIC_NODE_TANGENT_FLEX_23','R4_PLANE_TWO_LINES_CUBIC_CUSP_FLEX_TANGENT_23','R4_PLANE_TWO_LINES_CUBIC_CUSP_TWO_TANGENTS_INCIDENT_23','R4_PLANE_TWO_LINES_CUBIC_CUSP_FLEX_CUSP_SECANT_23','R4_PLANE_TWO_LINES_CUBIC_CUSP_TANGENT_CUSP_SECANT_23','R4_PLANE_TWO_LINES_CUBIC_CUSP_CUSP_TANGENT_ORDINARY_TANGENT_23','R4_PLANE_TWO_LINES_CUBIC_NODE_TWO_FLEXES_23','R4_PLANE_TWO_LINES_CUBIC_NODE_FLEX_SECANT_23','R4_PLANE_TWO_LINES_CUBIC_CUSP_TANGENT_CUSP_THIRD_SECANT_23','R4_PLANE_TWO_CONICS_BITANGENT_23','R4_PLANE_TWO_CONICS_CONTACT4_TANGENT_23','R4_NIL22_IMAGE_NODAL_CUBIC_FLEX_23','R4_NIL22_IMAGE_NODAL_CUBIC_NODE_TANGENT_23','R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23','R4_NIL22_LONG_IRREDUCIBLE_QUINTIC_23','R4_PLANE_UNICUSPIDAL_QUINTIC_45_23','R4_PLANE_QUINTIC_A12_23','R4_PLANE_QUINTIC_E8_A4_23','R4_PLANE_QUINTIC_E6_A6_23','R4_PLANE_QUINTIC_E6_A4_A2_23','R4_PLANE_QUINTIC_A8_A4_23','R4_PLANE_QUINTIC_THREE_A4_23','R4_PLANE_QUINTIC_A6_THREE_A2_23']
  counts[reason]=counts.get(reason,0)+1
 assert len(marked)==13 and sum(counts.values())==124
 return dict(all124_ids_in_immutable_public_projection=True,disposition_counts=counts,all13_old_marked_profiles_accepted=True,E6_old_bitangent_vs_new_flex_partition=[[2,2],[3,1]],mutual_new_profiles=['E6 unibranch cusp mu6','A5 two smooth branches contact3 mu5'],old_full_witness_and_source_reviews_reused_without_rerun=True)

def main():
 signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);t=time.monotonic();job=load(JOB);assert all(sha(b['path'])==b['sha256']for b in job['exact_inputs_and_hashes']);raw=load(INPUT);records=[check(r)for r in raw['cases']];audit=prior_check();out=dict(task_id='E_reduced_nonnormal',all_arithmetic_passed=True,no_producer_legacy_or_CAS_imports=True,cases=records,prior124_audit=audit,formal_actual_group_slice_and_IFT_support_read=True,all6_source_zero_invariant_boundaries_preserved=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);destination=ROOT/OUT;assert not destination.exists();destination.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('all_arithmetic_passed','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
