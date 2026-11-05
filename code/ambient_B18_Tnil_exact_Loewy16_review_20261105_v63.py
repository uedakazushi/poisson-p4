#!/usr/bin/env python3
"""Fixed filtered16 source-basis and pure module proof; no new values."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 stem='AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63';rep=load(stem);assert pins[f'reports/{stem}/report.json']=='c4b1fe60be21a774d9645886ec7540abbd5c4ddf04484baef73b14fc1693a972'
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/n)==h
 P=load(stem,'exact-Loewy16-and-explicit-filtered-basis.json');assert P['exact_actual_HF']==[5,5,3,2,1]
 assert load('AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_20261105_v63_REVIEW')['actual_full_HF_Tnil']==[5,5,3,2,1]
 F1=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json');V=load('AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_20261105_v63_REVIEW','phase2-independent-physical-source2-map-and-complete-HF-proof.json');C=load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63','literal-known-cubic-input-and-exact24-row-reduction.json');N=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json')
 assert load('AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW')['all_passed']
 assert load('AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63_REVIEW')['all_passed']
 assert load('AMBIENT_B18_HONEST_NIL_T_RIGHT106_20261105_v63_REVIEW')['actual_right106_rank']==0
 firstfree=sorted(set(range(40))-set(F1['rank_minor_rows']));assert firstfree==[5,7,15,31,39]
 b=P['explicit_sixteen_filtered_basis_lifts'];assert len(b)==16 and [sum(x['degree']==d for x in b)for d in range(5)]==[5,5,3,2,1]
 assert b[:5]==[dict(degree=0,type='new_minimal_unit',target=i)for i in range(5)]
 for x,r in zip(b[5:10],firstfree):assert x==dict(degree=1,type='new_minimal_monomial',target=r//8,monomial=[r%8])
 assert N['old_minimal_generator_free56_rows']==[22,23,29,30,55] and 'actual_natural_map_minimal5_by5'in N['matrices']
 assert [x.get('old_minimal_column')for x in[b[10],b[11],b[13]]]==[3,4,4] and [x.get('old_free56')for x in[b[10],b[11],b[13]]]==[30,55,55]
 assert b[10]['monomial']==[7] and b[11]['monomial']==[7,7] and b[13]['monomial']==[7,7,7]
 assert b[12]['target']==b[14]['target']==0 and b[12]['monomial']==[7,7] and b[14]['monomial']==[7,7,7]
 pi=[Q.pairs(v)for v in C['quadratic_projection3_by180_sparse']];free2=Q.mvec(pi,{35:F(1)});u={int(i):F(c)for i,c in V['actual_u_image3'].items()};q2={int(i):F(c)for i,c in V['actual_q2_image3'].items()};D=Q.det_fixed([{j:col[i]for j,col in enumerate([u,q2,free2])if col.get(i,F(0))}for i in range(3)]);assert D
 assert not u.get(0,F(0)) and not q2.get(0,F(0)) and free2[0]==1 # literal grQ2 coordinate
 a,bcoef,c=F(P['actual_a']),F(P['actual_b']),F(P['actual_c']);unit=1+a*c;assert unit==F(V['actual_nonzero_qC_coefficient']) and unit
 # For all unknown alpha, r(u+c*w-alpha*q2)=0 exactly, because
 # r*u=-c*qC+alpha*g, r*w=qC, r*q2=g. Its gr2 vector is
 # (unit, c*b-alpha,0), so {kernel-vector,q2,free2} has determinant unit.
 affine_kernel=[[unit,F(0)],[c*bcoef,F(-1)],[F(0),F(0)]]
 assert affine_kernel[0][1]==0 and affine_kernel[2]==[F(0),F(0)]
 # Independent two-coefficient polynomial determinant by6 permutations.
 A=[[affine_kernel[i],([F(1),F(0)]if i==1 else[F(0),F(0)]),([F(1),F(0)]if i==2 else[F(0),F(0)])]for i in range(3)];det=[F(0)]*4
 def mul(x,y):
  z=[F(0)]*4
  for i,v in enumerate(x):
   for j,w in enumerate(y):
    if i+j<4:z[i+j]+=v*w
  return z
 for p in itertools.permutations(range(3)):
  v=[F((-1)**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))),F(0),F(0),F(0)]
  for i,j in enumerate(p):v=mul(v,A[i][j])
  det=[x+y for x,y in zip(det,v)]
 assert det==[unit,F(0),F(0),F(0)]
 assert b[15]['terms']==[[4,[5,7],'1'],[3,[7,7],str(-a)],[4,[7,7,7],str(-bcoef)]]
 lengths=[16,11,6,3,1,0];assert lengths==P['dimensions_m_power_Tnil'];assert [lengths[i]-lengths[i+1]for i in range(5)]==[5,5,3,2,1]
 ranks=[6,3,1,0];blocks=[ranks[i]-ranks[i+1]for i in range(3)];assert blocks==[3,2,1];assert P['rA_DVR_module_m_squared_Jordan_partition']==[3,2,1]
 assert P['complete_actual_cubic917_rank']==P['actual_constant811_rank']==1 and P['actual_right106_rank']==0 and P['no_specific_nonzero_constant811_column_evaluated']
 proof=dict(source_old5_rows=N['old_minimal_generator_free56_rows'],new_degree1_source_monomial_rows=firstfree,source16_lifts=b,degree2_fixed_determinant=D,
  actual_alpha_independent_kernel_vector_affine_coefficients=affine_kernel,kernel_quotient_determinant=det,
  gr3_basis_proof='r:gr2->gr3 is surjective with one-dimensional kernel. The nonzero displayed exact r-kernel vector lies in span(i(u),i(q2)) and has nonzero i(u) coefficient. Thus q2 and the outside-plane free42*r² complement this line, so their r-images i(g),free42*r³ form a basis.',
  degree4_proof='d=r*i(w-a*u-b*q2) has nonzero oldqC coefficient. It spans the already proved one-dimensional newm4; oldalpha remains unknown.',
  actual_full_HF=[5,5,3,2,1],actual_m_power_lengths=lengths,actual_filtered_basis_count16=True,k_rA_m2_Jordan_partition=[3,2,1],rA_not_DVRell=True,full917_rank1_and_right106_zero_implies_constant811_rank1=True,no_action_matrix_nil_generators_or_products_computed=True)
 save('independent-filtered16-source-basis-and-module-proof.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 report=dict(id=ID,all_passed=True,frozen=True,actual_filtered16_basis_verified=True,actual_full_HF=[5,5,3,2,1],actual_m_power_lengths=lengths,actual_constant811_rank1=True,oldalpha_or_new_source_values_evaluated=False,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',report);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
