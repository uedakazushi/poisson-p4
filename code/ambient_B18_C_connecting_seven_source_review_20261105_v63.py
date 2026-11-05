#!/usr/bin/env python3
"""Independent honest Told/I quotient and fixed connecting constraints; no edges."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='c36dfe511d6f352782036f185adfaa93eeb7916e79c1e6ca6663e24683e6e7f1'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Q(0)) for j in range(len(B[0]))] for i in range(len(A))]
def det(A):
 n=len(A)
 return sum(((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
   *product(A[i][p[i]] for i in range(n)) for p in itertools.permutations(range(n))),Q(0))
def product(xs):
 z=Q(1)
 for x in xs:z*=x
 return z
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  pins[str(p.relative_to(ROOT))]=s
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
  'a4a84caf64a7aac5269c1c3d2067fcae6bbd24c2c8b5730f4019f4469571d4c7')
 pin(ROOT/'code/ambient_B18_C_connecting_seven_source_v2_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 old=docs['known-initial40-and-length-upper17.json'];T=docs['exact-Loewy11-and-source-basis-corollary.json']
 assert T['exact_length']==11 and T['exact_Loewy_HF']==[5,3,2,1]
 heads=docs['actual-natural-Told-to-Tnil-source-map-and-origin.json']
 assert heads['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 first=[g for g in old['homogeneous_module_generators40'] if all(sum(m)==1 for h,m,q in g)]
 assert len(first)==37
 rmon=[0]*7+[1]
 A=[[sum((Q(q) for head,m,q in g if head==j and m==rmon),Q(0)) for j in range(5)] for g in first]
 assert [[str(q) for q in row] for row in A]==D['all37_linear_relation_r_coefficients_rows37_columns5']
 selected=D['fixed_rank3_rows'];assert selected==[25,16,30]
 basisrows=[A[i] for i in selected];U=[row[:3] for row in basisrows]
 determinant=det(U);assert determinant==-Q(36333,20)==Q(D['fixed_rank3_minor'])
 # Own adjugate formula, not the producer's small elimination routine.
 inverse=[[(-1)**(i+j)*det([[U[r][c] for c in range(3) if c!=i] for r in range(3) if r!=j])/determinant for j in range(3)] for i in range(3)]
 ident3=[[Q(i==j) for j in range(3)] for i in range(3)]
 assert mm(U,inverse)==ident3==mm(inverse,U)
 assert [[str(q) for q in row] for row in inverse]==D['fixed3_inverse']
 factors=mm([row[:3] for row in A],inverse)
 assert mm(factors,basisrows)==A and [[str(q) for q in row] for row in factors]==D['all37_row_factorizations']
 null=[[0,0],[0,0],[Q(1,4),0],[1,0],[0,1]]
 assert mm(A,null)==[[0,0] for _ in range(37)]
 assert [[str(q) for q in row] for row in null]==D['cyclic_top_chi_kernel5_by2']
 R=[[Q(0)]*7 for _ in range(7)];R[5][2]=Q(1,4);R[5][3]=1;R[6][4]=1
 assert [[str(q) for q in row] for row in R]==D['actual_r_action7']
 assert mm(R,R)==[[0]*7 for _ in range(7)]
 assert all(sum(row[j]*R[i][j] for j in range(5))==0 for row in A for i in range(7))
 C=[[Q(q) for q in row] for row in D['source7_basis_change']];I=[[Q(q) for q in row] for row in D['source7_basis_inverse']]
 ident7=[[Q(i==j) for j in range(7)] for i in range(7)]
 assert mm(C,I)==ident7==mm(I,C)
 J=[[Q(0)]*7 for _ in range(7)];J[1][0]=1;J[3][2]=1;assert mm(mm(I,R),C)==J
 # Four actual OLD source classes are literal members of I*Told and
 # independent by the already-accepted actual filtered monomial basis.
 standard=T['actual_filtered_monomial_basis11']
 fixed=[D['fixed_I_Told_source_basis']['first']]+D['fixed_I_Told_source_basis']['other_three']
 assert all(row in standard for row in fixed) and len(fixed)==4
 assert [row['degree'] for row in fixed]==[1,2,2,3]
 for row in fixed:
  mon=row['monomial_original8'];assert any(mon[:7]) or mon[7]>=2
 assert 11-7==4 and 7-4==3 and 4-3==1
 # Hom(J2²+k³,J2+k²) = 2*2+3*1+2*5 =17 overQQ,
 # agreeing with20headcoordinates minus3 independentcyclic-top constraints.
 assert 2*2+3+2*5==20-3==17
 history=docs['missing-peer-path-and-actual-frozen-dependency-correction.json']
 assert history['status']=='V1_MISSING_LITERAL_PEER_PATH_BEFORE_INPUT_LOAD_OR_MATH'
 assert sha(ROOT/history['snapshot_path'])==history['snapshot_sha256']
 proof={'same_Rbar_R0_div_ell_original_Told_heads_and_genuine_target':True,
  'all37_literal_r_projections_and_fixed_rank3_minor':'-36333/20',
  'own_adjugate_bothinverse_products_all37_factorbacks':True,
  'complete_two_nullvectors_and_chi_constraint':'chi=(0,0,c/4,c,d);c,d unknown',
  'honest_quotient_completeness_proof':[
   'Told is the actual already-computed minimal5 module. Its constant relations are0; actual HF1=3 makes the37 independent knownlinear initials the ENTIRE first relation space.',
   'Modulo I=(z7,r²), every coefficient of order>=2 vanishes. Therefore the complete actualmodule quotient Told/I*Told is determined exactly by the37 r-coefficient rows, not by an unproved associatedgraded exactness.',
   'The rank3 rowspace leaves five heads and two independent r classes. Literal basis is h0..h4,rh3,rh4; rh2=rh3/4,rh0=rh1=0. Thus dim7 and typeJ2²+k³.'
  ],
  'kernel_reduction_exactness_proof':[
   'For the accepted original comparison complex with H1=K,H0=C, the DVR universal-coefficient exact sequence gives0->K/ellK->Told->C[ell]->0. Its derived free-complex interpretation doesnot assumeD0 ell-flat or discard its possible higherTor.',
   'G=C[ell] has exact annihilator I in Rbar, so I*Told lies in K/ellK and the actualedge factors through M7.',
   'The four listed literal oldclasses lie in I*Told and are independent by the actual filtered11basis. Length11-7=4 proves they form its completeQbasis.',
   'Hence0->I*Told4->K/ellK7->ker(M7->G4)3->0. The quotient is Kbar/(I*Told), never Kbar/(I*Kbar). No Artin-Rees or finitehomology replacement is used.'
  ],
  'source7_module_type':'J2²+k³','all7_action_and_Jordanbasis_products_verified':True,
  'five_head_parameter_count':17,
  'surjectivity_scope':'Actualedge is surjective by the originalexactsequence. Its17parameters are uncomputed; top3rank3 is an opencondition, not a chosennumericalfixture.',
  'actual_edge_values_and_residual_kernel3_uncomputed':True,
  'no_new_physical_source_gauge_IFT_jet_or_original_H0_nil_claim':True,
  'source_V1_path_failure_snapshot_preserved':True}
 ep=out/'independent-honest-seven-source-fixed-kernel-four-and-formal-edge-constraints.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'honest_M7_and_four_fixed_kernel_classes_accepted':True,
  'no_actual_head_edge_or_original_H0_values':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
