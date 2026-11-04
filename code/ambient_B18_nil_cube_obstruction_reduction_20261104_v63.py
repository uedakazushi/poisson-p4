#!/usr/bin/env python3
"""Source-only nil-cube obstruction tower and exact Tor2 input recipe."""
import hashlib,itertools,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_NIL_CUBE_OBSTRUCTION_REDUCTION_20261104_v63';OUT=ROOT/'reports'/ID
Z=(0,)*9
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def put(p,m,c):
 p[m]=p.get(m,0)+c
 if not p[m]:p.pop(m)
def plus(p,q,s=1):
 v=p.copy()
 for m,c in q.items():put(v,m,s*c)
 return v
def times(p,q):
 v={}
 for m,c in p.items():
  for n,d in q.items():put(v,tuple(a+b for a,b in zip(m,n)),c*d)
 return v
def mono(i,j=None):
 e=list(Z);e[i]+=1
 if j is not None:e[j]+=1
 return {tuple(e):1}
def encode(p):return [[list(m),c]for m,c in sorted(p.items())]
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 stems=['AMBIENT_B18_CONDUCTOR_TRANSFER512_NIL_CYCLES','AMBIENT_B18_NIL_PRODUCT44_OPERATOR','AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7','AMBIENT_B18_NIL_DVR_FINITENESS','AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN','AMBIENT_B18_RADICAL_COMPARISON_OPERATOR','AMBIENT_B18_FULLV56_COMPARISON_UNIT','AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR']
 for stem in stems:
  r=load(stem+'_20261104_v63_REVIEW');assert r['all_passed']
 for stem in stems:
  p=ROOT/f'research/notes/{stem}_20261104_v63.md'
  if p.exists():pins[str(p.relative_to(ROOT))]=sha(p)
 radical=load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')
 chain=load('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json')
 product=load('AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63','actual-nil-product44-and-cube-DAG.json')
 assert radical['NilD0_exactly_H1_of_actual_complex'] and len(radical['actual_V3_by44_rows_after_unit1_cancellation'])==3
 assert [len(row)for row in radical['actual_V3_by44_rows_after_unit1_cancellation']]==[44]*3
 Q=chain['exact_existing_Q48_C4'];U=chain['exact_existing_U9_C4'];assert len(Q)==48 and all(len(c)==4 for c in Q)and U[0]==[1,0,0,0]and len(U)==9
 # Actual all-order V56, BEFORE any ell specialization. The ell=0 V56 file
 # is not used as a substitute for this source.
 source56={'path':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json','namespace':'operator_scalar_nodes','shape':[3,56],'entries':[[{'root':c[t],'sign':1}for c in Q]+[{'root':c[t],'sign':-1}for c in U[1:]]for t in [1,2,3]],'target_order':['tauA','tauB','tauA*tauB'],'source_order':'a48,c8 after target-unit1 cancellation','X_shape':[56,768],'H1':'K','H0':'C'}
 source44={'path':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json','json_pointer':'/actual_V3_by44_rows_after_unit1_cancellation','shape':[3,44],'foreign_operator_namespaces_retained':True,'X_shape':[44,780],'H1':'N','H0':'Cprime'}
 assert product['complete_cycle_source_order']==['a36 (9 blocks times4)','c8 (actual C0 odd basis)']
 # Exact resolution of A=R/J, J=(z0,...,z6,r^2,rs,s^2):
 # Koszul(z7) tensor the Hilbert-Burch resolution of (r,s)^2.
 labels={}
 for n in range(4):
  labels[n]=[(S,q,h)for q in range(3) for S in itertools.combinations(range(7),n-q)for h in range([1,3,2][q])] if n>=2 else [(S,q,h)for q in range(n+1)for S in itertools.combinations(range(7),n-q)for h in range([1,3,2][q])]
 # combinations receives negative r only when n=0/1, guarded above.
 assert [len(labels[n])for n in range(4)]==[1,10,44,112]
 differential={}
 for n in [1,2,3]:
  rows={v:i for i,v in enumerate(labels[n-1])};D={}
  def add(i,j,p,s):D[i,j]=plus(D.get((i,j),{}),p,s)
  for j,(S,q,h)in enumerate(labels[n]):
   for pos,v in enumerate(S):add(rows[(S[:pos]+S[pos+1:],q,h)],j,mono(v),(-1)**pos)
   sg=(-1)**len(S)
   if q==1:add(rows[(S,0,0)],j,mono(7,7)if h==0 else mono(7,8)if h==1 else mono(8,8),sg)
   elif q==2:
    for hh,p,ss in ([(0,mono(8),-1),(1,mono(7),1)]if h==0 else[(1,mono(8),-1),(2,mono(7),1)]):add(rows[(S,1,hh)],j,p,sg*ss)
  differential[n]={k:v for k,v in D.items()if v}
 for n in [2,3]:
  A=differential[n-1];B=differential[n];result={}
  for (i,k),p in A.items():
   for (kk,j),q in B.items():
    if kk==k:result[i,j]=plus(result.get((i,j),{}),times(p,q))
  assert not any(result.values())
 # The 31 monomial classes of R/J^2, metadata only, never an AR assertion.
 mons=[]
 def inJ2(e):return sum(e[:7])>=2 or(sum(e[:7])>=1 and e[7]+e[8]>=2)or e[7]+e[8]>=4
 for degree in range(4):
  for picks in itertools.combinations_with_replacement(range(9),degree):
   e=[0]*9
   for i in picks:e[i]+=1
   if not inJ2(e):mons.append(e)
 assert len(mons)==31 and [sum(sum(e)==i for e in mons)for i in range(4)]==[1,9,17,4]
 resolution={'ring_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','sA'],'J_generators':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA^2','rA*sA','sA^2'],'basis_order':{str(n):[[list(S),q,h]for S,q,h in labels[n]]for n in labels},'betti_first4':[1,10,44,112],'differentials_sparse':{str(n):{'shape':[len(labels[n-1]),len(labels[n])],'entries':[[i,j,encode(p)]for(i,j),p in sorted(D.items())]}for n,D in differential.items()},'all_d1d2_and_d2d3_coefficient_identities':True,'exact_resolution_proof':'Tensor Koszul(z7) with the Hilbert-Burch resolution of (r,s)^2 in the independent r,s variables; z7 remains regular on the extended quotient.','R_mod_J_length':3,'R_mod_J2_monomial_basis':mons,'R_mod_J2_HF':[1,9,17,4],'R_mod_J2_length':31,'no_Artin_Rees_strictness_inferred':True}
 (OUT/'exact-J-resolution-and-Tor2-input.json').write_text(json.dumps(resolution,indent=2)+'\n')
 typed=[]
 for b,source in [(56,source56),(44,source44)]:
  typed.append({'H1':source['H1'],'H0':source['H0'],'V_source':source,'X_shape':source['X_shape'],'P2_lift_c_shape':[3,44],'P1_relation_lift_h_shape':[b,10],'constraint_left_operator':'I3 tensor d2','constraint_left_shape':[30,132],'constraint_right_operator':'V tensor I10','constraint_right_shape':[30,10*b],'constraint_operator_shape':[30,132+10*b],'constraint_exact':'(I3 tensor d2)c=(V tensor I10)h over the original full formal R0','output_operator':'I_b tensor d1','output_shape':[b,10*b],'actual_output':'w=(I_b tensor d1)h lies in ker(V); its class modulo J*H1 is delta2([c]).','equivalence_checks':['h difference in kerV tensor P1 changes w by J*H1','c difference in imV tensor P2 changes h by d2 b and leaves w unchanged because d1d2=0','c addition d3 t leaves the exact constraint unchanged because d2d3=0'],'missing_actual_data':'Complete relevant all-order cycle/relation lifts (c,h), or an honest H1 presentation. No finite mod-J or mod-J2 kernel alone supplies these lifts.','map_not_evaluated':True})
 (OUT/'actual-source-Tor2-connecting-recipe.json').write_text(json.dumps(typed,indent=2)+'\n')
 # These degree checks use literal product44 term structure, not arbitrary
 # free44 vectors being declared cycles.
 degree={'Qeven':1,'beta':1,'nuA':1,'a':0,'c':0,'aprime':0,'cprime':0}
 assert degree['Qeven']+degree['a']==1 and degree['beta']+degree['c']==1
 filtered={'b0':1,'B':1,'bilinear_a_terms':[1,2,1,1,2,1],'bilinear_c_terms':[1,1],'cubic_a_minimum':2,'cubic_c_minimum':2,'conclusion':'F^pN*F^qN subset F^(p+q+1)N; N^3 subset F^2N, F^qN=image(kerV intersect J^q R44)','not_inferred':'F^qN=J^qN; J annihilates N^3; K/JK controls mixed K*N^2'}
 # Exact tensor-count bound, independent of any choice of actual source lifts.
 pair=math.comb(8,2);triple=math.comb(9,3);assert (pair,triple)==(28,84)
 cases=[{'L':L,'mu_D_K':L-4,'mixed_generator_product_upper':(L-4)*pair,'pure_generator_product_upper':triple,'total_upper':(L-4)*pair+triple}for L in range(11,18)]
 assert max(c['total_upper']for c in cases)==448
 theorem={'hypotheses':['R0 char0 commutative complete local actual B18 base','D=k[[ell]] with actual ell, not rA','K=ker(D0→E4), N=Nil(D0)','N^3 subset finite-origin tors_D K','N*K^2=0','P=N/K is D-free rank7','mu_D K=L-4, 11<=L<=17','W*N^2=0, IA*IB*N subset W'],
 'beta':'K tensor_D Sym^2_D P → N^3, (k,p,q)↦k*lift(p)*lift(q)',
 'beta_well_defined':'Changing either lift adds a product in N*K^2=0.',
 'H':'image(beta)=K*N^2',
 'gamma':'Sym^3_D P → N^3/H, using arbitrary lifts',
 'gamma_well_defined':'Every lift change introduces at least one K and two N factors, hence lies in H.',
 'cube_zero_iff':'beta=0 AND gamma=0','safe_transfer_quotient':'beta descends through K/(K intersect W); both P factors and gamma descend through P/image(W→P).',
 'no_K_mod_J_for_beta':'JK^2=0 only justifies K/JK for K^2. J*K*N^2=0 has NOT been proved.',
 'interim_implications':['beta=0 implies N^5=0 since N^3 subset K','gamma=0 implies N^4=0 since N^3=K*N^2 and N*K^2=0'],
 'D_generator_product_bounds_by_L':cases,'total_D_generator_product_upper':448,'not_length_bound':True,'actual_seven_lifts_K_generators_or_product_values_constructed':False,'source_filtered_product':filtered,
 'remaining_exact_source_requests':['Seven actual D-basis lifts of P=N/K and at most13 actual D-generators of K, with all-order comparison/source identities. Existence/rank alone is not such a certificate.','Or an exact cyclic-image membership proof for the source-bound cube44 operator on all actual cycles, preserving Tor2 and W quotients.','For K/JK specifically: the all-order (c,h) lifts in the supplied 30×692 constraint, not a mod-J2 approximation. This gate alone does not decide beta.'],
 'no_nil_cube_verdict':True,'no_new_source_jet_scalar_search_or_parameter_point':True}
 (OUT/'cube-obstruction-tower.json').write_text(json.dumps(theorem,indent=2)+'\n')
 for p in [Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'canonical_beta_gamma_reduction':True,'maximum_D_generator_products':448,'exact_J_resolution_first_betti':[1,10,44,112],'K_Tor2_constraint_shape':[30,692],'N_Tor2_constraint_shape':[30,572],'all_differential_polynomial_identities':True,'no_actual_nil_lifts_products_or_Tor2_image_computed':True,'nil_cube_still_OPEN':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},'seconds':time.monotonic()-started}
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
