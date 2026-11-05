#!/usr/bin/env python3
"""Unpaired leading K-square SOURCE tensor; no U²/US/S² value or rank."""
import hashlib,itertools,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63'
PATHS={
 'point':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63/five-cache-functionals-and-genuine-fourteen-kernel-source.json',
 'point_peer':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_20261105_v63_REVIEW/report.json',
 'Pi':'reports/AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63/actual512-origin-values-and-Pi0.json',
 'Pi_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/literal-known-cubic-input-and-exact24-row-reduction.json',
 'quad_peer':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW/report.json',
 'zero3':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json',
 'zero3_peer':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW/report.json',
 'seven':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63/actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json',
 'q2_peer':'reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63_REVIEW/report.json',
 'H_order':'reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63/specific-H-template-order-and-raw-domain-companion.json',
 'product':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/actual-nil-product44-and-cube-DAG.json',
 'product_peer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
 'upper':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/actual-trace-first-gate-and-entire-conductor-proof.json',
 'Loewy':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/exact-Loewy16-and-explicit-filtered-basis.json',
 'Loewy_peer':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json',
 'F2':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63/actual-head4-z5-map-and-F2-kernel-coefficient.json',
 'F2_peer':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_20261105_v63_REVIEW/report.json',
 'operator':'reports/AMBIENT_B18_C_KSQUARE_THREE_PAIR_OPERATOR_SOURCE_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ht(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
PRE=0;STORED=0
def add(a,b):
 global PRE,STORED
 h=2*max(ht(a),ht(b))+1;PRE=max(PRE,h);assert h<=1153
 v=a+b;STORED=max(STORED,ht(v));assert ht(v)<=384;return v
def mul(a,b):
 global PRE,STORED
 h=ht(a)+ht(b);PRE=max(PRE,h);assert h<=1153
 v=a*b;STORED=max(STORED,ht(v));assert ht(v)<=384;return v
def put(d,k,v):
 if v:
  w=add(d.get(k,Q(0)),v)
  if w:d[k]=w
  elif k in d:del d[k]
def mon(xs):
 e=[0]*8
 for x in xs:e[x]+=1
 return tuple(e)
def shift(e,i):v=list(e);v[i]+=1;return tuple(v)
M2=[mon(x)for x in itertools.combinations_with_replacement(range(8),2)];I2={m:i for i,m in enumerate(M2)}
PAIRS=list(itertools.combinations(range(8),2));IP={p:i for i,p in enumerate(PAIRS)}
def H3(v):
 """Degree2 of actual whole-vector H on a homogeneous degree3 input."""
 out={}
 for sec,n in enumerate([1,1,2]):
  rows=[{m:c for (a,m),c in v.items()if a==8*sec+i}for i in range(8)]
  hs=[{}for _ in range(7)]
  for e,c in rows[7].items():
   active=[i for i in range(7)if e[i]]
   if not active:continue # actual cycle makes the excluded pure-r part0
   i=active[0];f=list(e);f[i]-=1;put(hs[i],tuple(f),c)
  w=[row.copy()for row in rows[:7]]
  if n==1:
   for i in range(7):
    for e,c in hs[i].items():put(w[i],shift(e,7),c)
  for i in range(7):
   for e,c in hs[i].items():put(out,(28*sec+IP[i,7],e),c)
  for i,j in itertools.combinations(range(7),2):
   for row,var,sign in[(w[j],i,1),(w[i],j,-1)]:
    for e,c in row.items():
     if e[var]:
      f=list(e);f[var]-=1;f=tuple(f)
      put(out,(28*sec+IP[i,j],f),mul(c,Q(sign*e[var],sum(f[:7])+2)))
 return out
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'point':['source_m0_24_by44','head42_origin_cycle44_sparse','all24_m0_e43_entries_zero'],
  'Pi':['constant_minimal_projection5_by128_sparse'],
  'quad':['quadratic_projection3_by180_sparse','quadratic_NF_basis180_indices'],
  'zero3':['actual_three_minimal_first_columns40','actual_old_free56_labels'],
  'seven':['actual_correction_coefficients'],
  'product':['product_source_formula','mu_trace_order','actual_beta8_C0_node_ids','actual_Qrad_even36_R0_component_root_ids'],
  'upper':['all_Qrad_even_coefficients_in_J','actual_joint_ideal','J_annihilates_K_squared'],
  'Loewy':['explicit_sixteen_filtered_basis_lifts'],
  'F2':['known_actual_u_equals_e3r_image3','known_actual_q2_equals_e4r2_image3','fixed_pair_determinant'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 assert S.count(d)+40000<150000
 m=[[Q(c)for c in row]for row in d['point']['source_m0_24_by44']]
 v42={i:Q(c)for i,c in d['point']['head42_origin_cycle44_sparse']}
 assert v42.get(43,Q(0))==0 and all(row[43]==0 for row in m)
 for row in m:
  total=Q(0)
  for j,v in v42.items():total=add(total,mul(row[j],v))
  assert total==0
 P=[{i:Q(c)for i,c in row}for row in d['Pi']['constant_minimal_projection5_by128_sparse']]
 assert P[0]=={42:Q(1)}and P[1]=={43:Q(1)}
 A=[{i:Q(c)for i,c in row}for row in d['quad']['quadratic_projection3_by180_sparse']]
 assert d['quad']['quadratic_NF_basis180_indices']==[35,143,179]and A[0]=={35:Q(1)}
 # Source-derived universal fraction majorant: every tensor entry sums at
 # most28224 paths, each with one m0,Pi0,pi3 coefficient, at most TWO
 # v42 coefficients (the B term contains beta(v42_b)*T(v42_a)), factor2.
 arrs=[[x for row in m for x in row],[x for row in P for x in row.values()],[x for row in A for x in row.values()],list(v42.values())]
 env=[]
 for ar in arrs:
  den=math.lcm(*(x.denominator for x in ar));num=max(abs(x.numerator)*(den//x.denominator)for x in ar);env.append([num,den])
 DEN=12*math.prod(v[1]for v in env)*env[-1][1];NUM=2*28224*math.prod(v[0]for v in env)*env[-1][0]
 H=max(NUM.bit_length(),DEN.bit_length());assert H<384 and 2*H+1<1153
 Wq={i:mul(Q(2 if i>=36 else 1),x)for i,x in v42.items()}
 Wb={}
 for j in range(9):
  for a,b in[(4*j,4*j+1),(4*j+2,4*j+3)]:
   if a in v42:Wb[b]=v42[a]
 beta={mon([i-36]):x for i,x in v42.items()if 36<=i<43}
 assert all(sum(e)==1 for e in beta)
 def extract(raw):
  vv={}
  for (j,e),c in raw.items():
   for i in range(24):
    if m[i][j]:put(vv,(i,e),mul(m[i][j],c))
  hh=H3(vv);hh2={}
  for (j,e),c in hh.items():
   for t in range(5):
    if 44+j in P[t]:put(hh2,(t,e),mul(P[t][44+j],c))
  out=[Q(0)]*3
  for (t,e),c in hh2.items():
   ix=36*t+I2[e]
   for a in range(3):
    if ix in A[a]:out[a]=add(out[a],mul(A[a][ix],c))
  assert out[0]==0
  return out
 tensor=[];maxscratch=0
 for mi,e in enumerate(M2):
  for qvar in range(9):
   raw={}
   if qvar<8:
    for j,c in Wq.items():put(raw,(j,shift(e,qvar)),c)
   else:
    for f,b in beta.items():
     ex=tuple(x+y for x,y in zip(e,f))
     for j,c in Wb.items():put(raw,(j,ex),mul(b,c))
   values=extract(raw);assert all(ht(x)<=H for x in values)
   # These are unpaired linear extraction coefficients on S_ab and formal
   # q_k. L_U,L_S are deliberately NOT substituted; no pair value exists.
   tensor.append({'quadratic_monomial_index':mi,'formal_q_index_or_B_constant':qvar,'rows3':[str(x)for x in values]})
   maxscratch=max(maxscratch,S.count(raw)+24*120+84*36+5*36)
 data={
  'scope':'Unpaired QQ-linear SOURCE tensor for leading degree2 of the three genuine products. NO L_U/L_S substitution, NO Qeven1 coefficient readback, NO pair value/rank.',
  'source_namespace_bindings':{k:{'path':p,'sha256':pins[p]}for k,p in PATHS.items()},
  'actual_S_frame':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
  'quadratic_monomials36':[list(x)for x in M2],
  'literal_v42':[[i,str(x)]for i,x in v42.items()],
  'Wq44_sparse':[[i,str(x)]for i,x in Wq.items()],
  'Wbeta44_sparse':[[i,str(x)]for i,x in Wb.items()],
  'linear_B_h0':[[list(e),str(x)]for e,x in beta.items()],
  'all324_unpaired_extraction_coefficients':tensor,
  'QQ_scalar_symbols_not_Rbar_coefficients':'For qvar0..7 use formal q_k as scalar coefficient of x_k in Qeven1*v42_a. For qvar8 use the B term. Each quadratic monomial is multiplied by x_k or the full linear B polynomial BEFORE whole H, then Pi0 and pi3. NO S-linearity of H.',
  'future_pair_reconstruction':'For pairU/V: form homogeneous quadratic L_U*L_V over S, then contract its36 coefficients with this tensor and scalar q0..7. qvar8 is the already-known B term. No d_U/d_S contribution since whole m0*e43=0.',
  'leading_input_formula':'c_X1=L_X*v42+d_X*r*e43 forX=U,S; Gamma0=0 killssection1. L_U=oldZERO3h3_head0-a3*oldZERO3h1_head0; L_S=oldZERO3h2_head0-oldZERO3h3_head0/4-aS*oldZERO3h1_head0.',
  'all_higher_leakages_excluded':'Productc>=3; raw128 degree2=(0,H(m0*c3)_2). m1*c3/inputc4 giveorder>=3 afterH, and Pi1 onhonestdegree2 starts3. pi3 is thecompleteactualgr2projection. NOGamma2/m1/Pi1/unknownlowercarry enters leadingclass.',
  'pure_plane_first_NF_coordinate_zero':'Pi0firstrow=[e42,0] and pi3row0=delta(head0*r²), so all324 tensor first coordinates arezero BEFOREanyrequestedpairing.',
  'sufficient_rank2_criterion':'If actualthreeleadingcolumns have rank2, theyspan the plane<i(v),i(q2)>. r onthisplane hasrank1 andr*i(q2)=i(g) nonzeroingr3. Choosepairlinearcombinationswithindependentdegree2initials; one r-multiple hasnonzero degree3initial. TheseTHREEactualclassesareindependent byfiltration, no strictembedding assumption. SinceallareinH0bar dimension3, H0=D*K². ThisdoesNOTclaimthisantecedent.',
  'available_actual_first_cache':d['zero3'],
  'actual_correction_coefficients':d['seven']['actual_correction_coefficients'],
  'q_callback_FRONTIER':'The eight q_k remainFORMAL. NextSOURCE must eitherprovetherequiredsensitivitiescancel afterL substitution or bindtheactualQevenfirstcombo toexistingnativeproperC36/radicalpolynomial cache BEFOREsignedpairvalues.',
  'fraction_majorant_before_arithmetic':{'input_envelopes_U_D':env,'common_denominator':DEN,'numerator_upper':NUM,'H':H,'preoperation':2*H+1,'path_count_upper':28224},
  'specific_homotopy':'Exactactualn=(1,1,2)Hdegree3→2,orderedoriginalzHadamard,fullw=vi+r^nHi,positiveEuler_z+2 inverse. Bothr^n andderivativesperformedinsideH.',
  'all_product_values_ranks_H0_nil_UNCOMPUTED':True,'no_signed_value_GO':True,'selective_receipts':receipts,
 }
 owned=S.count(d)+S.count(tensor)+S.count(data)+maxscratch+20000;assert owned<150000
 data['whole_SOURCE_owned_upper']=owned;data['max_stored_H']=STORED;data['max_before_operation_H']=PRE
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'unpaired-leading-three-pair-extraction-tensor-and-rank2-sufficient-criterion.json';p.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'UNPAIRED_SOURCE_LEADING_EXTRACTION_AND_CONDITIONAL_RANK2_CRITERION','input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'tensor_entries':len(tensor),'max_stored_H':STORED,'max_before_operation_H':PRE,'source_majorant_H':H,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
