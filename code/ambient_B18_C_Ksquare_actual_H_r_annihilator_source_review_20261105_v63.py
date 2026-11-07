#!/usr/bin/env python3
"""Independent SOURCE audit: literal max-prefix H, raw Pi1, true clearing.

No LU/LS/q coefficient is substituted and no product value is requested.
Only the standard-library arithmetic and our previous selective reader are used.
"""
import hashlib,itertools,json,math,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C_KSQUARE_ACTUAL_H_R_ANNIHILATOR_SOURCE_20261105_v63_REVIEW'
PINS={
 'reports/AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63/report.json':'e494a6c3d59d86626e3ca8c0a668748313ebda9b9487b4b47ab7c359ed325d95',
 'reports/AMBIENT_B18_C_KSQUARE_LEADING_LINEAR_SPAN_SOURCE_20261105_v63/report.json':'f3e62e5f9f0dadd30ae001268259270c50df6b431714d2a5d7d5a48e9f1ce8bf',
 'reports/AMBIENT_B18_C_KSQUARE_R_DUAL_C1_DOMAIN_DIAGNOSTIC_20261105_v63/report.json':'d89c88433ab109553c7f26608f9906cb16d31a302c416dd36400162ee064b434',
 'reports/AMBIENT_B18_C_KSQUARE_R_ANNIHILATOR_CORRECTED_SOURCE_20261105_v63/report.json':'b2c7149ae43bd2e006e1af411af97cba39ccb8ab038031f95d382c99aa1ea595',
 'reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63_REVIEW/report.json':'caed49bc0049782a568e30b738068869d864395d8ef754be6edfdb4d69606f6d',
}
FILES={
 'point':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63/five-cache-functionals-and-genuine-fourteen-kernel-source.json',
 'Pi':'reports/AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63/actual512-origin-values-and-Pi0.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/literal-known-cubic-input-and-exact24-row-reduction.json',
 'linear':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'zero3':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json',
 'q35':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/phase1-all35-nine-convolution-projections.json',
 'rightsource':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63/right106-full-first-span-and-cache-footprint.json',
 'mu':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63/complete-twelve-generator-kernel-and-two-cache-callback-gate.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'X':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'Hcode':'code/ambient_B18_honest_nil_T_operator_20261104_v63.py',
 'quartic':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/actual-eight-relations-and-complete-quartic-dual.json',
 'Loewy':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/report.json',
 'compact':'reports/AMBIENT_B18_C_KSQUARE_R_PAIR_DUAL_SOURCE_20261105_v63/compact-two-term-linear-clearing-and-sector0-row7-zero-factorization.json',
 'leading':'reports/AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63/unpaired-leading-three-pair-extraction-tensor-and-rank2-sufficient-criterion.json',
 'linear_result':'reports/AMBIENT_B18_C_KSQUARE_LEADING_LINEAR_SPAN_SOURCE_20261105_v63/unpaired-known174-membership-and-pure-r-H2-source.json',
}
PRE=STORED=0
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def count(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.keys()))+sum(map(count,x.values()))
 return 0
def height(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def op(a,b,multiply=False):
 global PRE,STORED
 if not a:return Q(0) if multiply else b
 if not b:return Q(0) if multiply else a
 na,da=abs(a.numerator).bit_length(),a.denominator.bit_length()
 nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
 h=max(na+nb,da+db) if multiply else max(max(na+db,nb+da)+1,da+db)
 PRE=max(PRE,h);assert h<=1153
 v=a*b if multiply else a+b
 STORED=max(STORED,height(v));assert height(v)<=384
 return v
def add(a,b):return op(a,b)
def mul(a,b):return op(a,b,True)
def total(xs):
 z=Q(0)
 for x in xs:z=add(z,x)
 return z
def put(d,k,c):
 if c:
  v=add(d.get(k,Q(0)),c)
  if v:d[k]=v
  elif k in d:del d[k]
def env(vals):
 D=U=1
 for v in vals:
  D=math.lcm(D,v.denominator);U=max(U,(abs(v.numerator)+v.denominator-1)//v.denominator)
 return U,D
def prior(groups,paths):
 es=[env(g) for g in groups];D=math.prod(d for u,d in es);U=paths*math.prod(u for u,d in es)
 h=max((U*D).bit_length(),D.bit_length());assert h<=384 and 2*h+1<=1153
 return {'group_U_D':es,'absolute_path_bound':paths,'stored_H':h,'pre_H':2*h+1}
def get(k,f):return field(ROOT/FILES[k],f)
def sparse(row):return {i:Q(c) for i,c in row if Q(c)}
def mon(xs):
 e=[0]*8
 for x in xs:e[x]+=1
 return tuple(e)
M2=[mon(x) for x in itertools.combinations_with_replacement(range(8),2)]
I2={e:i for i,e in enumerate(M2)}
PAIR=list(itertools.combinations(range(8),2));PI={p:i for i,p in enumerate(PAIR)}
def actual_H2(v):
 """Literal max-index prefix, retaining n=2 correction before truncation."""
 out={}
 for s,n in enumerate((1,1,2)):
  h=[{} for _ in range(7)]
  last={e:c for (j,e),c in v.items() if j==8*s+7}
  for e,c in last.items():
   active=[i for i in range(7) if e[i]]
   # The whole actual input is a cycle. Its last component at z=0 is0.
   if not active:continue
   i=max(active);f=list(e);f[i]-=1;put(h[i],tuple(f),c)
  w=[{e:c for (j,e),c in v.items() if j==8*s+i} for i in range(7)]
  for i in range(7):
   for e,c in h[i].items():
    f=list(e);f[7]+=n;put(w[i],tuple(f),c)
    if sum(e)==2:put(out,(s*28+PI[i,7],e),c)
  for i,j in itertools.combinations(range(7),2):
   for p,deriv,sgn in ((w[j],i,1),(w[i],j,-1)):
    for e,c in p.items():
     if not e[deriv]:continue
     f=list(e);f[deriv]-=1;f=tuple(f)
     value=mul(c,Q(sgn*e[deriv],sum(f[:7])+2))
     if sum(f)==2:put(out,(s*28+PI[i,j],f),value)
 return out
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 pins={p:sha(ROOT/p) for p in list(PINS)+list(FILES.values())}
 for p,h in PINS.items():assert pins[p]==h
 sourcebytes=sum((ROOT/p).stat().st_size for p in pins);assert sourcebytes<32*1024**2
 results={};phase=[]
 # Stage1: raw174 string inverse alone; retain only four required rows.
 labels=get('linear','all280_source_labels_relation_and_variable')
 ids=[labels[8*i][0] for i in range(35)];assert ids==get('rightsource','selected35_actual_first_relations')
 cp=get('linear','fixed_rank_pivot_columns');rr=get('linear','fixed_rank_pivot_rows')
 iv=get('linear','fixed_minor_inverse');assert len(iv)==len(rr)==174
 parse174=count((iv,cp,rr,labels))+18000;assert parse174<150000
 wanted=(623,627,673,675)
 chosen={i:sparse(zip(rr,iv[cp.index(8*i+7)])) for i,j in enumerate(ids) if j in wanted}
 del iv
 q35=[Q(row[0]) for row in get('q35','all35_projected_M2_lower_representatives')]
 assert [ids[i] for i,c in enumerate(q35) if c]==list(wanted)
 assert [r['source_relation'] for r in get('q35','all35_source_input')]==ids
 assert get('rightsource','all_nonzero_convolution_duals')[0]==[0,7,[[35,'1']]]
 assert 'separately' in get('mu','callback_means')
 assert 'No retained-wedge' in get('mu','individual_callback_U_sensitivity_zero_reason')
 assert 'sector0' in get('mu','individual_callback_V_sensitivity_zero_reason')
 clear_prior=prior([[c for row in chosen.values() for c in row.values()],q35],174*35)
 weighted={}
 for i,row in chosen.items():
  for j,c in row.items():put(weighted,j,mul(q35[i],c))
 assert weighted=={7:Q(50,3),29:Q(125,3)}
 assert [[j,str(c)] for j,c in sorted(weighted.items())]==get('compact','weighted_M2_head0r2_times_fixed174_r_section')
 results['true_M2_weighted_fixed174_r_section']=[[j,str(c)] for j,c in sorted(weighted.items())]
 results['selected35_source_ids']=ids
 results['individual_direction0_cached_q35_actual']='Literaldual[[35,1]], separateU/V sensitivities0; relation-side -8/3 remains inside its accepted relation domain.'
 phase.append(['174inverse_string_parse_and_four_rows',parse174])
 del chosen,labels,cp,rr,q35
 # Stage2: universal actual-max homotopy tensor, never substitute L/q.
 m=[[Q(c) for c in row] for row in get('point','source_m0_24_by44')]
 v=sparse(get('point','head42_origin_cycle44_sparse'));assert not v.get(43)
 P=list(map(sparse,get('Pi','constant_minimal_projection5_by128_sparse')))
 A=list(map(sparse,get('quad','quadratic_projection3_by180_sparse')))
 B=list(map(sparse,get('linear','exact_quotient_duals_rows_sparse')))
 assert P[0]=={42:Q(1)} and P[1]=={43:Q(1)} and A[0]=={35:Q(1)}
 assert get('quad','quadratic_NF_basis180_indices')==[35,143,179]
 assert all(not row[43] for row in m)
 assert all(total(mul(row[j],c) for j,c in v.items())==0 for row in m)
 Wq={j:mul(c,Q(2 if j>=36 else 1)) for j,c in v.items()}
 Wb={j+1:v[j] for j in range(0,36,2) if j in v}
 beta={mon([j-36]):c for j,c in v.items() if 36<=j<43}
 assert all(total(mul(m[7][j],c) for j,c in W.items())==0 for W in (Wq,Wb))
 hcode=(ROOT/FILES['Hcode']).read_text();assert 'i=max(ids)' in hcode
 tensor_prior=prior([[x for r in m for x in r],[c for r in P for c in r.values()],
                     [c for r in A+B for c in r.values()],v.values(),v.values(),[Q(1,60)]],2**30)
 data=[];maxscratch=0;rawhash=hashlib.sha256();pureH27=True
 for n,e in enumerate(M2):
  for channel in range(9):
   raw={}
   if channel<8:
    ex=list(e);ex[channel]+=1
    for j,c in Wq.items():put(raw,(j,tuple(ex)),c)
   else:
    for f,b in beta.items():
     ex=tuple(x+y for x,y in zip(e,f))
     for j,c in Wb.items():put(raw,(j,ex),mul(b,c))
   vv={}
   for (j,f),c in raw.items():
    for i,row in enumerate(m):
     if row[j]:put(vv,(i,f),mul(c,row[j]))
   # Strong universal factorization: whole sector0 last component0.
   assert not any(j==7 for j,f in vv)
   hh=actual_H2(vv);assert not any(j==27 for j,f in hh)
   mini={}
   for (j,f),c in hh.items():
    for t,row in enumerate(P):
     if row.get(44+j):put(mini,36*t+I2[f],mul(row[44+j],c))
   assert not any(j<36 for j in mini)
   va=[total(mul(c,row[j]) for j,c in mini.items() if j in row) for row in A]
   vb=[total(mul(c,row[j]) for j,c in mini.items() if j in row) for row in B]
   assert va==[0]*3 and vb==[0]*6
   assert total(mul(c,weighted[j]) for j,c in mini.items() if j in weighted)==0
   rawhash.update((json.dumps([[j,str(c)] for j,c in sorted(mini.items())],separators=(',',':'))+'\n').encode())
   data.append([n,channel,[str(c) for c in va],[str(c) for c in vb]])
   maxscratch=max(maxscratch,count((raw,vv,hh,mini,va,vb))+16000)
 owned_tensor=count((m,v,P,A,B,Wq,Wb,beta,weighted,data))+maxscratch+30000
 assert owned_tensor<150000;phase.append(['max_active_unpaired324_tensor',owned_tensor])
 zero3=get('zero3','actual_three_minimal_first_columns40');assert all(not Q(c.get('7','0')) for c in zero3)
 results['actual_prefix']='Literal original-z prefix chooses MAX active index; n=(1,1,2); n2 w correction retained before degree2 truncation.'
 results['all324_actual_max_final3_and_known174_quotient6_zero']=True
 results['all324_actual_max_columns']=data
 results['actual_max_raw_minimal_degree2_stream_sha256']=rawhash.hexdigest()
 results['literal_H27_degree2_whole_input_zero']=True
 results['LU_LS_r_coefficients_zero_from_corrected_ZERO3']=True
 results['convention_scope']='Frozen leading/min-active intermediate is not asserted literal. This independent MAX computation supplies the actual source result; no source bytes changed.'
 del m,v,P,A,B,Wq,Wb,beta,zero3,data,raw,vv,hh,mini
 # Stage3: reconstruct rawPi1 from original780X and81 unitblock, not effective relation-side C1.
 assert 42 not in get('m0','selected_source_columns24')
 X=get('X','all780_source_columns44_sparse_dual9');xp=get('X','X0_pivot_columns')
 fac=get('X','all780_X0_factorization_coefficients')
 wr=get('W','W_pivot_rows');wp=get('W','W_pivot_residual_columns');res=get('W','residual_original_X_columns')
 inv=get('W','W_fixed_minor_inverse');wall=get('W','all_W84_columns_sparse')
 owned_raw=count((results,X,xp,fac,wr,wp,res,inv,wall))+30000
 assert owned_raw<150000;phase.append(['780original_X_and81_raw_unitblock',owned_raw])
 ar=[]
 for col in X:
  c=next((dict(a) for j,a in col if j==42),{})
  assert not Q(c.get(0,'0'));ar.append(Q(c.get(8,'0')))
 # Store strings once, convert entries only when read. No expanded full918 or inverse construction.
 raw_prior=prior([ar,(Q(c) for row in fac for c in row),(Q(c) for row in inv for c in row)],81*19)
 sel=[]
 for k in wp:
  j=res[k];sel.append(add(ar[j],-total(mul(Q(c),ar[p]) for c,p in zip(fac[j],xp))))
 actual={}
 for k,t in enumerate(wr):put(actual,t,-total(mul(c,Q(inv[j][k])) for j,c in enumerate(sel)))
 assert actual=={27:Q(8,3)}
 # Independently multiply every saved inverse row/column against the ORIGINAL selected W minor.
 Wsel=[sparse(wall[k]) for k in wp]
 for i in range(81):
  for j in range(81):
   left=total(mul(Wsel[k].get(wr[i],Q(0)),Q(inv[k][j])) for k in range(81))
   right=total(mul(Q(inv[i][k]),Wsel[j].get(wr[k],Q(0))) for k in range(81))
   assert left==right==int(i==j)
 results['actual_raw128_Pi1_head0_r_wedge_row']=[[j,str(c)] for j,c in sorted(actual.items())]
 results['all81_selected_original_unit_inverse_both_products']=True
 results['all780_X0_raw_row42_zero']=True
 results['selected_residual_first_r81']=[str(c) for c in sel]
 results['raw_domain_proof']='Direct n=E(mE)^-1 has row42 identically0. First24 cancellation on raw free44zero uses coefficients0; second18 unitrows lie in free20 and likewise uses0. Third81 head0 coupling0 kills inverse derivatives; rawPi1=-residualX1_row42*W0inv, giving +8/3 H27.'
 del X,xp,fac,wr,wp,res,inv,wall,ar,sel,Wsel
 assert get('quartic','actual_complete_quartic_dual1650_sparse')==[[329,'1']]
 L=json.loads((ROOT/FILES['Loewy']).read_text())
 assert L['exact_actual_HF']==[5,5,3,2,1] and L['dimensions_m_power_Tnil']==[16,11,6,3,1,0]
 results['true_relation_clearing_and_r_dual_proof']=[
  'Actual products c>=3; raw128 degree2=(0,H(m0*c3)_2). Max-prefix tensor is in known174 before actual3 quadratic pivots. Clear with TRUE all-order35 relations and linear coefficients u1.',
  'c3=LU*LV*(qWq+B Wbeta)+e43cross. LU,LV have no r component, so c42,r3=0. m0 e43=0 removes cross from leading H.',
  'The only rawPi1 contribution to the head0 r3 coefficient is +(8/3)H27,r2=0, since WHOLE sector0 last row is identically0. RawPi0 head0=e42 kills all higher wedge terms.',
  'TRUE relation M2 clearing contributes weighted fixed174 section = (50/3)head0,z0r+(125/3)head0,z4r=0: all raw degree2 minimal head0 coefficients are0. Relation-side -8/3 is not reused on raw128.',
  'After clearing the genuine product starts m3. Its r multiple starts m4; complete1650 dual is delta(head0 r4) and m5Tnil=0. Thus r times each U2,US,S2 is0. No Tnil algebra or unknown higher jet set0 is assumed.',
 ]
 results['actual_bar_consequences']={'Kbar_squared_subset_m3_Tnil':True,'r_Kbar_squared_zero':True,'QQ_dimension_Kbar_squared_upper':2}
 results['dimension_and_H0_scope']='Accepted q2 annihilator plus J-closure generates entire Kbar² over Rbar by U²,US,S², including r closure. In m3Tnil(length3), r image=m4(length1), so ker dimension2. H0/ellH0 injects Kbar (D-finite torsion saturation); Kbar injects Nbar since N/K embeds D-free normalizer nil. Since ell²K²=0, K² lies H0. Therefore its bar image dimension<=2 cannot D-generate H0 of generator number>=3. No H0 length/ell annihilator/original cube/index is inferred.'
 results['no_new_product_values_or_jets']=True
 bounds={'prior_parse174':parse174,'tensor_fraction_prior':tensor_prior,'true_clearing_fraction_prior':clear_prior,'raw_C1_fraction_prior':raw_prior,
         'phase_numeric_bounds':phase,'whole_numeric_upper':max(n for _,n in phase),'stored_H_observed':STORED,'pre_H_observed':PRE,
         'hard_numeric':150000,'hard_seconds':60,'hard_AS':2000000000,'RAW_input_bytes':sourcebytes}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 ev=out/'independent-actual-max-H324-raw-C1-and-r-annihilator-source.json'
 ev.write_text(json.dumps({'input_sha256':pins,'results':results,'bounds':bounds},indent=2)+'\n')
 rep={'id':ID,'all_passed':True,'status':'SOURCE_ACTUAL_MAX_H_RAW_PLUS_C1_AND_BAR_SQUARE_R_ANNIHILATOR',
      'input_sha256':pins,'code_sha256':sha(pathlib.Path(__file__)),'evidence_sha256':{ev.name:sha(ev)},
      'whole_numeric_upper':bounds['whole_numeric_upper'],'all324_actual_max_quotient6_zero':True,
      'raw_C1_sparse':[[27,'8/3']],'r_Kbar_squared_zero':True,'bar_square_QQ_dimension_upper':2,
      'source_pair_values_UNCOMPUTED':True,'stored_H_observed':STORED,'pre_H_observed':PRE,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
