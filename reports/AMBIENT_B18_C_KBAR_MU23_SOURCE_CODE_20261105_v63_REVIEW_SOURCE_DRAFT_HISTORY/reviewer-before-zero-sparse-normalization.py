#!/usr/bin/env python3
"""Independent muK source/custody audit. NO requested rr/action dot."""
import ast,hashlib,importlib.util,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C_KBAR_MU23_SOURCE_CODE_20261105_v63_REVIEW'
SRC='reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63'
INV='reports/AMBIENT_B18_C_KBAR_MU23_SOURCE_INVENTORY_20261105_v63'
CTL='reports/AMBIENT_B18_C_KBAR_MU23_VALUE_CONTROL_V2_CUSTODY_20261105_v63'
WORKER='code/ambient_B18_C_Kbar_mu23_cached_values_v2_20261105_v63.py'
PINS={SRC+'/report.json':'332217ddd9dbdde76278f4b9d5e3f44240b69f86cab3821c75a91a699009fd7e',
 SRC+'/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json':'317b2c5e491e7f8369e73a7578982731415ac45d77444219308dcc5cc34b0127',
 WORKER:'94d50c2461dc1f9d77e62ccbf547acd9c97f0bcd93b3bc829a2829e4ab77b7e7',
 CTL+'/report.json':'c8aadfb38eac43bd8a7391ee35267fb3a8726ce339313fecbb2d67dc8addafcc',
 CTL+'/v2-control-only-custody-full-prefix-toy-fixtures-and-whole-budget.json':'c01453f6e65a9213b0b85cfc1c74fec3911f5e80ef8c36a108d02b3389dbf551',
 'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63_REVIEW/report.json':'662785338116da297cfe47ea8b8e117fee0b4ee2834d808458600f59f95384c7',
 'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63_REVIEW/report.json':'a5ec34be9b391a8dc9b6a34ce242c962046249301faab7c43deea42cc5cfe71a'}
PRE=STORED=0
def parse(c):
 q=Q(c);assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=384;return q
def op(a,b,m=False):
 global PRE,STORED
 if not a:return Q(0) if m else b
 if not b:return Q(0) if m else a
 na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
 h=max(na+nb,da+db) if m else max(max(na+db,nb+da)+1,da+db)
 PRE=max(PRE,h);assert h<=1153
 v=a*b if m else a+b;STORED=max(STORED,max(abs(v.numerator).bit_length(),v.denominator.bit_length()));assert STORED<=384;return v
def qa(a,b):return op(a,b)
def qm(a,b):return op(a,b,True)
def sm(vals):
 v=Q(0)
 for c in vals:v=qa(v,c)
 return v
def put(d,k,c):
 if c:
  v=qa(d.get(k,Q(0)),c)
  if v:d[k]=v
  elif k in d:del d[k]
def mon(*a):return tuple(a.count(i) for i in range(8))
def poly(v):return {(t,tuple(e)):parse(c) for t,e,c in v}
def enc(p):return [[t,list(e),str(c)] for (t,e),c in sorted(p.items())]
def col(v):
 a=[Q(0)]*40
 for i,c in v:a[i]=parse(c)
 return a
def pfrom(v):return {(i//8,mon(i%8)):c for i,c in enumerate(v) if c}
def pcomb(gs,co):
 d={}
 for c,g in zip(co,gs):
  if c:
   for k,v in g.items():put(d,k,qm(c,v))
 return d
def envelope(xs):
 xs=list(xs);D=math.lcm(*(c.denominator for c in xs));N=max(abs(c.numerator)*(D//c.denominator) for c in xs)
 return N,D
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 pins={};raw=0
 def pin(p,h=None,read=False):
  nonlocal raw
  p=pathlib.Path(p);path=p if p.is_absolute() else ROOT/p;k=str(path.relative_to(ROOT));s=sha(path)
  if h:assert s==h,(k,s,h)
  if k not in pins:raw+=path.stat().st_size
  pins[k]=s
  return json.loads(path.read_text()) if read else path
 for p,h in PINS.items():pin(p,h)
 sr=json.loads((ROOT/(SRC+'/report.json')).read_text());ir=pin(INV+'/report.json',read=True)
 for rep in (sr,ir):
  assert rep['all_passed']
  for p,h in rep['input_sha256'].items():pin(p,h)
  for p,h in rep['evidence_sha256'].items():pin((SRC if rep is sr else INV)+'/'+p,h)
 compiled=ROOT/(SRC+'/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json')
 D=json.loads(compiled.read_text())
 inventory=ROOT/(INV+'/fourteen-z-actions-common37-M2-source-functionals-and-minimal-generator-criterion.json')
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 def f(s,k):return field(src(s),k)
 comp=pin(SRC+'/all37-raw-X0-zero-and-C0-cache-hash-source-companion.json',read=True)
 for p,h in comp['input_sha256'].items():pin(p,h)
 # Source37 origin phase is separate from all later decoded arrays.
 combos=f('actual119-cubic-X2-only-source-operator-preflight.json','actual_selected37_original768_constant_combinations')
 X0=f('actual-X0-and-unit-cancellation.json','full_X0_56_by768')
 owned_origin=count((D,comp,combos,X0,sr,ir,pins))+24000;assert owned_origin<150000
 for co in combos:
  assert all(0<=j<768 for j,c in co)
  for row in X0:assert sm(qm(parse(c),parse(row[j])) for j,c in reversed(co))==0
 del X0
 # Independent fixed linear system; no producer monomial division.
 cert=f('actual-linear-quadratic-T-certificate.json','complete_linear_rank_certificate')
 chosen=cert['selected_original_columns'];assert len(chosen)==37
 norms=list(map(col,cert['basis_normalized_columns']))
 actual=f('phase1-actual-five-target-linear-relations.json','actual40_by657_linear_relation_columns')
 original=[col(actual[j]) for j in chosen];del actual
 gs=list(map(poly,f('known-initial40-and-length-upper17.json','homogeneous_module_generators40')))
 assert [pfrom(v) for v in norms]==gs[:37]
 assert f('known-initial40-and-length-upper17.json','standard_degree1')==['e3*r','e4*z5','e4*r']
 assert f('known-initial40-and-length-upper17.json','standard_degree2')==['e4*z5*r','e4*r²']
 assert f('exact-Loewy11-and-source-basis-corollary.json','exact_Loewy_HF')==[5,3,2,1]
 T=[[parse(c) for c in row] for row in field(inventory,'normalized37_actual37_constant_lift')]
 assert len(T)==37 and all(len(row)==37 for row in T)
 for i,co in enumerate(T):
  for j in range(40):assert sm(qm(c,original[k][j]) for k,c in reversed(list(enumerate(co))))==norms[i][j]
 coeff=field(src('actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json'),'actual_correction_coefficients')
 aU,cU,aS,cS=[parse(coeff[k]) for k in ('a3=b3/b1','c3=(p3-a3*p1)/d','aS=(b2-b3/4)/b1','cS=(p2-p3/4-aS*p1)/d')]
 heads=[{3:Q(1),1:-aU},{2:Q(1),3:Q(-1,4),1:-aS}]
 vectors=[]
 for v in heads:
  a=[Q(0)]*40
  for t,c in v.items():a[8*t+7]=c
  vectors.append(a)
 for v in heads:
  for z in range(7):
   a=[Q(0)]*40
   for t,c in v.items():a[8*t+z]=c
   vectors.append(a)
 piv=[i for i in range(40) if i not in (31,37,39)]
 A=[[norms[j][i] for j in range(37)]+[v[i] for v in vectors] for i in piv]
 # Bottommost pivot order, different from producer's leading-monomial reduction.
 for j in range(37):
  ii=next(i for i in reversed(range(j,37)) if A[i][j]);A[j],A[ii]=A[ii],A[j]
  inv=1/A[j][j];A[j]=[qm(q,inv) for q in A[j]]
  for i in reversed(range(37)):
   if i!=j and A[i][j]:
    c=A[i][j];A[i]=[qa(x,-qm(c,y)) for x,y in zip(A[i],A[j])]
 cofs=[[A[j][37+k] for j in range(37)] for k in range(16)];rem=[]
 for v,cf in zip(vectors,cofs):
  r=[qa(x,-sm(qm(c,norms[j][i]) for j,c in enumerate(cf))) for i,x in enumerate(v)]
  assert all(not x for i,x in enumerate(r) if i not in (31,37,39));rem.append(r)
 assert rem[0][31]==1 and rem[0][39]==0 and rem[1][31]==rem[1][39]==0
 lam=rem[1][37];assert lam
 assert field(inventory,'known_rU_rS_firsts')=={'A':str(rem[0][37]),'lambdaS':str(lam)}
 assert all(r==[0]*40 for r in rem[2:])
 proofs=field(inventory,'all14_z_U_S_full_source_and_coherent_M2_functionals')
 forms=[]
 for k,cf in enumerate(cofs[2:]):
  weights=[sm(qm(cf[j],T[j][i]) for j in reversed(range(37))) for i in range(37)]
  assert [str(c) for c in weights]==proofs[k]['actual_selected37_M2_weights']
  assert pfrom(vectors[2+k])==poly(proofs[k]['degree1_original_source'])
  assert pcomb([pfrom(v) for v in original],weights)==pfrom(vectors[2+k])
  assert proofs[k]['source']==[k%7,'U' if k<7 else 'S'] and proofs[k]['first_NF_v_w']==['0','0']
  z=k%7;cx=cU if k<7 else cS;F2={(4,mon(z,7)):-cx}
  assert poly(proofs[k]['degree2_explicit_source'])==F2
  back=poly(proofs[k]['degree2_known_NF'])
  for j,sh,c in reversed(proofs[k]['quadratic_initial_membership']):
   for (t,e),v in gs[j].items():put(back,(t,tuple(a+b for a,b in zip(e,sh))),qm(parse(c),v))
  assert back==F2
  forms.append([Q(0)]+[-c for c in weights])
 assert [[str(c) for c in row] for row in forms]==field(inventory,'all14_affine_rows_constant_then37_unknown_t')
 assert field(inventory,'actual_pi_q2_degree2_sparse')==[[179,'1']]
 # Literal covector179 kills every degree2 first/quad relation.
 for g in gs:
  deg=sum(next(iter(g))[1])
  if deg>2:continue
  for sh in ([mon(i) for i in range(8)] if deg==1 else [(0,)*8]):
   assert not any(t==4 and tuple(a+b for a,b in zip(e,sh))==mon(7,7) and c for (t,e),c in g.items())
 source_table=field(inventory,'actual_selected37_sources')
 assert source_table==[[j,168+j,c] for j,c in zip(chosen,combos)]
 owned_linear=count((D,comp,combos,cert,norms,original,gs,T,heads,vectors,A,cofs,rem,proofs,forms,source_table,pins))+26000
 assert owned_linear<150000
 del cert,norms,original,gs,heads,vectors,A,cofs,rem,proofs,source_table
 # Full literal last-row frame and first-order Schur carries.
 unit=src('actual-fullV56-m0-and-typed-comparison.json')
 assert field(unit,'typed_n56_by24')[55]==[0]*24
 assert field(unit,'typed_regular_function_operators')[0]==['zero']
 schur=src('actual-minimal5-by825-Schur-operator.json')
 assert field(schur,'actual_free_five_source_rows')==[22,23,29,30,55]
 assert field(schur,'constant_free32_quotient_projection5_rows')[4]==[[31,'1']]
 assert field(schur,'actual_relation_order')==['168 literal d3 images','657 remaining actual source cycles']
 first=f('actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9')
 assert len(first)==640
 for v in first:
  for j,row in v:
   if j==55:assert not parse(dict(row).get(8,'0'))
 del first
 xpr=f('actual-projected-three-class-X2-operator-preflight.json','actual_pivot_X1_by27_sparse56_times8')
 assert len(xpr)==27 and all(not parse(dict(v).get(447,'0')) for v in xpr)
 ys=f('phase1-actual-five-target-linear-relations.json','normalized84_unitcycle_Q1_columns_sparse')
 assert len(ys)==84 and all(not parse(dict(v).get(39,'0')) for v in ys)
 assert f('actual119-cubic-X2-only-source-operator-preflight.json','all37_m0_Xj1_zero_and_pi_Xj1_equal_actual_Lj')
 P=f('honest-C0-D0-D4-DAG.json','honest_D0_relation_columns');assert len(P)==640 and all(len(c)==9 for c in P)
 w37=[]
 for c in combos:
  d={}
  for j,q in reversed(c):
   if j<640 and P[j][8]:put(d,P[j][8],parse(q))
  w37.append(d)
 assert [[[j,str(c)] for j,c in sorted(v.items())] for v in w37]==D['all37_Pg_source_weights']
 roots=sorted({j for v in w37 for j in v});assert len(roots)==21 and sum(map(len,w37))==98
 w14=[]
 for row in forms:
  d={}
  for q,v in zip(row[1:],w37):
   for j,c in v.items():put(d,j,qm(q,c))
  w14.append(d)
 assert [[[j,str(c)] for j,c in sorted(v.items())] for v in w14]==D['all14_direct_Pg_source_weights']
 assert [[str(c) for c in row[1:]] for row in forms]==D['all14_same37_weights']
 # Compare WHOLE provided polynomials, never select the requested rr leaf.
 cache=f('phase3-oldB-common-source2-values-and-symbols.json','C0_composite_values_and_opaque_second_gradients')
 supplied={i:(v,g) for i,v,g in cache if i in roots};del cache
 assert len(supplied)==21 and all(g==[] for v,g in supplied.values())
 assert [[i,supplied[i][0]] for i in roots]==D['all21_given_scalar_cache_rows_not_selected_rr']
 mon20=f('phase3-oldB-common-source2-values-and-symbols.json','base20_monomials_original8')
 assert mon20==D['actual_cache_monomials_original8'] and len(mon20)==20 and mon20[19]==[7,7] and mon20[5]==[7]
 cached_report=json.loads(src('_OLDB_AD/report.json').read_text())
 C0path=str(src('honest-C0-D0-D4-DAG.json').relative_to(ROOT))
 assert cached_report['input_sha256'][C0path]==pins[C0path]
 # Positive envelope of all provided entries; no source/function response dot.
 Ec=envelope(parse(c) for v,g in supplied.values() for k,c in v)
 E37=envelope(c for v in w37 for c in v.values());E14=envelope(c for v in w14 for c in v.values());Ef=envelope(c for row in forms for c in row[1:])
 bound=lambda a,b,n:(n*a[0]*b[0],a[1]*b[1])
 Et=bound(E37,Ec,21);Ea=bound(E14,Ec,21);E2=bound(Ef,Et,37);El=(abs(lam.numerator),lam.denominator)
 Em=(Ea[0]*El[0]**2,Ea[1]*El[1]**2)
 es={'cache':Ec,'weights37':E37,'weights14':E14,'weights14_by37':Ef,'all37t':Et,'all14direct':Ea,'all14two_stage':E2,'witness_unit_minor':Em}
 H=max(max(abs(n).bit_length(),d.bit_length()) for n,d in es.values());assert H==136
 assert all(list(v)==D['positive_whole_future_height'][k] for k,v in es.items())
 assert D['positive_whole_future_height']['committed']==H and D['positive_whole_future_height']['temporary']==2*H+1==273
 owned_cache=count((D,comp,combos,T,forms,xpr,ys,P,w37,w14,supplied,mon20,pins))+27000;assert owned_cache<150000
 # Literal worker/control audit and unrelated synthetic guard test only.
 worker=ROOT/WORKER;text=worker.read_text();tree=ast.parse(text)
 assert {n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)}=={'fractions','pathlib'}
 assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ('eval','exec','compile') for n in ast.walk(tree))
 assert text.index("SCI['phase']='mandatory_lower_readbacks'")<text.index("SCI['phase']='raw37_rr'")<text.index("SCI['phase']='fourteen_action_folds'")
 for s in ["SCI['pending']['computed_candidate']=str(q)","SCI['tterms']","except BaseException as exc:","failure_prefix(exc);raise","assert (COMMITTED,TEMPORARY)==(136,273)","assert count(peer)<=1000","SCI['rank_result']","terms.append([int(root),str(w),str(v),str(p)]);z=add(z,p)"]:
  assert s in text,s
 spec=importlib.util.spec_from_file_location('mu_source_guard_fixture_only',worker);w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
 assert w.COMPILED_SHA==pins[str(compiled.relative_to(ROOT))]
 assert w.count(D)==1047 and D['future_numeric_profile']['static_owned_upper']==14076
 # Both synthetic tests are unrelated to cache roots and never invoke main/dot.
 assert w.mul(Q(2,7),Q(3,5))==Q(6,35)
 old=w.COMMITTED;w.COMMITTED=3
 try:w.mul(Q(3),Q(3));raise AssertionError('missing candidate guard')
 except AssertionError as e:assert w.SCI['pending']['computed_candidate']=='9' and not w.SCI['pending']['completed']
 w.COMMITTED=old
 ctl=json.loads((ROOT/(CTL+'/v2-control-only-custody-full-prefix-toy-fixtures-and-whole-budget.json')).read_text())
 addendum=pin(CTL+'/controller-candidate-inclusive-conservative-count-companion.json',read=True)
 assert addendum['SCI_numeric_upper']==2200 and addendum['save_numeric_upper']==8100 and addendum['failure_numeric_upper']==10200
 assert 10200<14076<60000 and ctl['shape_only_RAW_proof']['whole_upper']==510000<800000
 toy=pin(CTL+'/toy-prefix-fixture/STOP-scientific-prefix.json',read=True)
 assert toy['status']=='STOP_SCIENTIFIC_PREFIX_NO_AUTOMATIC_RETRY' and toy['scientific_state']['pending']['operation']=='parse'
 assert not (ROOT/'reports/AMBIENT_B18_C_KBAR_MU23_CACHED_VALUES_20261105_v63').exists()
 assert raw+400000<32*1024**2
 evidence={'scope':'SOURCE/code/custody only; NO selected rr coefficient,37t,14action or minor value evaluated.',
  'input_sha256':pins,'independent_fixed_linear_system_16_RHS':True,'all1480_true37_lift_multiplybacks':True,
  'all37_full56_X0_source_combinations_zero':True,'all14_z_first_normal_forms_zero':True,'q2_projection_literal_index179':True,
  'all640_raw55_r_all27_Xp55r_all84_Ys1_head4r_zero':True,'all21_literal_P9_roots_same_OLD640_namespace':roots,
  'all37_98_root_weights_and14_two_stage_source_weights_match':True,'all21_full_provided_polynomials_match_and_opaque_gradients_empty':True,
  'same_actual_C0_hash_in_old_AD_cache':True,'positive_envelopes':{k:list(v) for k,v in es.items()},'future_H136_pre273':True,
  'mu_criterion_proof':'Known rank(r|Kbar)=4 gives Kbar/rKbar basis U,S,q2. All14 zU/S firstnormalforms vanish; their class lies in oldm² and modulo rK only q2 survives. Every ziq2 is in oldm³=kg with g=rq2. The fourteen scalarfunctionals therefore exhaust m action on this3fiber: some nonzero gives mu2, allzero mu3. Finite ell-Nakayama transfers to SAME original K over R0, not Dmu or H0.',
  'true_M2_proof':'Constant true37 OLD768 combinations have full56 X0=0. Direct section nrow55=0 identically; last Pi0=e55. FirstSchur raw55_r rows (640 and27) and secondSchur Ys1head4r allzero, so requested true M2 lasthead rr is literal OLD P9_g second coefficient. Full C0 nonlinear CB transport is retained in accepted all21 cachepolynomials; no homogeneous M2 substitution or unprovided higher extension zero.',
  'worker_V2_math_control_review':'Exact compiled pin, mandatory37 origin/r first readbacks precede rr; all14 direct21root vs via37 compare; witness lambdaS²*action onlyif actual nonzero; phase1 saved beforephase2. SCI pending operation/operands/computedcandidate/partial terms/cursors survive all BaseException and wallTimeout, no retry. V1 unrun bytes retained.',
  'synthetic_only_candidate_guard_passed':True,'candidate_inclusive_SCI2200_save8100_failure10200_below_static14076':True,
  'whole_source_numeric_upper':max(owned_origin,owned_linear,owned_cache),'SOURCE_phase_bounds':[owned_origin,owned_linear,owned_cache],
  'VALUE_static14076_hard60000_RAW510000_hard800000':True,'RAW_source_with_outputreserve':raw+400000,
  'actual_SOURCE_stored':STORED,'actual_SOURCE_pre':PRE,'no_actual_rr_or_action_dot':True}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);ep=out/'independent-true37-source-lastrow-cache-and-V2-code-custody.json'
 ep.write_text(json.dumps(evidence,indent=2)+'\n')
 report={'id':ID,'all_passed':True,'SOURCE_code_custody_only':True,'input_sha256':pins,'code_sha256':sha(pathlib.Path(__file__)),
  'evidence_sha256':{ep.name:sha(ep)},'actual_requested_values_UNCOMPUTED':True,'whole_source_numeric_upper':evidence['whole_source_numeric_upper'],
  'worker_code_sha256':PINS[WORKER],'value_static_owned':14076,'future_H136_pre273':True,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'whole_source':evidence['whole_source_numeric_upper'],'seconds':report['seconds']}))
if __name__=='__main__':main()
