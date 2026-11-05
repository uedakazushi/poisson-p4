#!/usr/bin/env python3
"""Canonical finite33→B17/C68 SOURCE gate only; no new G8 response values."""
from pathlib import Path
from fractions import Fraction as Q
import itertools,json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_20261105_v63'
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):
  (o/n).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in o.iterdir())<4000000
 def accepted(p,h=None):
  r=rd(p);assert r['all_passed']
  if h:assert pins[p]==h
  return r
 la='reports/AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_SOURCE_20261105_v63/';accepted(la+'report.json','d490c810468a9b5fc3011e762406e7342fbc8dbdc1978f9837e235664eb00c3b');L=rd(la+'actual-ten-nonzz-functionals-and-all660-source-NF-descent.json')
 accepted('reports/AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_HEIGHT_20261105_v63/report.json','fe9a9594f49cc6d97ddf0ebde9da5e954bea9c77ea4922a604702be0804d257a')
 accepted('reports/AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_SOURCE_20261105_v63_REVIEW/report.json','79a6ba0ff2ccc0190e8ac4d8a46515514349a75fa62428e28713a12a91a99ff6')
 ms='reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/';mr=accepted(ms+'report.json','cffac50257d126e9f0b4cae908228ec7d1b8ea2f547d5765cf331bdfd6633782')
 for p,h in mr['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in mr['evidence_sha256'].items():assert sha(ROOT/ms/p)==h
 accepted('reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json','92d5c943b12da089eb504862043a6e811ad17bef9d07e2b693e4eaf07a539b30')
 gauge=rd(ms+'phase2-canonical30-exact-values.json');native=rd(ms+'phase3-actual-kappa9-and-normalized-nilframe.json');monic=rd(ms+'phase1-source150-readback-fingerprints.json');kappa=native['kappa9_C132'];actual=native['actual14_C132'];nu=monic['nu_A_B33']
 assert len(kappa)==9 and all(len(v)==132 for v in kappa) and len(actual)==14 and all(len(v)==132 for v in actual) and len(nu)==2 and all(len(v)==33 for v in nu)
 assert all(Q(v[0])==0 for v in kappa+nu) and gauge['all30_constraints_exact'] and gauge['all30_gauge_B11_projections_exact'] and Q(gauge['fixed_J30_determinant'])==Q(4,9)
 interface='reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json';I=fld(interface,['retained_source9_indices','retained_source9_order','actual_group15_positions']);assert [actual[i] for i in I['retained_source9_indices']]==kappa
 g1=fld('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase3-gauge30-firsts.json',['even15_first9','odd15_first9']);ev=gauge['actual_even15_B33'];od=gauge['actual_odd15_B33'];assert len(ev)==len(od)==15
 assert [v[0] for v in ev]==[v[0] for v in g1['even15_first9']] and [v[0] for v in od]==[v[0] for v in g1['odd15_first9']] and all(Q(v[0])==0 for v in ev)
 J30=fld('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',['det_J30']);assert Q(J30['det_J30'])==Q(4,9)
 tablepath='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json';T=rd(tablepath);assert T['m5_zero'];basis33=[]
 for sec,r,e in T['basis33']:
  m=[0]*9;m[7]=r;m[8]=e
  if sec:m[sec-1]=1
  basis33.append(tuple(m))
 assert basis33==list(map(tuple,L['raw33_basis_monomials']))
 active=[0,1,4,5,7,8];pairs=list(map(tuple,L['pair_order10']));zero=(0,)*9;B=[zero]
 for i in active:m=[0]*9;m[i]=1;B.append(tuple(m))
 for i,j in pairs:m=[0]*9;m[i]+=1;m[j]+=1;B.append(tuple(m))
 assert len(B)==17 and len(set(B))==17;bi={m:i for i,m in enumerate(B)};projection=[({bi[m]:Q(1)} if m in bi else {}) for m in basis33];selection=[basis33.index(m) for m in B]
 for t0,row in enumerate(L['lambda10_by33_sparse_rows']):assert sp(row)=={selection[7+t0]:Q(1)}
 btable={}
 for i,m in enumerate(B):
  for j,n in enumerate(B):
   s=tuple(x+y for x,y in zip(m,n))
   if s in bi:btable[i,j]=bi[s]
 # Preoperation bound: projections have disjoint onehot supports. All ring/NF
 # checks copy or drop cached coefficients; coefficientproductsareby0/1only.
 NFp='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';F=fld(NFp,['monomials660_degree_le4','all660_actual_raw33_normal_forms'])
 sourcecount=count([L,mr,gauge,native,monic,I,g1,J30,T,basis33,B,projection,selection,btable,F,pins]);assert sourcecount+40000<150000
 ringprior=max([1]+[max(abs(Q(c).numerator).bit_length(),Q(c).denominator.bit_length()) for row in T['table_columns_for_each_left_basis'] for p in row for i,c in p]);assert ringprior<=384
 put('phase1-source-precheck.json',{'phase':'canonical33_input_before_known_ring_and_overlap_checks','source_owned_prior':sourcecount+40000,'ring_copy_prior_bits':ringprior,'no_G8_native_functionresponses':True})
 def proj(v):
  q={}
  for i,x in v.items():
   for j,c in projection[i].items():assert j not in q;q[j]=x
  return q
 for i in range(33):
  for j in range(33):
   lhs=proj(sp(T['table_columns_for_each_left_basis'][i][j]));a=projection[i];b=projection[j]
   rhs={} if not a or not b or (next(iter(a)),next(iter(b))) not in btable else {btable[next(iter(a)),next(iter(b))]:Q(1)};assert lhs==rhs
 for m,nf in zip(F['monomials660_degree_le4'],F['all660_actual_raw33_normal_forms']):assert proj(sp(nf))==({bi[tuple(m)]:Q(1)} if tuple(m) in bi else {})
 # EXACT cache-overlap maps; no new nonzzG8 response is evaluated.
 bp='reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json';b11=fld(bp,['B11_basis_monomials8','nu_A_B_actual_B11','gauge_source_convention']);b11idx=[basis33.index(tuple(m+[0])) for m in b11['B11_basis_monomials8']]
 oldb=fld('reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',['actual14_C_B44','kappa9_C_B44','actualG8_C_B44']);gb=fld('reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase2-gauge30-B11-correction.json',['even15_B11','odd15_B11'])
 def trace_select(v,ids):return [v[4*i+m] for i in ids for m in range(4)]
 assert all(trace_select(v,b11idx)==w for v,w in zip(kappa,oldb['kappa9_C_B44'])) and all(trace_select(v,b11idx)==w for v,w in zip(actual,oldb['actual14_C_B44']))
 assert all([v[i] for i in b11idx]==w for v,w in zip(ev,gb['even15_B11'])) and all([v[i] for i in b11idx]==w for v,w in zip(od,gb['odd15_B11']))
 for v,w in zip(nu,b11['nu_A_B_actual_B11']):assert {j:Q(v[i]) for j,i in enumerate(b11idx) if Q(v[i])}==sp(w)
 c16p='reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63/actual-native-Rell4-cache-source-input-and-adapter-gate.json';c16=fld(c16p,['Rell4_basis','raw33_to_Rell4_columns','actual_even15_Rell4','actual_odd15_Rell4','actual14_C16','actual_kappa9_C16','actual_nu_A_B_Rell4','normalizer_source_references'])
 reidx=[0,8,9,25];assert all(sp(c16['raw33_to_Rell4_columns'][i])=={j:Q(1)} for j,i in enumerate(reidx))
 assert all(trace_select(v,reidx)==w for v,w in zip(kappa,c16['actual_kappa9_C16'])) and all(trace_select(v,reidx)==w for v,w in zip(actual,c16['actual14_C16']))
 assert all([v[i] for i in reidx]==w for v,w in zip(ev,c16['actual_even15_Rell4'])) and all([v[i] for i in reidx]==w for v,w in zip(od,c16['actual_odd15_Rell4'])) and all([v[i] for i in reidx]==w for v,w in zip(nu,c16['actual_nu_A_B_Rell4']))
 for p,h in [('reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json','bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693'),('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63_REVIEW/report.json',None)]:accepted(p,h)
 K=rd('reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json');prep=rd('reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json');assert K['source9_variable_order']==I['retained_source9_order'] and prep['all7_unit_stripping_exact']
 # FormalG8polynomialsinκ,m^7truncated. P=5+d² inversebinomial4terms.
 # OnlysmallQQsourcecoefficientsarecombined,neverκfunctionvalues.
 roots=K['five_actual_rational_roots'];assert all(v['D_power']==0 and 0<=v['P_power']<=3 for v in roots.values());program=[];raw_terms=[];core=[0,1,2,6]
 for name in K['selected5_pivot_order']:
  p=roots[name];terms=[];e=[0]*9;e[K['source9_variable_order'].index(name)]=1;terms.append((tuple(e),Q(1)))
  for ex,n,d in p['numerator']:
   for k in range(4):
    if p['P_power']==0 and k:continue
    coeff=Q(n,d)*(Q((-1)**k*math.comb(p['P_power']+k-1,k),5**(p['P_power']+k)) if p['P_power'] else Q(1));e=[0]*9
    for i,c in zip(core,ex):e[i]=c
    e[1]+=2*k
    if sum(e)<=6:terms.append((tuple(e),-coeff))
  raw_terms.append(terms)
 for poly in prep['generators3']:
  terms=[]
  for ex,n,d in poly:
   e=[0]*9
   for i,c in zip(core,ex):e[i]=c
   if sum(e)<=6:terms.append((tuple(e),Q(n,d)))
  raw_terms.append(terms)
 formal_constD=math.lcm(*(q.denominator for terms in raw_terms for ex,q in terms));formal_constU=max(sum(abs(q.numerator)*(formal_constD//q.denominator) for ex,q in terms) for terms in raw_terms)
 assert max(formal_constD.bit_length(),formal_constU.bit_length())<128
 for terms in raw_terms:
  poly={}
  for e,q in terms:poly[e]=poly.get(e,Q(0))+q
  program.append([[list(e),str(q)] for e,q in sorted(poly.items()) if q])
 assert len(program)==8
 # SOURCE-only coefficient budget. Each nonzeroC68monomialuses<=6κfactors
 # and<=2νcarries;commoninputD^8 withformalcoeffD boundsallscalarterms.
 operands=[Q(v[4*i+m]) for v in kappa for i in selection for m in range(4)]+[Q(v[i]) for v in nu for i in selection];inputD=math.lcm(*(q.denominator for q in operands));inputU=max(abs(q.numerator)*(inputD//q.denominator) for q in operands)
 pathbits=100;coarse=max(8*inputD.bit_length()+formal_constD.bit_length(),8*max(inputD.bit_length(),inputU.bit_length())+max(formal_constD.bit_length(),formal_constU.bit_length())+pathbits);assert coarse<=384 and 3*coarse+1<=1153
 # Retain originalcachedarrays: futureprojectionandG8evaluationremainUNRUN.
 payload={'physical_R0_variables':L['physical_R0_variable_order'],'basis_B17':B,'raw33_basis':basis33,'raw33_to_B17_columns':[[[j,'1'] for j in v] for v in projection],'B17_raw33_selector_indices':selection,'B17_multiplication_nonzero_entries':[[i,j,k] for (i,j),k in sorted(btable.items())],'free4_basis_mask_order':['1','tauA','tauB','tauAtauB'],'C68_index_rule':'4*B17index+mask','actual_kappa9_C132_cached_UNEVALUATED':kappa,'actual_nu_A_B33_cached_UNEVALUATED':nu,'canonical30_branch':{'group15_positions':I['actual_group15_positions'],'same_even15_origin':[v[0] for v in ev],'same_odd15_origin':[v[0] for v in od],'fixed_J30_det':'4/9','aux_etaA_then_physical_tauA':'Canonicalscalar30solvedinbase[etaA]/etaA²=nuA;embedsamecoefficientsusingphysicaltauA. No etaA=tauA shortcutduringimplicitdifferentiation.'},'source_G8_order':K['selected5_pivot_order']+['Q0','Q1','s*r_core²'],'core_kappa_indices':core,'G8_formal_kappa9_polynomials_mod_radical7':program,'G8_formula_source':'first5κselected−exact5rationalroots(core);last3preparationgenerators3. OnlyP=5+d²units occur, expandedvia4termbinomial;nohonest/nonunitdivision.','G8_polynomial_truncation_proof':'B17m³=0;weightbasevar2,tau1;νinbase_m⇒carryweight>=2,maxbasisweight6. κresidue0⇒κ^7=0,(d²)^4=0. ThereforeP^-n=5^-n Σk0..3(-1)^k binom(n+k−1,k)d^(2k)/5^k. Droppingκdegree>=7exact.','response_inventory':{'nonzz_pairs10':pairs,'B17_pair_indices':list(range(7,17)),'G8_scalar_components':8*4,'G8_requested_second_QQ':320,'monic_requested_second_QQ':20,'point_active_first_overlap_QQ':8*4*7,'full_C68_G8_QQ_if_retained':544,'full_B17_monics_QQ_if_retained':34},'old_B11_projection_raw33_indices':b11idx,'old_B11_expected_G8_C44':oldb['actualG8_C_B44'],'old_C16_projection_raw33_indices':reidx,'old_C16_expected_G8_reference':{'path':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json','field':'actualG8'},'native_individual_pivot_L2S0_still_retained':True,'pureplane3_4_P9_is_SEPARATE':'B17killsm³;noG8r²ellorotherthirdresponseprovided.','genuine_honest_quotients_not_in_program':'No g/h/coordinatequotient callbackoccursinsideκ/G8unitformula. AbstractdescentdoesnotcomputeimagesofouterOLDquotients;theyremainA0OLD2frontier.','no_new_G8_kappa_monic_response_values_computed':True}
 numeric=count(payload);future=2*numeric+9*68+2*17+8*68*2+6000+4000+5000;assert future<60000
 opaque_path='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';opaque=fld(opaque_path,['typed_actual_C0_function_classification','newg4_second_symbols_not_zero_filled']);only79=[i for i,k in opaque['typed_actual_C0_function_classification'] if k=='C0_79_only'];assert len(only79)==36 and opaque['newg4_second_symbols_not_zero_filled']==[566400,566401,566402,566403]
 # Remainingwholecomparisonfrontiersareexplicit,notusedbythefutureG8worker.
 payload['outside_native_opaque_OLD2_source']={'path':opaque_path,'sha256':pins[opaque_path],'C0_79_only_ids_still_unresolved_or_order_gated_in_OTHER_lane':only79,'new_gA4_nonunit_quotient_ids':[566400,566401,566402,566403],'native_program_honest_dependencies':[]}
 sourceowned=count([L,mr,gauge,native,monic,kappa,actual,nu,I,g1,J30,T,basis33,B,projection,selection,btable,F,b11,oldb,gb,c16,K,prep,roots,raw_terms,program,payload,opaque,pins])+25000;assert sourceowned<150000
 rawinputs=sum((ROOT/p).stat().st_size for p in pins);assert rawinputs+4000000<33554432
 ep='canonical33-to-B17-C68-G8-source-input-and-response-budget.json';put(ep,payload)
 proof={'all1089_raw33_to_B17_ringproducts_and660_NFchecked':True,'all10_nonzzfunctionals_preserved':True,'canonical_actual14_to_kappa_retained_map_exact':True,'all30_origins_B11_Rell4_gauge_overlaps_exact':True,'all_kappa_actual14_B11_C16_overlaps_exact':True,'both_monic_B11_C16_overlaps_exact':True,'same_CA_F9_Q9_physicalR0_M0_frame_inherited':True,'only_fixed_unitIFT_cache_no_original105_or_gauge_solve':True,'G8_formal_unit_recipe_no_signed_response':True,'coarse_scalar_arithmetic_proof':'EachsurvivingC68producthas<=6κfactorsand<=2νcarries. Basisproducts0/1andatmost68^6*17^2*1000formalterms*4unitexpansionpaths<2^100. CommoninputdenominatorD^8andformalconstantdenominatorboundallpartialabsolute-sums.','input_LCM_bits':inputD.bit_length(),'input_scaled_absmax_bits':inputU.bit_length(),'formal_constant_LCM_bits':formal_constD.bit_length(),'formal_constant_abs_sum_scaled_bits':formal_constU.bit_length(),'native_future_coarse_committed_bound':coarse,'native_future_coarse_temporary_bound':3*coarse+1,'source_owned':sourceowned,'future_numeric_bound':future,'future_signed_response_worker_NOT_run':True,'native_valueGO_required_after_independent_source_peer_and_immutable_worker':True}
 pp='source-branch-descent-height-and-full-response-applicability-proof.json';put(pp,proof)
 report={'id':STEM,'all_passed':True,'SOURCE_cache_binding_descent_program_budget_only':True,'no_new_native_functionresponse_values':True,'B17_dimension':17,'C68_dimension':68,'requested_G8_second320_and_monic20':True,'all_cached_branch_monic_and_kappa_overlaps':True,'coarse_future_committed_bits':coarse,'coarse_future_temporary_bits':3*coarse+1,'future_numeric_upper':future,'source_owned_upper':sourceowned,'whole_SOURCE_RAW_upper':rawinputs+4000000,'future_value_worker_prepared_or_authorized':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep),pp:sha(o/pp)},'seconds':time.monotonic()-start};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'B17selection':selection,'coarse':coarse,'future':future,'owned':sourceowned}))
if __name__=='__main__':main()
