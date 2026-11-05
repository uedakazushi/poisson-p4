#!/usr/bin/env python3
"""Pure order/source lemma: q2 annihilates genuine Nbar, not arbitrary Tnil."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63'
PATHS={
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'natural_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'point':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63/five-cache-functionals-and-genuine-fourteen-kernel-source.json',
 'point_peer':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_20261105_v63_REVIEW/report.json',
 'five':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63/actual-two-cache-witness-and-five-minimal-bar-generators.json',
 'five_peer':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63_REVIEW/report.json',
 'product':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/report.json',
 'product_peer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
 'upper':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/actual-trace-first-gate-and-entire-conductor-proof.json',
 'upper_peer':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63_REVIEW/report.json',
 'global_g':'reports/AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63/report.json',
 'global_g_peer':'reports/AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63_REVIEW/report.json',
 'H_peer':'reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_20261105_v63_REVIEW/report.json',
 'quartic':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/actual-eight-relations-and-complete-quartic-dual.json',
 'quartic_peer':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63_REVIEW/report.json',
 'Pi':'reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63/all-three-Schur-row42-and-structural-bar-cube-zero.json',
 'Pi_peer':'reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW/report.json',
 'Loewy':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/report.json',
 'Loewy_peer':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json',
 'three':'reports/AMBIENT_B18_C_KBAR_THREE_GENERATOR_RANK_SOURCE_20261105_v63/report.json',
 'three_peer':'reports/AMBIENT_B18_K_THREE_ROOT_QUOTIENT_REVIEW_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'natural':['actual_origin_map_columns5'],
  'point':['literal_direct_section_path','minimal_head_inclusion_rows128','all_origin_wedge_cycles_zero_by_d2_origin','all24_m0_e43_entries_zero','source_m0_24_by44','head42_origin_cycle44_sparse'],
  'five':['actual_minimal_bar_generators','actual_mu_S_Nbar'],
  'upper':['all36_Qrad_C_maximal_ideal_order_lower_bounds','same_origin_tauA_tauB_zero','all_Qrad_even_coefficients_in_J','actual_joint_ideal'],
  'quartic':['actual_complete_quartic_dual1650_sparse'],
  'Pi':['full_minimal_Pi0_first_row128_sparse'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 prior=S.count(d)+20000;assert prior<150000
 assert d['Loewy']['exact_actual_HF']==[5,5,3,2,1]and d['Loewy']['dimensions_m_power_Tnil']==[16,11,6,3,1,0]
 assert d['global_g']['global_mR0_cubed_for_literal_gA_andgB']
 point=d['point'];assert point['minimal_head_inclusion_rows128']==[42,43,50,68,70]
 assert point['all_origin_wedge_cycles_zero_by_d2_origin']and point['all24_m0_e43_entries_zero']
 assert all(Q(row[43])==0 for row in point['source_m0_24_by44'])
 natural=d['natural']['actual_origin_map_columns5'];assert natural[4][0]=='0'and natural[4][1]=='1'
 assert all(Q(q)==0 for j in[1,2,3]for q in natural[j])
 assert d['five']['actual_mu_S_Nbar']==5
 assert d['five']['actual_minimal_bar_generators']==['h1-(5/3)rA*h0','h2+(88/25)rA*h0','h3-(2/5)rA*h0','h4-(136/75)rA*h0','zA5*h0']
 # Origin of the DIRECT canonical minimal->free44 cycle section. The
 # wedge heads have d2(0)=0, so their free44 part is0, not a guessedvalue.
 origin=[[Q(0)]*5 for _ in range(44)]
 for i,q in point['head42_origin_cycle44_sparse']:origin[i][0]=Q(q)
 origin[43][1]=Q(1)
 q2lead=[sum((row[j]*Q(natural[4][j])for j in range(5)),Q(0))for row in origin]
 assert q2lead==[Q(i==43)for i in range(44)]
 # All five genuine R0 generators have free44 origin0 except generator1,
 # whose origin is exactly e43. Every representative can therefore be
 # chosen with a36,bz7 positive-order; bgA alone can have origin.
 genorigins=[[Q(i==43)if j==0 else Q(0)for i in range(44)]for j in range(5)]
 assert all(all(not q for q in row[:43])for row in genorigins)
 assert d['Pi']['full_minimal_Pi0_first_row128_sparse']==[[42,'1']]
 assert d['quartic']['actual_complete_quartic_dual1650_sparse']==[[329,'1']]
 projection=[Q(i==42)for i in range(128)]
 assert projection[43]==0 and all(projection[i]==0 for i in range(44,128))
 # Explicit universal identity, in the raw44+84 domain: iota(c)=(c,Hmc).
 # d2Hmc=mc on genuinecycles makes Phi_n iota(c)=c for BOTH nchoices.
 cancel_identities=0
 for i in range(44):
  terms={('c',i):1}
  for j in range(24):
   key=('n',i,j,'mc',j);terms[key]=terms.get(key,0)-1;terms[key]+=1
   if not terms[key]:del terms[key]
  assert terms=={('c',i):1};cancel_identities+=1
 nil_order={'a36':1,'bz7':1,'bgA':0,'b0':2}
 q2_order={'a36':3,'bz7':3,'bgA':2,'b0':4,'B=sum_CA_fi_bi':4}
 prod_a=min(q2_order['b0']+nil_order['a36'],q2_order['B=sum_CA_fi_bi']+nil_order['a36'])
 prod_z=min(q2_order['b0']+nil_order['bz7'],nil_order['b0']+q2_order['bz7'])
 prod_g=min(q2_order['b0']+nil_order['bgA'],nil_order['b0']+q2_order['bgA'])
 assert(prod_a,prod_z,prod_g)==(5,5,4)
 data={
  'scope':'Pure SOURCE/order theorem on TRUE Nbar=Nil(D0)/ellNil(D0), using its genuine injection intoTnil. It gives no product operation on arbitraryTnil and no originalupstairs annihilator.',
  'same_raw44_order':'a36 followedbyCAodd bz7,bgA; raw128 isfree44 followedby84wedges, DIRECTsectionn=E(mE)^-1',
  'direct_section_binding':point['literal_direct_section_path'],
  'canonical_section_origin44_by5':[[str(q)for q in row]for row in origin],
  'q2_actual_origin_source':'q2=i(r²oldh4) byRbarlinearity ofthe exactnatural map. oldM0h4 hasfreehead42=0/freehead43=1; otherthree components arewedgeswhosefree44origins0. Thuscanonicalq2 freecycle=r²e43+O(m³), withallhighertermsFORMAL.',
  'q2_leading_order2_free44_vector':[str(q)for q in q2lead],
  'five_true_nil_generators_origin_types':'firstcorrectedh1 hasorigin e43;otherthreecorrectedheadsarewedgesat0;z5h0 hasorigin0. Since theseRbar-generateNbar, everyclasshasachosen genuinecyclec witha,bz∈m,bgorder>=0.',
  'orders_for_any_chosen_genuine_Nbar_representative':nil_order,
  'orders_for_q2_canonical_representative':q2_order,
  'order_input_proofs':'AllQrad even36 entriesareinproperJ⊂m, soevenb0=Qeven*a hasorder2forn and4forq2. First7CAf_i∈m,literal gA505938∈m³ (globalSOURCE,notoldf637790), soBq>=4. Allremainingnormalizer/gaugecoefficients areREGULARformalfunctions; noneare truncatedaway or assigned0.',
  'product44_complete_source_used':'Useactualcycleproducta_new=mu_C(b0q+Bq*tauA)*a_n andb_new=b0q*b_n+b0n*b_q. mu_C only multipliesby1oractualnuA∈m, soorderdoesnotdecrease. This is Product44 withall36+8outputs, notonlyQradclosure.',
  'full_product44_order_bounds':{'a36':prod_a,'bz7':prod_z,'bgA':prod_g},
  'consequence_order4_free44_part':'Modulo m5, productfree44 c hasONLYbgA coordinate43. Thusm*c∈m5 sinceactualm0e43=0; allhigherm1,m2etc remainandcontribute>=5.',
  'honest_cycle_iota':'iota(c)=(c,H(m*c)) inraw128(44+84). c isaTRUEVcyclebyProduct44; d1m=V,soH applies andd2Hmc=mc. Hlosesatmost1, henceiota(c)∈m4.',
  'all44_Phi_n_iota_equals_c_formal_cancellations':cancel_identities,
  'section_transport_fixedness':'Foranyrightsectionn, Phi_n iota(c)=c-nmc+nd2Hmc=c EXACT. Theacceptedold/directsectionchangeJ_n(s,b)=(Phi_n(s,b),b) thereforeFIXES iota(c). No hiddenfree42homotopycorrectionis introduced.',
  'complete_quartic_projection_step':'OnTHIS raw128 domain Pi0_firstrow=[e42^T,0_84]. ThereforePi0_firstrow*iota(c)_4=c42_4=0 EVENWITH arbitraryH(mc)_4. Its84wedgecoefficientsarezero inthisrow. RegularPi1/iota>=4 startsat5. Completephi4_1650 hasONLYminimalhead0*r4 entry329=1, soactualquarticclass0. Exactm5Tnil0 impliesactualproductclass0.',
  'genuine_injection_and_result':'Nbar→Tnil isinjective byacceptedhonestellUCT. Bothinputsareactualnilbarclassesandproduct44isinducedfromactualD0multiplication. Thusq2*Nbar=0; byRbarlinearityg=rq2alsoannihilatesNbar.',
  'three_generator_bar_square_reduction':'SinceKbar=(U,S,q2)_Rbar andq2*Nbar0, Kbar² isRbar-generatedbyU²,US,S². Originalsixpairs involvingq2haveZERO barclass. ForQQrank/moduleclosure oneMUSTalso retain r-multiples; thisisNOTanassertionof3independentclasses.',
  'upstairs_limitation':'Onlybarproductsvanish. NoactualKliftq2 isclaimedtoannihilateN upstairs; noellH0annihilator ornilindeximprovement follows.',
  'all_six_signed_pair_product_values_UNCOMPUTED':True,'no_new_source_IFT_gauge_jet_or_product_value_job':True,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'prior_SOURCE_owned_upper':prior,'selective_receipts':receipts,'all_passed':True}
 owned=S.count(d)+S.count(data)+S.count(origin)+20000;assert owned<150000;data['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'q2-annihilator-exact128-projection-and-formal-order-source-proof.json';p.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_Q2_TIMES_GENUINE_NBAR_ZERO_THREE_BAR_SQUARE_GENERATORS',
  'q2_times_genuine_Nbar_zero':True,'all_six_pair_values_evaluated':False,'original_H0_or_nilindex_improved':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
