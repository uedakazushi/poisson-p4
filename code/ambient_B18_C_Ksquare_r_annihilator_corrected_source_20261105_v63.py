#!/usr/bin/env python3
"""Corrected pure SOURCE proof: raw+8/3, H27=0; no product-value job."""
import hashlib,json,signal,time
from pathlib import Path
from fractions import Fraction as Q
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KSQUARE_R_ANNIHILATOR_CORRECTED_SOURCE_20261105_v63'
PATHS={
 'leading':'reports/AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63/unpaired-leading-three-pair-extraction-tensor-and-rank2-sufficient-criterion.json',
 'linear':'reports/AMBIENT_B18_C_KSQUARE_LEADING_LINEAR_SPAN_SOURCE_20261105_v63/unpaired-known174-membership-and-pure-r-H2-source.json',
 'old_rdual':'reports/AMBIENT_B18_C_KSQUARE_R_PAIR_DUAL_SOURCE_20261105_v63/report.json',
 'hold':'reports/AMBIENT_B18_C_KSQUARE_R_PAIR_DUAL_SOURCE_20261105_v63/RAW_C1_SIGN_HOLD_ERRATUM.json',
 'raw_domain':'reports/AMBIENT_B18_C_KSQUARE_R_DUAL_C1_DOMAIN_DIAGNOSTIC_20261105_v63/literal-raw128-domain-and-all81-origin-unit-derivative-reconstruction.json',
 'compact':'reports/AMBIENT_B18_C_KSQUARE_R_PAIR_DUAL_SOURCE_20261105_v63/compact-two-term-linear-clearing-and-sector0-row7-zero-factorization.json',
 'zero3':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json',
 'q2_peer':'reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63_REVIEW/report.json',
 'quartic':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/actual-eight-relations-and-complete-quartic-dual.json',
 'Loewy':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/report.json',
 'Loewy_peer':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json',
 'C1_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fs={'leading':['all324_unpaired_extraction_coefficients'],'linear':['all324_quotient6_zero'],
 'raw_domain':['actual_raw128_wedge84_Pi1_head0_r_sparse','first24_domain_proof','second18_domain_proof','third81_domain_proof'],
 'compact':['weighted_M2_head0r2_times_fixed174_r_section','m0_row7_Wq_and_Wbeta_both_zero','all_raw_degree2_head0_coefficients_literal_zero'],
 'zero3':['actual_three_minimal_first_columns40'],
 'quartic':['actual_complete_quartic_dual1650_sparse']}
 d={k:S.selective(ROOT/p,fs[k],receipts)if k in fs else json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k,v in d.items():
  if k.endswith('_peer'):assert v['all_passed']
 assert all(v['rows3']==['0','0','0']for v in d['leading']['all324_unpaired_extraction_coefficients'])
 assert d['linear']['all324_quotient6_zero']
 assert d['raw_domain']['actual_raw128_wedge84_Pi1_head0_r_sparse']==[[27,'8/3']]
 assert d['compact']['weighted_M2_head0r2_times_fixed174_r_section']==[[7,'50/3'],[29,'125/3']]
 assert d['compact']['m0_row7_Wq_and_Wbeta_both_zero']==['0','0']
 assert all(not Q(col.get('7','0'))for col in d['zero3']['actual_three_minimal_first_columns40'])
 assert d['quartic']['actual_complete_quartic_dual1650_sparse']==[[329,'1']]
 assert d['Loewy']['dimensions_m_power_Tnil']==[16,11,6,3,1,0]
 evidence={
 'status':'CORRECTED_PURE_SOURCE_NO_PRODUCT_PAIRINGS_PENDING_INDEPENDENT_ADOPTION','input_sha256':pins,
 'old_rdual_formula_not_adopted':'edeef usedrelation-side−8/3 asrawΠ1. PreservedHOLD6855. Theproofhereusesactualraw+8/3 explicitlyanddoesNOTinferfromoldsignincorrecttensoralone.',
 'correct_raw_domain_C1':[[27,'8/3']],
 'actual_input_domain':'iota2=(free44ZERO,H2), directnrow42identically0; first24/18unitcancellationpreservesthiswedgeinput; third81unitderivativegivesliteralrawΠ1above. Allunknownunitderivatives multiplyoriginfree42coupling0.',
 'whole_H_input_specific_zero':'Productc3=L_U L_V(qWq+B Wbeta)+e43cross. m0e43=0 andm0row7Wq=m0row7Wbeta=0 implyWHOLE sector0v7_3=0. HenceH27_2=0 beforeLU/LS/qpairing, signindependent. ThisisactualsharedH, notR0linearcolumnextension.',
 'full_true_relation_clearing_zero':'AllleadingpairsinKNOWN174first35*S1. Choosefixed174sectionu1. ClearusingTRUEallorderM*u1, nothomogeneousM1alone. Pure-r cubic head0clearingsum_j q35_j*(u1_j)_r isonly(50/3)*w2_head0,z0r+(125/3)*w2_head0,z4r=0 becausePi0row0=e42andrawfree2zero. NoactualquadM3liftneeded.',
 'raw42_cubic_zero':'LUandLSrcoeff0, c42_3=2LU*LV*q*v42_b6; pure-r³0. e43crosscontributesnoraw42.',
 'higher_coherent_dependencies_retained_but_absent_for_this_dual':'Gamma2,section1,m1maychangeH3/wedgeorrawc4. Pi0head0killsALLwedgecomponents;Pi1onH3isdegree4,Pi1onrawc3is4,Pi2onH2is4. Thereforethispure-r cubic head0dualdependsonlyactualc3andPi1_r H2_rr; higherdatawereNOTassigned0.',
 'r_pair_class_zero_proof':'Aftertruefirst35clearingrepresentativestartsdegree3. Itsrmultiple starts4, andcompletephi4readshead0r4. AllthreecomputedSOURCEtermsaboveare0, soactualrproductclass0 byexactm5Tnil0. GenuineNbar→Tnilinjects; actualproductsfromD0, noTnilalgebraassumed.',
 'result_on_true_bar_objects':{'r_times_U_squared':0,'r_times_U_S':0,'r_times_S_squared':0,'r_times_whole_Kbar_squared':0,'Kbar_squared_subset_m3_Tnil':True,'QQ_dimension_Kbar_squared_upper':2},
 'dimension_proof':'m3Tnilhaslength3 andr*m3=m4length1 byacceptedrreduction/Loewy55321. Thuskerronm3hasdimension2. EntireKbar²liesinthiskernel usingq2ann+JK²moduleclosure.',
 'negative_H0_bridge_only':'H0barhasdimension3 andinjectsinNbar; Kbar²imagehasdimension≤2, soH0cannotbeD-generatedbyK². Thisdoesnotdeterminerank0/1/2, H0length/ellann, nilcubeororiginalindex.',
 'no_new_product_values_or_jets':True,'selective_receipts':receipts}
 owned=S.count(d)+S.count(evidence)+20000;assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'corrected-raw-C1-plus-sign-H27-zero-and-entire-bar-square-dimension-bound.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 r={'id':STEM,'all_passed':True,'status':'CORRECTED_PURE_SOURCE_BAR_SQUARE_R_ANNIHILATOR_DIMENSION_UPPER2','r_Kbar_squared_zero':True,'bar_square_QQ_dimension_upper':2,'original_H0_nil_changed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));signal.alarm(0)
if __name__=='__main__':main()
