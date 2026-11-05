#!/usr/bin/env python3
"""Typed source identity only: combine explicit A2/M1 into a known plane row."""
import json
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_C_COMPRESSION_SOURCE_20261105_v63'
def main():
 pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 def bind(p,field,shape):
  v=read(p,[field])[field];return {'path':p,'sha256':pins[p]['sha256'],'field':field,'shape':shape},v
 for p in ['reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER_FREE3_GAUGE_SOURCE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63_REVIEW/report.json']:assert read(p)['all_passed']
 invariants,inv=bind('reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63/actual-two29-pairings-and-recovered30-kernel-invariants.json','actual_all_three_point_kernel_direction_invariants',[3,2]);assert set(inv)=={'23','29','30'} and all(len(v)==2 for v in inv.values())
 chart,mini=bind('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json','fixed_bottom30_inverse_columns',[30,30]);assert len(mini)==30
 mp='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json';mi=read(mp,['free30_V_indices','original27_X_source_columns']);ids=mi['original27_X_source_columns'];assert len(ids)==27
 Gbind,G=bind('reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63/617-direct-old-highestthird-coherent-source-positive-input.json','Gamma_two',[2,56]);assert len(G)==2
 gamma_rows={i for row in G for i,x in row};assert max(gamma_rows)==55
 first0,A=bind('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9',[640,56,9]);assert len(A)==640
 firstN,N=bind('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9',[128,48,9]);assert len(N)==128
 firste,E=bind('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json','actual_first640_ell_columns',[640,56]);assert len(E)==640
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';nat=read(mug,['L0_Res4_Ntranspose_columns','native_S0_columns','native_S1ell_columns']);assert len(nat['L0_Res4_Ntranspose_columns'])==48 and len(nat['native_S0_columns'])==len(nat['native_S1ell_columns'])==128
 leb,Le=bind('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json','actual_L1ell_columns',[48,48]);assert len(Le)==48
 pivot_recipe=[]
 for p,j in enumerate(ids):pivot_recipe.append({'pivot27_index':p,'original_source768':j,'r_first':{'binding':first0 if j<640 else firstN,'column_index':j if j<640 else j-640,'dual_slot8':True,'source_row_limit':56 if j<640 else 48,'bottom8_native_zero':j>=640},'ell_first':{'binding':firste,'column_index':j} if j<640 else {'operation':'SUM','terms':['L0*S1ell_column','L1ell*S0_column'],'local128_column':j-640,'all_individual_LeS0_carries_retained':True}})
 node=[{'kind':'point_projections','P27':'first27rows(chart_inverse) on SAME free30 V indices','Pi3':'last3rows(chart_inverse), NOT sourceindices8/11/12','binding':chart,'free30_indices':mi['free30_V_indices']},
 {'kind':'known_plane_Xp1_callbacks27','callbacks':pivot_recipe},
 {'kind':'known_plane_M1H0_in_existential_coherent_frame','r_rows':['alpha23,alpha29,alpha30','0,0,0'],'ell_rows':['beta23,beta29,beta30','0,0,0'],'invariant_binding':invariants,'H0_original_source_axes':[23,29,30],'canonical_other_M1_values_not_assigned':True},
 {'kind':'shared_known_plane_linear_C','formula':'C_t= (lambda_t M1 H0)_plane Pi3 - Gamma_t (Xp1)_plane P27','Gamma_binding':Gbind,'shape':[2,2,56],'all_numeric_C_coefficients_uncomputed':True},
 {'kind':'shared99_quadratic_K','formula':'K_tj=Gamma_t Xstar_j2+C_t Xstar_j1','rr':'Gamma Xstar_rr+C_r Xstar_r','rell':'Gamma Xstar_rell+C_r Xstar_ell+C_ell Xstar_r','actual_K_values_uncomputed':True},
 {'kind':'compressed617_upper_source','formula':'upper_t=[Gamma_t Xc3-Gamma_t Xp2*a1+C_t*(Xc2-Xp1*a1)-sum_j K_tj*b1_j]_(r^2 ell)','source_cases617':'Useexactconstantlowerkernel617recipesBEFOREcoefficientqueries','higher_original_functions_not_replaced':True},
 {'kind':'b1_source_frontier','formula':'b1=B99_point*(Ybarc1-G3linear_J*q,0); q=F25 PiY0 Ybarc1','full_domain_slots':918,'unknown_Yp0_retained_inside_Ybarc1':True,'new_rowspan_required':'No inference that ker(G2linear)=im(G3linear); any adjoint/lift independence must be verified on the full918 domain.'}]
 scope={'algebraic_proof':'Substitute A2=P27P2 and B2=Pi3P2 into the peered upper identity, where P2=Xc2-Xp1a1-Xstar1b1. Collect the common row C=lambda M1H0 Pi3-Gamma Xp1 P27. Expansion gives C(Xc2-Xp1a1)-(Gamma Xstar2+C Xstar1)b1, with the highest/pivotsecond terms unchanged. This is a distributive source identity, not an assumed vanishing.',
 'plane_selection':'Only r/ell first components contribute to ordinary cubic output r^2ell BEFORE any nonhomogeneous conductor normal form. Unknown z-components of M1 are not assigned0; they cannot contribute to this pure output. No factorial is introduced.',
 'known_and_unknown':'C can be constructed from cached six invariant values, pointchart, Gamma and all27 actual plane firsts. Its numerical coefficient construction still needs a separate positive/count/source gate. K uses source seconds and b1 uses actual Ybar1/Yp0; those responses remain unknown. Only the explicit A2/M1H0 expression is compressed, not every possible hidden comparison choice.',
 'individual_pivots':'Native pivot ell first must be L0*S1ell+L1ell*S0 with the full SUM. Weighted source S0=0 for617 does not hold for individual pivots. Original fullR0 regular quotient functions remain the same before ell2 basechange.'}
 owned=count([mi,inv,G,A,N,E,nat,Le,node])+8000;assert owned<60000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'typed-shared-plane-C-and99-K-compressed617-source-identity.json';ep.write_text(json.dumps({'typed_nodes':node,'source_proof_and_counterconditions':scope,'known_C_plane_coefficient_shape':[2,2,56],'source_owned_count_bound':owned,'numeric_C_K_or_actual_upper_values_not_computed':True,'future_required':['positiveknownC/coefficientconstruction gate','shared99K source/cache/sensitivity gate','full918adjoint or exactlift-independence gate','wholecombinedvalueheight/ownership/GO']},separators=(',',':'))+'\n');r={'id':STEM,'all_passed':True,'typed_source_identity_only':True,'explicit_A2_M1H0_combined_without_zero_filling':True,'unknown_b1_Yp0_K_and_higher_responses_retained':True,'no_numeric_C_K_or_upper_values':True,'source_owned_count_bound':owned,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(r,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned}))
if __name__=='__main__':main()
