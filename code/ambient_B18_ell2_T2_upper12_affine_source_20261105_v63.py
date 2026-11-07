#!/usr/bin/env python3
"""Remaining12 affine24-symbol/factorgroup SOURCE, no response/rank query."""
import argparse,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER12_AFFINE_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return H.selected(p,k)
 def ck(q):
  nonlocal bits
  bits=max(bits,abs(q.numerator).bit_length(),q.denominator.bit_length());assert bits<=354;return q
 def add(v,k,q):
  z=ck(v.get(k,F(0))+q)
  if z:v[k]=z
  else:v.pop(k,None)
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';s=read(sp)
 fp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json';f=read(fp);allformal=dict(f['all617_combined_OLD2_independent_symbol_coefficients'])
 frontier=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63/literal12-opaque-function-and24-plane-second-frontier-with-existing-rr-labels.json')
 peer=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63_REVIEW_V2/report.json');assert peer['all_passed']
 shared=read('reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63/complete605-shared-source-linearity-and-global24-opaque-frontier.json');cases=shared['opaque_cases12'];assert cases==[579,580,593,594,597,598,605,606,611,613,615,616]
 Gamma=[{i:F(q) for i,q in row} for row in s['outer_four_pivot_source']['Gamma_two']]
 groups=[[330164,330184,566400,330168],[330289,330309,566401,330293],[330524,330544,566402,330528],[330594,330614,566403,330598]]
 opaque={fn for g in groups for fn in g[:3]};assert len(opaque)==12
 rawplane={i:dict(v) for i,v in s['raw619_old4_plane_a1_nonempty28coefficients']};recipes={r['free_case619']:r['subtract_actual_basis563609'] for r in s['kernel617_source_recipes']}
 plane=[];matrix=[];latent=[];known_rr_offset=[]
 assert count([s,f,frontier,peer,shared,rawplane,recipes])+25000<100000
 for case in cases:
  a={i:F(q) for i,q in rawplane.get(case,{}).items()}
  for base,q in zip([563,609],recipes[case]):
   if F(q):
    for slot,w in rawplane.get(base,{}).items():add(a,slot,ck(-F(q)*F(w)))
  assert all(slot//9<4 and slot%9 in [7,8] for slot in a);plane.append([case,[[i,str(q)] for i,q in sorted(a.items())]])
  for t in range(2):
   expected={};A={};B={};rr={}
   for k,(gold,hrest,gnew,rawr) in enumerate(groups):
    cr=a.get(9*k+7,F(0));ce=a.get(9*k+8,F(0));g4=Gamma[t].get(4,F(0));g28=Gamma[t].get(28,F(0));g55=Gamma[t].get(55,F(0));gg=ck(g28+g55)
    for pair,co in [(8,ce),(9,cr)]:
     for fn,g in [(hrest,g4),(gold,gg),(gnew,-g55)]:add(expected,(fn,pair),ck(-co*g))
    # Twelve shared unknown FUNCTIONALS after the8 cached rr values:
    # Hrr[k] and Ere[t,k]; this is a factorization, not a rank claim.
    add(B,k,ck(-ce*g4));add(B,4+4*t+k,-cr)
    add(rr,(gold,8),ck(-ce*gg));add(rr,(gnew,8),ck(ce*g55))
   for pair,fn,q in allformal[case][t]:
    if fn in opaque:add(A,(fn,pair),F(q))
   assert A==expected
   matrix.append([case,t,[[fn,pair,str(q)] for (fn,pair),q in sorted(A.items())]])
   latent.append([case,t,[[i,str(q)] for i,q in sorted(B.items())]])
   known_rr_offset.append([case,t,[[fn,pair,str(q)] for (fn,pair),q in sorted(rr.items())]])
 j1=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json');J={r['free_case619']:r for r in j1['all617_formal_directOLD_J0_Jr_Jell']}
 thirdcases={}
 for case in cases:
  supp=sorted({fn for row in J[case]['two_formal_J0_Jr_Jell_symbol_rows'] for fn,v in row if fn in [299650,505938] and v[0]!='0'})
  if supp:thirdcases[case]=supp
 assert set(thirdcases)=={579,613,616}
 kp='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json';aa=field(kp,'all_zero_bottom_and_Y0_source_combinations');rawcases=set(cases)
 for i in cases:
  for b,q in zip([563,609],recipes[i]):
   if F(q):rawcases.add(b)
 c0=[[i,aa[i]['original768_combination']] for i in sorted(rawcases)];fulla1=[[i,aa[i]['actual_a1_27_times9']] for i in sorted(rawcases)];del aa
 nm=read('reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/all617-native-constant-sources-factor-through-fixed105-kernel.json');nm={r['free_case619']:r for r in nm['all617_source_maps'] if r['free_case619'] in cases};assert set(nm)==set(cases)
 native_basis=sorted({k for r in nm.values() for k,q in r['kernel105_coordinates']})
 # Presence inventory only. Do not read/copy any existing scalar response.
 rrbind=frontier['accepted_existing_rr_only_cache_binding'];rr8=field(rrbind['data'],rrbind['field']);assert {r[0] for r in rr8}==set(rrbind['functionIDs']);del rr8
 slots=frontier['opaque24_scalar_symbols'];known8=[[i,8] for i in rrbind['functionIDs']];remaining16=frontier['remaining_individual_scalar_frontier_after_rr_cache_presence'];assert len(known8)==8 and len(remaining16)==16
 firstpoint=field('reports/AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63/pruned605-OLD757-C0532-B21-FIRST162-and-native92-source-input.json','individual_FIRST_column_recipes162');have={(r['pivot27'],r['variable9']) for r in firstpoint};need=set()
 for i,a1 in fulla1:
  for slot,q in a1:
   p,var=divmod(slot,9)
   for x,y in s['pair_order10']:
    if var==x:need.add((p,y))
    if var==y and x!=y:need.add((p,x))
 assert need<=have
 gbs=read('reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json');assert gbs['all_passed']
 data={'SOURCE_only_no_rank_adjoint_actual_pairing_or_scalarfunction_values':True,'cases12':cases,'corrected_OLD4_plane_a1_from_full563609_carry':plane,'source_groups_gold_restoredh_newg_rawr':groups,'Gamma_two4_28_55_exact':[[[row,str(Gamma[t].get(row,F(0)))] for row in [4,28,55]] for t in range(2)],'all24_unknown_symbol_columns':slots,'all24_equations_exact_opaque_coefficient_matrix_rows':matrix,'all24_equations_match_existing_independent_full617_formal_certificate':True,'known8_rr_value_REFERENCES_ONLY':rrbind,'known8_rr_symbols':known8,'remaining16_individual_second_symbols_not_unavoidable_newqueries':remaining16,'twelve_shared_unknown_functionals':[{'index':k,'kind':'Hrr','restoredh_function':g[1],'ordinary_pair8':True} for k,g in enumerate(groups)]+[{'index':4+4*t+k,'kind':'Ere','channel':t,'ordinary_pair9_terms':[[g[1],str(Gamma[t].get(4,F(0)))],[g[0],str(ck(Gamma[t].get(28,F(0))+Gamma[t].get(55,F(0))))],[g[2],str(-Gamma[t].get(55,F(0)))]]} for t in range(2) for k,g in enumerate(groups)],'factor24_by12_coefficient_rows':latent,'known_rr_offset_coefficient_rows_NOT_paired':known_rr_offset,'structural_unknown_response_dimension_upper12_NOT_matrixrank':True,'restoredh_second_exact_source_identity':'restoredh=rawr-h1*gold; h1=327554, both origins0. rr: rawr_rr-h1_r*gold_r. rell: rawr_rell-h1_r*gold_ell-h1_ell*gold_r. AllfirstsarecachedFULL9; doNOTassignrawr_second0.','literal_original_regular_division_ancestors':frontier['literal_original_C0_arithmetic_and_regular_quotient_nodes'],'literal_regular_quotient_frontier':frontier['literal_quotient_frontier'],'division_rule':'Performformalregular/orderedprefix division BEFOREanyfiniteNF. gdenominators299650/505938 order3; rdenominator27948 haspoint/r/ell0 andsourcez1 unitfirst. No quotient value followsfromfinite33projectionalone.','raw_required_cases_with_full_c0':c0,'raw_required_cases_with_full_a1':fulla1,'source563609_full_correction_for12':[[i,recipes[i]] for i in cases],'known_TRUE_conductor_THIRD_required_cases':thirdcases,'known_gA3_recipe':'505938=prefix0[7790]; onCAplane itsr2ell isliteral7790 raw33selector30','known_gB3_recipe':'literal28539_r2ell minusSUM7(zB_rr*h_ell+zB_rell*h_r); source acceptedCURRENTnegativeboundary. h0=zB0=zBplane1=0. No h2/h3 orfalsegBthird0.','physical_directOLD_thirds_and_combined_seconds_existing164_B21_known_only':True,'native12_constant_kernel_maps_SOURCE':[[i,nm[i]] for i in cases],'native105_basis_indices_needed':native_basis,'native_response_source_operator_reference':'reports/AMBIENT_B18_NATIVE_KERNEL105_FUNCTIONAL_SOURCE_20261105_v63/report.json','native_FRESH_RESPONSEgate_needed_notnewphysicalIFT':'AllcanonicalG8/nu nonzz10+r2ell andN0/Nfirst/Nplane2cached; unknownN2mixed/N3cancelONLYaftertruefullS0kernel. Noactualnative12responsesqueriedhere.','FIRST_needed_pivot_var_pairs':sorted(need),'FIRST_full_actual_cached162_reference':'reports/AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63/report.json','native_fullL2S0_all92_and_opaqueOLD4_source_distinct':True,'AFFINE_known_offset_K12_uncomputed':'ΓdirectOLD3 +DOLD2knownphysical/FIRSTcross +nativekernel −FIRSTcarry −native92outer +known8rrOFFSET +TRUEgA/gBthirds. All24unknownsecondsstayindependentuntilsourcecoherenceproof.','no_more_than12_shared_ancestor_functionals_suffice_but_none_claimednecessary':True,'future_rank_or_adjoint_or_knownK_pairing_requires_freshwhole_SOURCE_height_owner_code_custody_peer_ROOTGO':True}
 ep=out/'remaining12-affine24-second-symbols-factor12-and-cache-ancestor-source.json';H.dump(ep,data)
 owned=count([s,f,frontier,peer,shared,j1,J,rawplane,recipes,Gamma,plane,matrix,latent,known_rr_offset,data,firstpoint,nm,pins])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_exact_affine_and_group_factor_inventory_only':True,'24_symbol_matrix_matches_full617_formal_rows':True,'known8_rr_refs_leave16_individual_unknowns':True,'factor_through12_shared_functionals_NO_rank_claim':True,'gA_gBthird579613616_not_zeroed':True,'native105_FIRST92_known_source_only_no_pairing':True,'constant_coefficient_factorization_actual_maxbits':bits,'source_owned_upper':owned,'SOURCE_RAW_read_upper':raw,'all12_upper_values_UNRUN':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','constant_coefficient_factorization_actual_maxbits','source_owned_upper','SOURCE_RAW_read_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
