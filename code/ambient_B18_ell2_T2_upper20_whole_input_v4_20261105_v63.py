#!/usr/bin/env python3
"""Whole20 SOURCE/coupling/lifetime gate; no polynomial/response values."""
import argparse,json,math,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_INPUT_GATE_V4_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return json.loads((ROOT/p).read_text())
 lifecycle=read('reports/AMBIENT_B18_ELL2_T2_UPPER20_CODEC_LIFECYCLE_20261105_v63/report.json');assert lifecycle['all_passed'] and lifecycle['synthetic_nonmath_lifecycle_fixture_ONLY']
 released=read('reports/AMBIENT_B18_ELL2_T2_UPPER20_CODEC_LIFECYCLE_20261105_v63/synthetic-codec-frame-and-context-release.json');assert released['active_exception_traceback_none'] and released['codec_local_objects_alive_after_release']==[False,False]
 withdrawn=read('reports/AMBIENT_B18_ELL2_T2_UPPER20_VALUES_V2_BEFORE_ARITHMETIC_WITHDRAWAL_20261105_v63/report.json');assert withdrawn['all20responses_UNRUN']
 vpath='code/ambient_B18_ell2_T2_upper20_values_v4_20261105_v63.py';pins[vpath]=H.sha(ROOT/vpath);raw+=2*(ROOT/vpath).stat().st_size
 op='reports/AMBIENT_B18_ELL2_T2_UPPER20_PRUNED_SOURCE_20261105_v63/';orr=read(op+'report.json');assert orr['all_passed'];o=read(op+'pruned20-exceptional-OLD-FIRST-outer-source-and605-coverage.json')
 np='reports/AMBIENT_B18_NATIVE20_PRUNED_SOURCE_20261105_v63/';nr=read(np+'report.json');assert nr['all_passed'];n=read(np+'pruned-native20-formal16-toactual15-and-explicit10zero-kernel-source.json')
 mark=read('reports/AMBIENT_B18_ELL2_T2_UPPER605_POSITIVE_WHOLE_20261105_v63_REVIEW/report.json');assert mark['all_passed']
 nm=read('reports/AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW/report.json');assert nm['all_passed']
 single=read('reports/AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63_REVIEW/report.json');assert single['all_passed']
 assert o['exceptional_cases20']==n['exceptional_cases20'];assert len(n['raw_auxiliary_program'])==16 and len(n['actual_kernel_restrictions'])==15 and len(n['explicit_zero_kernel_certificates'])==10
 p=n['native_primitive_known_cache'];assert p['Gamma_two56']==o['Gamma_two'] and p['D_two_by10_by56']==o['D_two_by10_by56'];assert p['N_known_positions']==[0,1,2,3,4,5,6,15,16]
 oldbasis=[tuple(x) for x in o['ordinary_B21_basis']];nb=[tuple(x) for x in p['basis18']];assert len(oldbasis)==21 and len(nb)==18
 for i,var in enumerate([0,1,4,5,7,8]):assert oldbasis[1+var]==nb[1+i]
 assert oldbasis[10:20]==nb[7:17] and oldbasis[20]==nb[17]
 rec=dict(o['corrected20_source_recipes']);coords={r['case619']:r for r in n['all20_constant_native_coordinate_recipes']};assert set(rec)==set(coords)
 assert all(coords[i]['source_subtract563609']==rec[i] for i in rec) and rec[617]==['-20','0']
 # Literal CONSTANT source reconstruction only: no function or response.
 defs={r['kernel105']:r['original_native_source_combination'] for r in n['actual_kernel_restrictions']+n['explicit_zero_kernel_certificates']};c0=dict(o['original_source_allc0_raw21']);cbits=0
 def add(v,k,q):
  nonlocal cbits
  z=v.get(k,F(0))+q;cbits=max(cbits,abs(z.numerator).bit_length(),z.denominator.bit_length());assert cbits<=354
  if z:v[k]=z
  else:v.pop(k,None)
 for i,r in coords.items():
  left={};right={}
  for k,q in r['full_kernel105_coordinate_pairs']:
   assert k in defs
   for j,w in defs[k]:add(left,j,F(q)*F(w))
  for source,q in c0[i]:
   if source>=640:add(right,source-640,F(q))
  for base,q in zip([563,609],rec[i]):
   if F(q):
    for source,w in c0[base]:
     if source>=640:add(right,source-640,-F(q)*F(w))
  assert left==right
  assert sorted(r['full_kernel105_coordinate_pairs'])==sorted(r['kernel105_coordinate_pairs']+r['explicit_proven_zero_kernel_terms'])
 zeroids=[r['kernel105'] for r in n['explicit_zero_kernel_certificates']]
 for r in n['explicit_zero_kernel_certificates']:assert all(int(q[1])==0 for row in r['two_channels_Gamma_D_sum_DU'] for q in row)
 primitive={k:v for k,v in p.items() if k not in ['Gamma_two56','D_two_by10_by56']}
 kernel_fine=[{'kernel105':r['kernel105'],'two_channels_Gamma_D_sum_DU':r['two_channels_Gamma_D_sum_DU']} for r in n['actual_kernel_fine_Gamma_D_sum_DU']+n['explicit_zero_kernel_certificates']]
 native={'primitive':primitive,'raw_program':n['raw_auxiliary_program'],'restrictions':n['actual_kernel_restrictions'],'kernel_fine':kernel_fine,'source_zero_kernel_ids':zeroids,'raw_point_S0_columns48':n['raw_point_S0_columns48'],'case_coordinates':[[i,coords[i]['full_kernel105_coordinate_pairs']] for i in o['exceptional_cases20']],'actual_kernel_and_formal_raw_scope':n['formal_auxiliary_rule'],'native_case617_ALREADY_corrected_not_added_twice':True}
 # Exact Boolean support counts for the executor's NEW-state ledger.
 bidx={x:i for i,x in enumerate(oldbasis)};table=[[bidx.get(tuple(x+y for x,y in zip(a,b))) for b in oldbasis] for a in oldbasis]
 def prod(a,b):return {table[i][j] for i in a for j in b if table[i][j] is not None}
 leaves={i:{k for k,q in enumerate(v) if F(q)} for i,mode,v in o['boundary_cache']};C=dict(leaves);T={}
 for i,r in o['literal_C0_program']:
  if r[0]=='c':v={0} if r[1] else set()
  elif r[0]=='n':v=C[r[1]]
  elif r[0]=='a':v=C[r[1]]|C[r[2]]
  else:assert r[0]=='m';v=prod(C[r[1]],C[r[2]])
  C[i]=v
 for i,r in o['typed_program']:
  if r[0]=='actual_C0_scalar_function':v=leaves[r[1]]
  elif r[0]=='sum':v=set().union(*(T[k] for k,a,b in r[1] if a)) if r[1] else set()
  else:
   assert r[0]=='product';v={0} if r[1] else set()
   for k in r[3]:v=prod(v,T[k])
  T[i]=v
 L0={j:{i for i,q in row if F(q)} for j,row in o['native_L0_selected_columns']};Le={j:{i for i,q in row if F(q)} for j,row in o['native_Le_selected_columns']};first={}
 for rr in o['FIRST_column_recipes']:
  r=rr['first_column_recipe']
  if r['kind']=='literal56':v={i for i,q in r['column'] if F(q)}
  else:
   assert r['kind']=='native_ell_SUM';v=set()
   for name,mat in [('Se',L0),('S0',Le)]:
    for j,q in r[name]:
     if F(q):v|=mat[j]
  first[rr['pivot27'],rr['variable9']]=v
 factor_count=0
 for i,a1 in o['original_full_a1_raw21']:
  for slot,q in a1:
   var=slot%9
   for a,b in o['pair_order10']:factor_count+=(var==a)+(var==b and a!=b)
 lc=sum(map(len,leaves.values()));cc=sum(len(C[i]) for i,r in o['literal_C0_program']);tc=sum(map(len,T.values()));fc=sum(map(len,first.values()));lentries=sum(len(v) for j,v in o['native_L0_selected_columns']+o['native_Le_selected_columns'])
 oldcols=len(o['OLD_Gamma_columns']);assert oldcols==len(o['OLD_D_columns'])
 # All shapes refer to the actual stdlib executor, including repeated
 # parser forms, source indices, current vectors and final certificates.
 ledger={'parsed_nativeN8x12x9x4':3456,'parsed_G576_nu36':612,'parsed_GammaD_sparse_indexvalue':2*(sum(map(len,o['Gamma_two']))+sum(sum(F(q)!=0 for q in row) for ch in o['D_two_by10_by56'] for row in ch)),
 'OLD_known_leaf_and_phasepoly_shapes':4*lc+2*len(leaves)+2*max(cc,tc)+max(len(o['literal_C0_program']),len(o['typed_program'])),
 'parsed_L0Le':2*lentries+len(L0)+len(Le),'actual_FIRST62_sparse':2*fc+len(first),'FORMAL_OLDtables':6*oldcols,
 'raw21partials_and_canonical_FIRSTfactors':21*9+2*factor_count+21,'native_aux_point_kernel_responses':16*5+200+25*3,
 'current_nativeScNX_sparse':2*2*72+7*72*2+9,'all20_results_classes_and_rankcertificate':2000,'monomial_index_tables_and_route_scratch':2000,'source_gate_peer_custody_metadata':2000}
 derived=sum(ledger.values());known_NEW=2*max(cc,tc)+max(len(o['literal_C0_program']),len(o['typed_program']))+2*fc+len(first)+6*oldcols+21*9+2*factor_count+21+16*5+200+25*3+2*2*72+7*72*2+9+2000+500
 # Two extra codec copies beyond the already counted NEW state, plus
 # pending operands/cursor, with a proven ROUND-UP reserve.
 guard_extra=2*known_NEW+1000
 E={(name,t):tuple(map(int,q)) for name,t,q in o['positive_global_envelopes']};U0=E['upper',0][1];U1=E['upper',1][1];detU=2*U0*U1;checkU=2*max(U0,U1)*detU
 assert detU.bit_length()<=334 and checkU.bit_length()<=502
 data={'OLD':o,'NATIVE':native,'positive_global_envelopes':o['positive_global_envelopes'],'committed_prior':354,'temporary_prior':1063,'hard_guards':[384,1153],'metadata_reserve':2000,'all20_values_UNRUN':True,'fresh_CODE_custody_independent_SOURCEpeer_and_ROOTGO_required':True,'native_full128_constant_reconstructions_all20_checked':True,'source_coefficient_reconstruction_maxbits':cbits,'source617_plus20_raw563_OLD_FIRST_outer_and_corrected_native_once':True,'physical_coordinate_and_monic_branch_coupling_literal':True,'source_guard_known_input_by_hash_only_NEWstates_encoded':True,'rank_integer_minor_bound':str(detU),'rank_TEMP_comparison_bound':str(checkU),'derived_workspace_ledger':ledger,'NEWstate_numeric_upper_for_two_extra_codec_copies':known_NEW,'guard_extra_codec_reserve':guard_extra}
 success_codec=16*5+25*3+21*9+(2*factor_count+21)+20*21+count(o['all605_coverage'])+10+50
 assert success_codec<4000
 rank_codec=40+20+2+1+2+2+1+40+40+120+30;assert rank_codec<1000
 codec_extra=max(guard_extra,2*success_codec+1000,2*rank_codec+1000,2000)
 assert codec_extra==guard_extra
 data.update({'shared_codec_reserve_rule':'SUCCESS evidence/rank/report codecs are sequential and mutually exclusive withguard. ORIGINAL_NEWstate is alreadyinderived; extra2encodings max(SUCCESS,GUARD), neverSUMplusduplicateOUTPUT4000.','largest_success_codec_numeric_upper':success_codec,'rank_codec_numeric_upper':rank_codec,'shared_codec_extra_numeric_upper':codec_extra,'preserved_before_value_forecast_failure_code_sha256':'5c85b0dd81ac24fa98bf687ac381093f2b3303078799e637a0a1fd3345028fc1','preserved_failure_receipt_sha256':'8400e1321a3e5d148cc80b5d26215bf7c3b805b52de30f4f9fa10532e22a26c9'})
 data.update({'VALUE_V3_traceback_frames_context_cause_and_unreachable_cycles_released_BEFORE_GUARD':True,'VALUE_V3_worker_sha256':pins[vpath],'bounded_nonmath_codec_lifecycle_PASS':True,'computation_success_window_seconds':55,'exclusive_guard_window_seconds':5,'whole_seconds_limit':60,'oldV2_value_GO_premise_WITHDRAWN_before_values':True,'VALUE_V4_env_Cramer_kernel_and_TEMPdebug_prefix_complete_within_work500':True})
 base=count(data);cap=2*(base+1000)+derived+codec_extra;rawcap=2_000_000+cap*(2*math.ceil(354*math.log10(2))+50)
 assert cap<60000 and cap<150000 and rawcap<32*1024*1024
 data.update({'input_numeric_before_finalmetadata':base,'metadata_reserve_per_input_copy':1000,'whole_numeric_cap':cap,'whole_RAW_cap':rawcap,'RAW_codec_rule':'Two input copies plus NEW codec only; atmostwhole_numeric_cap records at2*ceil(354log10(2))+50bytes,2MBsource/rawreserve. SingleTEMPoffender<=1063bits chargedwithinmetadata.'})
 ep=out/'whole20-coupled-pruned-cache-input-and-literal-owner-height.json';H.dump(ep,data)
 assert count(data)-base<1000
 owned=count([o,n,orr,nr,mark,nm,single,data,defs,coords,C,T,first,leaves,L0,Le,ledger,pins])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_coupling_and_lifetime_only_no_actual20responses':True,'native_true_restrictions_not_individual_raw_aux':True,'584sourcezero_plus554actualzero_plus20_complete605_coverage':True,'20literal_full128_constant_native_reconstructions_verified':True,'fine_global_mapped_Q_prior354_temporary1063':True,'derived_whole_numeric_cap':cap,'derived_whole_RAW_cap':rawcap,'source_owned_numeric_upper':owned,'SOURCE_RAW_read_upper':raw,'FIRST_factor_count_per_one_channel_raw21':factor_count,'Boolean_cells_C0_typed_FIRST':[cc,tc,fc],'value_input_path':str(ep.relative_to(ROOT)),'value_input_sha256':H.sha(ep),'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','derived_whole_numeric_cap','derived_whole_RAW_cap','FIRST_factor_count_per_one_channel_raw21','Boolean_cells_C0_typed_FIRST','source_owned_numeric_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
