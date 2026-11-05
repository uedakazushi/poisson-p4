#!/usr/bin/env python3
"""ONE549/J18_with_r2ell native SOURCE only, no response products."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_ONE549_FRAME_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple,set)):return sum(map(count,x))
 return int(isinstance(x,(int,float)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def fields(p,names):
 s=p.read_text();d=json.JSONDecoder();i=s.index('{')+1;z={}
 while True:
  while s[i].isspace() or s[i]==',':i+=1
  if s[i]=='}':break
  k,i=d.raw_decode(s,i)
  while s[i].isspace() or s[i]==':':i+=1
  v,i=d.raw_decode(s,i)
  if k in names:z[k]=v
  del v
 assert set(z)==set(names),(p,set(names)-set(z));return z
def main(out):
 t0=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 origp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';orig=read(origp,['original619_c0_source_combinations','kernel617_source_recipes','D_two_by10_by56','pair_order10'])
 a1p='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json';a=read(a1p,['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'][549]
 assert a['original768_combination']==orig['original619_c0_source_combinations'][549]==[[298,'1/4'],[649,'1']]
 recipe=next(r for r in orig['kernel617_source_recipes'] if r['free_case619']==549);assert recipe['subtract_actual_basis563609']==['0','0']
 aa=a['actual_a1_27_times9'];assert len(aa)==32
 plane=[[i,q] for i,q in aa if i%9 in [7,8]];assert plane==[[61,'-79/24'],[62,'-11/30']]
 ccpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json';cc=read(ccpath,['Gamma_two','pivot27_source_indices','P27_columns56','Pi3_columns56','native_S0'])
 assert cc['pivot27_source_indices'][6]==648 and cc['native_S0'][9]==[] and cc['native_S0'][8]==[[15,'-1/2']]
 opath='reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json';op=read(opath,['native_S12_by32','basis','structure16','chain_identity1'])
 ss=op['native_S12_by32'];assert [v for v in ss['koszul28_in_first8_rows'] if v[1]==2]==[[0,2,3,1],[3,2,0,-1]]
 mup='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';mu=read(mup,['native_S0_columns']);assert mu['native_S0_columns'][9]==cc['native_S0'][9] and mu['native_S0_columns'][8]==cc['native_S0'][8]
 c68p='reports/AMBIENT_B18_NATIVE_CANONICAL68_VALUES_20261105_v63/actual-G8-C68-and-monic-second-values.json';c68=read(c68p)
 c28p='reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_VALUES_20261105_v63/actual32-G8-r2ell-and-two-monic-cached-coefficients.json';c28=read(c28p)
 for p in ['reports/AMBIENT_B18_NATIVE_CANONICAL68_DIRECT_REVIEW_20261105_v63/report.json','reports/AMBIENT_B18_NATIVE_C28_DIRECT_REVIEW_20261105_v63/report.json']:assert read(p)['all_passed']
 cinp='reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/canonical33-to-B17-C68-G8-source-input-and-response-budget.json';cin=read(cinp,['physical_R0_variables','basis_B17','actual_nu_A_B33_cached_UNEVALUATED','canonical30_branch','source_G8_order','free4_basis_mask_order'])
 basis=[tuple(v) for v in cin['basis_B17']]+[(0,0,0,0,0,0,0,2,1)];bi={m:i for i,m in enumerate(basis)};assert len(bi)==18
 selector=[0,7,6,3,2,8,9,22,23,20,21,14,15,12,13,24,25,30]
 tablep='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json';table=read(tablep,['basis33','table_columns_for_each_left_basis'])
 rb=[]
 for z,r,e in table['basis33']:
  m=[0]*9
  if z:m[z-1]=1
  m[7]=r;m[8]=e;rb.append(tuple(m))
 images=[bi.get(m) for m in rb];assert selector==[images.index(i) for i in range(18)]
 nfp='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';nf=read(nfp,['monomials660_degree_le4','all660_actual_raw33_normal_forms'])
 sourcepre=count([orig,a,aa,cc,op,mu,c68,c28,cin,basis,table,nf,pins])+40000;assert sourcepre<150000 and raw+2_000_000<32*1024*1024
 dump('phase.json',{'phase':'literal_current549_S0_column_and_preoperation_source_owner_PASS','sourcepre':sourcepre,'no_native_response_values':True})
 for m,v in zip(nf['monomials660_degree_le4'],nf['all660_actual_raw33_normal_forms'],strict=True):
  projected={images[k]:q for k,q in v if images[k] is not None and q!='0'};assert projected==({bi[tuple(m)]:'1'} if tuple(m) in bi else {})
 for i in range(33):
  for j in range(33):
   p={images[k]:q for k,q in table['table_columns_for_each_left_basis'][i][j] if images[k] is not None and q!='0'};m=tuple(a+b for a,b in zip(rb[i],rb[j]));idx=bi.get(m) if images[i] is not None and images[j] is not None else None;assert p==({idx:'1'} if idx is not None else {})
 G72=[v+q for v,q in zip(c68['full_G8_C68'],c28['actual_G8_r2ell_eight_by4'],strict=True)]
 nu18=[[q[k] for k in selector] for q in cin['actual_nu_A_B33_cached_UNEVALUATED']]
 assert c68['G8_order']==c28['source_G8_order']==cin['source_G8_order'] and [v[17] for v in nu18]==c28['actual_nu_A_B_r2ell']==['0','0']
 assert G72[0][:4]==['0','0','0','1/2'] and G72[3][:4]==['0','0','0','0'] and all(v[0]=='0' for v in nu18)
 # The factor tauA on inputmask1 kills G0(0)=tauAtauB/2 in C4;
 # compare the entire original48-component Res4 column, not a G-entry test.
 fullS0_zero=mu['native_S0_columns'][9];assert fullS0_zero==[]
 # All monic pair products vanish by support in this specific J18 ring. This
 # fact is source-only, does not assign N mixed/higher coefficients values.
 ns=[{i for i,q in enumerate(v) if q!='0'} for v in nu18]
 for x in ns:
  for y in ns:
   assert all(tuple(a+b for a,b in zip(basis[i],basis[j])) not in bi for i in x for j in y)
 nrp='reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json';nr=read(nrp,['N_C_B44'])['N_C_B44']
 nep='reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json';ne=read(nep,['N'])['N']
 nfp0='reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json';nfirst=read(nfp0,['N'])['N']
 known_positions=[0,1,2,3,4,5,6,15,16];Nknown=[]
 for i in range(12):
  row=[]
  for j in range(12):
   b=nr[i][j];e=ne[i][j];f=nfirst[i][j]
   for rawidx,dualidx in [(0,0),(3,1),(4,2),(5,5),(6,6),(1,8)]:assert b[4*rawidx:4*rawidx+4]==f[4*dualidx:4*dualidx+4]
   assert b[:8]==e[:8]
   row.append([b[:4],b[12:16],b[16:20],b[20:24],b[24:28],b[4:8],e[8:12],b[8:12],e[12:16]])
  Nknown.append(row)
 for p in ['reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63_REVIEW/report.json']:assert read(p)['all_passed']
 external='reports/AMBIENT_B18_ELL2_T2_UPPER549_OLD_CARRY_SOURCE_20261105_v63/report.json';er=read(external);assert er['all_passed']
 ev={'selected_case619':549,'literal_full_original768_source':a['original768_combination'],'native_original768_source649_coefficient':'1','native128_source9_Koszulsource2_inputmask1':True,
  'full_a1_27times9_32_entries_UNEVALUATED':aa,'same617_lower_basis_recipe':recipe,
  'full_intrinsic_Res4_S0_column9_exactempty':fullS0_zero,'comparison_pivot648_mask0_S0_column8_nonzero':cc['native_S0'][8],
  'canonical_G8_point_C4_rows0_and3':[G72[0][:4],G72[3][:4]],'full_original_source_S12_by32':ss,
  'Gamma_two_current56':cc['Gamma_two'],'D_two_by10_by56':orig['D_two_by10_by56'],'D_ordinary_pair_order10':orig['pair_order10'],
  'P27_columns56_current':cc['P27_columns56'],'Pi3_columns56_current':cc['Pi3_columns56'],
  'ring_name':'J18_with_r2ell; distinct from older A18/C72(m3)','ordinary_basis_J18':basis,'raw33_selector18':selector,
  'all660_nonhomogeneous_NF_and1089_products_checked':True,
  'fiber_product_source':'J18_with_r2ell=B17 times_P5 P6, P6=P7/(r^3); dimensions17+6-5=18; same raw33 R0 quotient.',
  'trace_basis_order':cin['free4_basis_mask_order'],'canonical_branch':cin['canonical30_branch'],
  'already_accepted_G8_J18_C72_images_COPIED_NOT_REEVALUATED':G72,'already_accepted_nu_two_J18_images':nu18,
  'G8_order':c68['G8_order'],'all_monic_pair_products_zero_by_literal_support_not_cancellation':True,
  'N_original_row_column_order':'N[source12][physical12], transpose12 only; coefficient Res4 uses varyingnu and never transposes mu4.',
  'N_known_J18_positions':known_positions,'N_known12_by12_by9_by4_cached_operands_ONLY':Nknown,
  'N_unknown_J18_positions_remain_FORMAL':list(range(7,15))+[17],
  'unknown_second_native_D_cancellation':'D degree2: all unknown L2_pair multiply ENTIRE Res4(S0)*e9=0; no M0-only inference.',
  'unknown_third_native_Gamma_cancellation':'Gamma pure r²ell cubic: unknown L3_r²ell multiplies ENTIRE Res4(S0)*e9=0; all unknown nonplane L2 have z-factor and cannot enter this pureplane coefficient.',
  'known_native_cubic_factored_template':'Gamma[L0 S3_r²ell +Lr S2_rell+Lell S2_rr+Lrr Sell+Lrell Sr]; all varyingmonic terms inside fullRes4 L/S retained.',
  'known_native_D_second_template':'for all10 ordinary pairs D_pair[L0 S2_pair + sum firstLi*firstSj]; L2*S0 col9=0 ONLYfor this source column.',
  'outer_native_pivot648_template':'-a1[6,ell]*GammaXp648_rr-a1[6,r]*GammaXp648_rell; BOTH fullL2S0 target references retained. No numerical target values read here.',
  'only_pureplane_a1_positions':plane,'all_other_a1_z_factors_do_not_enter_Gamma_pureplane':True,
  'external_OLD298_quarter_and_full32_knownFIRSTcarry_lane':{'source_report':external,'source_sha256':pins[external],'independent_peer_not_assumed_by_this_packet':True},
  'no_scalar_N2_N3_zero_fill_or_admissible_gauge_assumption':True,
  'no_native_S2_S3_Gamma_D_response_products_computed':True,'future_positive_height_whole_value_gate_and_worker_are_SEPARATE':True,'future_value_authorized':False}
 owner=count([orig,a,aa,cc,op,mu,c68,c28,cin,basis,table,nf,nr,ne,nfirst,Nknown,G72,nu18,ev,pins])+20000
 assert owner<150000 and raw+2_000_000<32*1024*1024
 ep='ONE549-J18-with-r2ell-native-full-S0-column-and-formal-N2N3-frontier.json';dump(ep,ev)
 rep={'id':STEM,'all_passed':True,'source_only_no_native_response_values':True,'J18_with_r2ell_and_C72_dimension':[18,72],
  'full_source649_intrinsic_S0_column_zero_pivot648_distinguished':True,'all_unknown_N2_N3_dependence_factorized_not_zero_filled':True,
  'all_accepted_canonical_G8_monic_and_N_known_jet_branch_bindings_checked':True,
  'whole_SOURCE_owned_upper':owner,'whole_SOURCE_RAW_upper':raw+2_000_000,'future_whole_positive_value_gate_still_required':True,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-t0}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'RAW':raw+2_000_000,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
