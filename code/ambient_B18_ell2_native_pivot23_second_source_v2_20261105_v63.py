#!/usr/bin/env python3
"""Source-only binding of canonical C68 to individual native pivot23.

No signed L2, intrinsic S2, Xp2, Gamma-Xp2 or response value is constructed.
Existing coefficient arrays are copied and positive input-height bounds only
are formed. The factored expressions retain the complete L2*S0 term.
"""
import argparse, hashlib, json, math, resource, signal, time
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_NATIVE_PIVOT23_SECOND_SOURCE_20261105_v63_V2'
PATHS={
 'canonical_report':'reports/AMBIENT_B18_NATIVE_CANONICAL68_VALUES_20261105_v63/report.json',
 'canonical_values':'reports/AMBIENT_B18_NATIVE_CANONICAL68_VALUES_20261105_v63/actual-G8-C68-and-monic-second-values.json',
 'canonical_peer':'reports/AMBIENT_B18_NATIVE_CANONICAL68_DIRECT_REVIEW_20261105_v63/report.json',
 'canonical_source':'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/report.json',
 'Nrr_report':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/report.json',
 'Nrr_peer':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW/report.json',
 'Nrr':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
 'Nre_report':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/report.json',
 'Nre_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63_REVIEW/report.json',
 'Nre':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',
 'operator_report':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/report.json',
 'operator':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'chart_report':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/report.json',
 'chart':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',
 'frame_report':'reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/report.json',
 'frame':'reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json',
 'Gamma_source':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json',
 'origin_report':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/report.json',
 'origin_adapter':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
}
PINS={'canonical_report':'3d24183702989cc457b87adec810fe95e980b24a372450a76af0d90637da5ccf',
 'canonical_peer':'881c090220b39ecd7aac04eead7190630f3840fe34f575c9ee4dadbb8e339297',
 'canonical_source':'48ecdc989f97b4d59a383d87feddaa59b671c86a6391b2a13c9177a8137820cd',
 'Nrr_report':'1f8cdb89b65d2d564977b8339d4d81bd9b49a4147c079c46c8f00273786d54f3',
 'Nre_report':'6b873099a3a407716c8a30c2c74e4c3036880870e2c78c0a06fbaa24df960d9b',
 'frame_report':'7ac51097df48c4416f377fa94380f3abb910f6f10bae77945f96e93ad76b30b4'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(tuple,list)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def fields(path,names):
 """Stream top-level JSON values; retain only selected fields, no producer."""
 text=path.read_text();dec=json.JSONDecoder();i=text.index('{')+1;out={}
 while True:
  while text[i].isspace() or text[i]==',':i+=1
  if text[i]=='}':break
  key,i=dec.raw_decode(text,i)
  while text[i].isspace() or text[i]==':':i+=1
  value,i=dec.raw_decode(text,i)
  if key in names:out[key]=value
  del value
 assert set(out)==set(names),(path,set(names)-set(out))
 return out
def fractions(v):
 if isinstance(v,dict):
  for x in v.values():yield from fractions(x)
 elif isinstance(v,(list,tuple)):
  for x in v:yield from fractions(x)
 else:yield F(v)
def DU(v):
 xs=list(fractions(v));d=1
 for x in xs:d=math.lcm(d,x.denominator)
 u=max([d]+[abs(x.numerator)*(d//x.denominator) for x in xs])
 return {'D':str(d),'U':str(u),'D_bits':d.bit_length(),'U_bits':u.bit_length(),'entries':len(xs)}

def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={};reports={};raw=0
 def read(k,names=None):
  nonlocal raw
  p=ROOT/PATHS[k];pins[PATHS[k]]=sha(p);raw+=p.stat().st_size
  if k in PINS:assert pins[PATHS[k]]==PINS[k],k
  return json.loads(p.read_text()) if names is None else fields(p,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 for k in PATHS:
  if k.endswith('report') or k.endswith('peer') or k=='canonical_source':
   reports[k]=read(k);assert reports[k]['all_passed'],k
 for key,report in [('canonical_values','canonical_report'),('Nrr','Nrr_report'),('Nre','Nre_report'),('operator','operator_report'),('chart','chart_report'),('frame','frame_report')]:
  name=Path(PATHS[key]).name
  assert reports[report]['evidence_sha256'][name]==sha(ROOT/PATHS[key]),(key,report)
 canonical=read('canonical_values')
 nr=read('Nrr',['N_C_B44','TH12_equals_G8_zero4','NF12_equals_dG8_zero4','all_actual_constants_and_firsts_match_frozen_source'])
 ne=read('Nre',['N','actualG8','TH12_G8_zero4','NF12_dG8_zero4','all_N0_rfirst_matches','all_G8_firstell_matches'])
 op=read('operator',['basis','structure16','native_S12_by32','matrix_Res4','chain_identity1','actual_Q4_by48_column_node_ids','all48_actual_Qcolumns_independently_bound_to_mu_F'])
 chart=read('chart',['original27_X_source_columns','free30_V_indices','remaining_free3_originalV_rows'])
 frame=read('frame',['Gamma_two','pivot27_source_indices','pivot27_ell_recipes','P27_columns56','Pi3_columns56','native_pivot_ell_is_SUM_L0Se_LeS0','full_coordinate_not_pureU'])
 gamma=read('Gamma_source',['two_Gamma_on56_actual_X2'])
 origin=read('origin_adapter',['N','G8_C0','all_source_H12_to_G8_identities','d_not_inverted'])
 # Before any coefficient-copy/positive products, charge all held source arrays,
 # their JSON-decoded current duplicates and the future payload/template reserve.
 initial=count([reports,canonical,nr,ne,op,chart,frame,gamma,origin,pins])
 assert initial+40000<150000,(initial,'source owner')
 assert raw+4_000_000<32*1024*1024
 dump('phase.json',{'phase':'all_literal_hashes_and_preoperation_source_owner_PASS',
  'initial_source_numeric':initial,'reserve':40000,'raw_read_plus_output':raw+4_000_000,
  'no_signed_native_pivot_response':True})
 assert canonical['all352_B11_and128_C16_overlap_coefficients_verified']
 assert op['all48_actual_Qcolumns_independently_bound_to_mu_F']
 assert nr['TH12_equals_G8_zero4'] and nr['NF12_equals_dG8_zero4']
 assert ne['TH12_G8_zero4'] and ne['NF12_dG8_zero4']
 assert nr['all_actual_constants_and_firsts_match_frozen_source'] and ne['all_N0_rfirst_matches'] and ne['all_G8_firstell_matches']
 assert origin['all_source_H12_to_G8_identities'] and origin['d_not_inverted']
 piv=[604,613,622,631,640,644,648,652,656,660,664]+list(range(752,768))
 assert chart['original27_X_source_columns']==frame['pivot27_source_indices']==piv
 assert frame['Gamma_two']==gamma['two_Gamma_on56_actual_X2']
 assert frame['native_pivot_ell_is_SUM_L0Se_LeS0'] and frame['full_coordinate_not_pureU']
 native=[j-640 for j in piv[4:]]
 assert frame['pivot27_ell_recipes'][4:]==[['native',j] for j in native]
 assert native==[0,4,8,12,16,20,24]+list(range(112,128))
 assert all(0<=i<56 for row in frame['Gamma_two'] for i,x in row)
 # Native X has zero bottom8; Gamma itself remains the full56 point covector.
 # Dense scalar normalizer coefficients already verified in the original branch.
 # N, not N^{-1}; native map is N^t S. Never transpose the 4x4 mu blocks.
 assert len(nr['N_C_B44'])==len(ne['N'])==12
 N=[]
 for i in range(12):
  row=[]
  for j in range(12):
   a=nr['N_C_B44'][i][j];b=ne['N'][i][j]
   assert len(a)==44 and len(b)==16 and a[:8]==b[:8]
   assert a[:4]==origin['N'][i][j]
   row.append([a[:4],a[4:8],b[8:12],a[8:12],b[12:16]])
  N.append(row)
 gs=canonical['full_G8_C68'];assert all(len(g)==68 for g in gs)
 plane_positions=[0,5,6,15,16]
 G=[[[g[4*b+m] for m in range(4)] for b in plane_positions] for g in gs]
 for i in range(8):
  assert G[i][0]==origin['G8_C0'][i]
  assert sum(G[i][:3],[])==ne['actualG8'][i][:12]
 # All selected monic coefficients are copied, not derived through nonunit
 # division. ν has plane order2, so at most one ν factor contributes to degree2.
 nu=[["0","0","0",canonical['monic_second10'][i][8],canonical['monic_second10'][i][9]] for i in range(2)]
 assert nu==[['0','0','0','-25/192','0'],['0','0','0','725/36','41/9']]
 S=op['native_S12_by32'];bycol={j:[] for j in range(32)}
 for i,j,g,sign in S['koszul28_in_first8_rows']:bycol[j].append(['G8',i,g,sign])
 for i,j,c in S['free4_zero_source_columns']:bycol[j].append(['constant',i,c])
 assert all(len(bycol[j])==2 for j in range(28))
 assert [bycol[j] for j in range(28,32)]==[[['constant',i,1]] for i in range(8,12)]
 recipes=[]
 for k,j in enumerate(native):
  source,mask=divmod(j,4)
  recipes.append({'pivot27_position':k+4,'original768_source':j+640,
   'native128_source':j,'S32_source':source,'input_trace_mask':mask,
   'S_source_entries':bycol[source],
   'target_Gamma_Xp_rr':['Gamma','L0*Srr','Lr*Sr','Lrr*S0'],
   'target_Gamma_Xp_rell':['Gamma','L0*Sre','Lr*Se','Le*Sr','Lre*S0'],
   'individual_S0_is_NOT_assumed_zero':True})
 # A finite universal trace product routing. Coefficient recipes, not values:
 # tau^a tau^c tau^b = tau^(a xor c xor b) nuA^pa nuB^pb.
 trace=[]
 for a in range(4):
  for c in range(4):
   for b in range(4):
    trace.append([a,c,b,a^c^b,((a&1)+(c&1)+(b&1))//2,((a>>1)+(c>>1)+(b>>1))//2])
 assert len(trace)==64 and all(pa<=1 and pb<=1 for a,c,b,o,pa,pb in trace)
 # Basis order uses ordinary, non-factorial coefficients. The tensor product
 # decomposition below comes from coefficients of L(R)S(R), not from derivative
 # product conventions. Positive common-denominator proof is deliberately loose.
 classes={'Gamma':DU([x for row in frame['Gamma_two'] for i,x in row]),
  'N_plane5':DU(N),'G8_plane5_with_one':DU([G,'1']),
  'nu_plane5_with_one':DU([nu,'1'])}
 d=math.prod(int(v['D']) for v in classes.values())
 u=math.prod(int(v['U']) for v in classes.values())*(1<<20)
 committed=max(d.bit_length(),u.bit_length());temporary=3*committed+1
 assert committed<=384 and temporary<=1153,(classes,committed,temporary)
 payload={'physical_R0_variables':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],
  'ell_definition':'5*(sA-rA)','ordinary_plane5_basis':['1','rA','ell','rA^2','rA*ell'],
  'trace_basis_masks':['1','tauA','tauB','tauA*tauB'],
  'source_composition':'Xp_native=Res4(N^t)*Res4(S12_by32(G8)); bottom8 rows zero',
  'normalizer_original_row_column_order':'N[source12][physical12]; transpose12 only, never transpose mu4',
  'original27_X_source_columns':piv,'native23_pivot_recipes':recipes,
  'Gamma_two':frame['Gamma_two'],'P27_columns56':frame['P27_columns56'],'Pi3_columns56':frame['Pi3_columns56'],
  'same_point_M0_Gamma_frame':'literal equality with accepted two_Gamma_on56_actual_X2; constant resolution/gauge transports already bound by CFULL source',
  'cached_N_12_by12_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE':N,
  'canonical_G8_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE':G,
  'canonical_nu_two_plane5':nu,'canonical_G8_order':canonical['G8_order'],
  'universal_trace_product_a_c_inputmask_output_nuA_nuB':trace,
  'full_L2S0_retained_in_both_second_recipes':True,
  'trace_carry_rule':'Each Lk and Sk is coefficient of its full Res4 matrix with varying monics; Lrr/Lre contain N0*nu2 collision terms. Direct product uses universal64 routing, retaining the same carries.',
  'no_S0_zero_weighted_kernel_shortcut_for_individual23':True,
  'source_positive_classes_D_U':classes,'source_path_count_upper':1<<20,
  'coarse_committed_bits':committed,'coarse_temporary_bits':temporary,
  'outside_scope_opaque':['all OLD79/nonunit quotients','outerOLD4 pivots604/613/622/631','other second directions of N','pureplane r^2ell third','617 upper images'],
  'no_signed_L2_S2_Xp2_or_Gamma_Xp2_response_evaluated':True,
  'future_value_worker_prepared_or_authorized':False}
 # Counts include input JSON+parsed copies, all 23x2x2 target outputs and their
 # at most4 term pieces, 3 current C20 vectors, serialized output and metadata.
 n=count(payload)
 future=2*n+count([N,G,nu,frame['Gamma_two']])+2*(92+368)+3*80+6000+3000
 assert future<60000,(n,future)
 payload['prospective_value_whole_numeric_bound']=future+2*count([future,n])+20
 assert payload['prospective_value_whole_numeric_bound']<60000
 source=count([reports,canonical,nr,ne,op,chart,frame,gamma,origin,N,G,nu,trace,classes,payload,pins])+20000
 assert source<150000 and raw+4_000_000<32*1024*1024
 ep='factored-native23-Gamma-pivot-plane2-source-and-full-L2S0.json';dump(ep,payload)
 report={'id':STEM,'all_passed':True,'literal_current27_and_native23_order_bound':True,
  'canonical_C68_G8_monic_original_N_branch_and_full_trace_carry_bound':True,
  'individual_native_L2S0_terms_retained':True,'new_physical_gauge_or_response_values':False,
  'source_owned_numeric_bound':source,'prospective_value_numeric_bound':payload['prospective_value_whole_numeric_bound'],
  'source_positive_committed_bits':committed,'source_positive_temporary_bits':temporary,
  'hard_guards':[384,1153],'whole_RAW_read_plus_output_bound':raw+4_000_000,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start,
  'scope':'SOURCE expression only. No native pivot response, honest quotient, upper617, HF or length value.'}
 dump('report.json',report);dump('phase.json',{'phase':'FINAL_SOURCE_ONLY','report_sha256':sha(out/'report.json')})
 print(json.dumps({k:report[k] for k in ['all_passed','source_owned_numeric_bound','prospective_value_numeric_bound','source_positive_committed_bits','source_positive_temporary_bits','seconds']}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
