#!/usr/bin/env python3
"""ONE554 weighted native SOURCE; no functional response or outer pairing."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_ONE554_FRAME_SOURCE_20261105_v63'
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
 assert set(z)==set(names);return z
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 bp='reports/AMBIENT_B18_NATIVE_ONE549_FRAME_SOURCE_20261105_v63/ONE549-J18-with-r2ell-native-full-S0-column-and-formal-N2N3-frontier.json';b=read(bp)
 pp='reports/AMBIENT_B18_NATIVE_ONE549_SOURCE_POSITIVE_20261105_v63_REVIEW/report.json';peer=read(pp);assert peer['all_passed'] and pins[pp]=='b02f0d2118aab1499477d6c05603e8e0a7f8738b83cc702795fd96224a234e18'
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_SOURCE_20261105_v63/report.json';sr=read(sp);assert sr['all_passed'] and pins[sp]=='1a5b1726d6592dbd67fe647cfbd9e65bd36997e615aa8f66696767d6a43ccc15'
 se=next(iter(sr['evidence_sha256']));s=read(str(Path(sp).parent/se));assert pins[str(Path(sp).parent/se)]==sr['evidence_sha256'][se]
 ap='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json';a=read(ap,['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'][554]
 source=[[289,'25/4'],[379,'-1'],[643,'29/4'],[645,'25'],[651,'-25/2'],[670,'1']];assert a['original768_combination']==source
 aa=a['actual_a1_27_times9'];assert len(aa)==29
 plane=[[i,q] for i,q in aa if i%9 in [7,8]];assert plane==[[52,'-175/3'],[53,'-377/60']]
 mp='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';mu=read(mp,['native_S0_columns'])['native_S0_columns']
 op='reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json';S=read(op,['native_S12_by32','basis','structure16'])
 # Replay fixed full48-component point inputs. All4 full columns happen to be
 # empty; this fact is checked rather than inferred from G8 or M0X0.
 native=[]
 for original,q in source:
  if original<640:continue
  j=original-640;col,mask=divmod(j,4);point=mu[j];assert point==[]
  terms=[['G8',row,g,sign] for row,c,g,sign in S['native_S12_by32']['koszul28_in_first8_rows'] if c==col]
  assert len(terms)==2
  native.append({'original768':original,'native128':j,'constant_source_coefficient':q,'S32_source':col,'input_trace_mask':mask,'full_point_S0_column48':point,'intrinsic_source_entries':terms})
 assert [(v['native128'],v['constant_source_coefficient']) for v in native]==[(3,'29/4'),(5,'25'),(11,'-25/2'),(30,'1')]
 assert mu[4]==[[11,'-1/2']] # actual outerpivot644 inputmask0 differs.
 # Full current covectors/a1 source in the external OLD+FIRST lane, no
 # branch/complement change. No OLD/FIRST coefficient is evaluated here.
 assert s['Gamma_two']==b['Gamma_two_current56'] and s['D_two_by10_by56']==b['D_two_by10_by56']
 # Current producer field names are validated explicitly before coefficient use.
 expected_a1=s.get('actual_full_a1_27times9',s.get('full_a1_27times9'))
 assert expected_a1==aa
 vp='reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json';v=read(vp)
 dp='reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json';direct=read(dp);assert direct['all_passed'] and pins[dp]=='10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92'
 k=v['original768_native23_sources'].index(644);assert k==1 and any(q!='0' for row in v['Gamma_Xp2_native23_twochannels_rr_rell'][k] for q in row)
 ev={'selected_case619':554,'literal_full_original768_source':source,'all29_a1_kept_UNEVALUATED':aa,'only_plane_a1_pivot5_original644':plane,
  'current_Gamma_two56':b['Gamma_two_current56'],'current_D_two_by10_by56':b['D_two_by10_by56'],'D_ordinary_pair_order10':b['D_ordinary_pair_order10'],
  'J18_with_r2ell_known_native_function_and_N_jets_source':{'path':bp,'sha256':pins[bp],'independent_peer':pp,'peer_sha256':pins[pp]},
  'parameterized_four_native_source_records':native,
  'weighted_FULL_Res4_S0_times_cN_exactzero_certificate':{'all48_components':'0','method':'sum of the4 literal fullinputmask columns, each verified empty in originalMUsource; constants29/4,25,-25/2,1','not_M0X0_or_G8_entrywise_inference':True},
  'outer_pivot644_mask0_FULL_S0_column_nonzero':mu[4],
  'outer_full_GammaXp644_rr_rell_actual_cached_twochannels':v['Gamma_Xp2_native23_twochannels_rr_rell'][k],
  'outer_source_formula':'SUBTRACT GammaXp644_rr*(-377/60) + GammaXp644_rell*(-175/3). KeepFULLL2S0; no signed pairing evaluated.',
  'native_known_coefficient_formula':'Gamma[L0S3_r²ell+LrS2_rell+LellS2_rr+LrrSell+LrellSr]+sum10Dpair[L0S2_pair+uniqueordinary_first_decompositions]; S means WEIGHTED Res4(S)*cN for all4 inputmasks/constants.',
  'unknown_mixedN2_N3_remain_FORMAL_then_factor_through_WEIGHTED_S0_or_zfactor':True,
  'no_individual_G8_S0_zero_assumption_or_old648_zero_reuse':True,
  'all_varying_nu_Res4_operator_first_and_planeL2_carries_retained':True,
  'external_OLD_and57_FIRST_lane':{'source_report':sp,'source_sha256':pins[sp],'separate_source_positive_peer_still_needed':True},
  'no_native_S2_S3_outer_pairing_FIRST_or_whole_response_computed':True,
  'whole_positive_height_owner_code_custody_and_ROOT_GO_still_required':True}
 owner=count([b,peer,sr,s,a,aa,mu,S,native,v,direct,ev,pins])+20000;assert owner<150000 and raw+2_000_000<32*1024*1024
 ep='ONE554-four-native-weighted-S0-currentframe-nonzero644-outer-source.json';dump(ep,ev)
 rep={'id':STEM,'all_passed':True,'SOURCE_ONLY_no_signed_responses':True,'full_weighted_native_S0_certificate_verified':True,
  'nonzero_outer644_prior_retained_fullL2S0':True,'shared_J18_and_true_Njet_branch_bound_to_independent_peer':True,
  'whole_source_numeric_upper':owner,'whole_RAW_upper':raw+2_000_000,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'RAW':raw+2_000_000,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
