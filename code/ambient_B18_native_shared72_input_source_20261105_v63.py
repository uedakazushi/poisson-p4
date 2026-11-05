#!/usr/bin/env python3
"""SOURCE-only shared75 formal auxiliaries ->72 actual native kernel responses."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0
 def read(path):
  nonlocal raw
  p=ROOT/path;pins[path]=sha(p);raw+=2*p.stat().st_size;assert raw<32*1024*1024;return json.loads(p.read_text())
 def ck(n):
  nonlocal bits
  bits=max(bits,abs(n).bit_length());assert bits<=384 and 3*bits+1<=1153;return n
 def dump(name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
 bp='reports/AMBIENT_B18_NATIVE_KERNEL105_FUNCTIONAL_SOURCE_20261105_v63/'
 br=read(bp+'report.json');assert br['all_passed'] and pins[bp+'report.json']=='be68ef59b96e555a8f9e8b0f6a5592648072f2904cb9b1fc943e8dd630ba5d95'
 be='factored-two-by105-native-kernel-functional-source-and-prior.json';b=read(bp+be);assert pins[bp+be]==br['evidence_sha256'][be]
 mp='reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/'
 mr=read(mp+'report.json');assert mr['all_passed'] and pins[mp+'report.json']=='faec4ff020e7f1af70cbf58d825ee09c85aeb31000d540b0f261c48ade511754'
 me='all617-native-constant-sources-factor-through-fixed105-kernel.json';m=read(mp+me);meta=read(mp+'metadata-distinct-basis-support-and-source-uniqueness.json')
 assert meta['all619_input_source_vectors_have_unique_indices'] and pins[mp+'metadata-distinct-basis-support-and-source-uniqueness.json']=='c07cdf2f87eb389532ee1dc3156f443cb86177d889c00234d705d024943d0faa'
 opaque=set(meta['opaque_second_case_indices12']);allmaps=m['all617_source_maps'];selected=[r for r in allmaps if r['free_case619'] not in opaque];assert len(selected)==605
 ids=sorted({i for r in selected for i,q in r['kernel105_coordinates']});assert ids==meta['opaque_free605_basis_indices'] and len(ids)==72
 restrictions=[b['kernel105_restrictions'][i] for i in ids];rawids=sorted({src for rec in restrictions for src,q in rec['original_native_source_combination']});assert len(rawids)==75
 actualraw={i for r in selected for i,q in r['literal_native128_constant_combination']};assert len(actualraw)==72 and set(rawids)-actualraw=={4,8,12}
 intrinsic=[rec for rec in b['intrinsic112_original_Koszul_columns'] if rec['native128'] in rawids];assert len(intrinsic)==75
 p=b['same_J18_with_r2ell_primitive_input'];assert p['N_known_positions']==[0,1,2,3,4,5,6,15,16]
 assert {row for r in intrinsic for row,gen,sign in r['Koszul_entries']}==set(range(8))
 pruned={'basis18':p['basis18'],'N_known_positions':p['N_known_positions'],'N_known_first8_by12_by9_by4':p['N_known12_by12_by9_by4'][:8],
  'G8_C72':p['G8_C72'],'nu_two_J18':p['nu_two_J18'],'Gamma_two56':p['Gamma_two56'],'D_two_by10_by56':p['D_two_by10_by56']}
 basis=[tuple(v) for v in p['basis18']];idx={v:i for i,v in enumerate(basis)};table=[[idx.get(tuple(x+y for x,y in zip(a,c))) for c in basis] for a in basis]
 nus=[[F(q) for q in row] for row in p['nu_two_J18']]
 assert all(table[i][j] is None for a in nus for c in nus for i,x in enumerate(a) if x for j,y in enumerate(c) if y)
 # No native response arithmetic: onlypropagate a previously source-proved
 # COMMON denominator/absolute numerator through sparseconstantmap.
 D0=ck(int(b['common_D']));U0=ck(int(b['common_U_with_paths2pow28']));assert max(D0.bit_length(),U0.bit_length())==341
 assert b['committed_prior_bits']==341 and b['temporary_prior_bits']==1024
 before=count([br,b,mr,m,meta,selected,ids,restrictions,rawids,intrinsic,pruned,table,nus,pins])+20000;assert before<150000
 dump('phase-before-positive-coordinate-composition.json',{'no_response_arithmetic':True,'basis_common_prior341_temp1024':True,'source_whole_upper':before})
 envelopes=[];mapping=[];maxcase=0;nnz=0
 for rec in selected:
  coeff=rec['kernel105_coordinates'];nnz+=len(coeff);assert all(i in ids for i,q in coeff);den=1
  for i,q in coeff:den=ck(math.lcm(den,F(q).denominator))
  u=0
  for i,q in coeff:
   q=F(q);u=ck(u+ck(abs(q.numerator)*(den//q.denominator)))
  if u:bound=[str(ck(D0*den)),str(ck(U0*u))]
  else:bound=['1','0']
  maxcase=max(maxcase,max(map(lambda q:int(q).bit_length(),bound)));envelopes.append([rec['free_case619'],[bound,bound]])
  mapping.append({'case619':rec['free_case619'],'kernel105_coordinate_pairs':coeff})
 assert nnz==120 and len(mapping)==605
 globalbits=max(341,maxcase);assert globalbits<=384 and 3*globalbits+1<=1153
 data={'SOURCE_ONLY_all_native_responses_UNRUN':True,'native_primitive_known_cache':pruned,
  'raw75_source_indices':rawids,'raw75_original_Koszul_program':intrinsic,
  'kernel72_original105_indices':ids,'kernel72_restrictions':restrictions,
  'all605_constant_native_coordinate_recipes':mapping,'all605_positive_native_sum_envelopes_same_two_channels':envelopes,
  'formal_auxiliary_rule':'Eachraw75 scalarrowfromknownpartialN is a FORMALextension, notactualunknownjet. Extras4,8,12 have nonzeroS0. Onlyfixed72kerS0 restrictions,then605constantmaps,haveactualprojectedGamma3+Dsecondmeaning.',
  'known_partialNbase9_formula':'N0+fouractivezfirst+Nr+Nell+Nrr+Nrell; missingN2mixed/N3remainformal and cancelONLYkerS0/pureplaneextraction.',
  'trace_rule':'Res4 in samebasis(1,tauA,tauB,tauAB), scalarordinarycoefficients; allnu_i*nu_j supportsareoutsideJ18. Rejectunsupportedbase before scalarproducts; keepall single-nu carry.',
  'kernel_functional_formula':'Gamma[L0S3r2ell+LrS2rell+LellS2rr+LrrSell+LrellSr]+Dnonzz[L0S2+uniqueLfirstSfirst]; rr once.',
  'no_OLD_FIRST_outer_or_honest_responses_or_605upper_total':True,
  'native_basis_common_D':str(D0),'native_basis_common_U':str(U0),'native_basis_prior':341,'native_basis_temporary':1024,
  'native605_coordinate_composition_prior':globalbits,'native605_coordinate_composition_temporary':3*globalbits+1,
  'fixed75to72to605_linear_maps_only_no_rank_search':True}
 n=count(data);workspace={
  'parsedN_first8_by12_by9_by4':8*12*9*4,
  'G576_nu36_Gamma18_D1120':576+36+18+1120,
  'streamed_one_intrinsicS_two_rows_by72_index_value':2*72*2,
  'streamed_partialNSc_seven_targets_by72_index_value':7*72*2,
  'two75formal_aux_values_index_value':2*75*2,
  'two72actual_restricted_values_index_value':2*72*2,
  'native605_two_values_andcaseID':3*605,
  'algebra_base_table_plusmonomialmetadata':324+18*10+1000,
  'source75_and72_and605_mapping_parse':count([intrinsic,restrictions,mapping]),
  'scalar_wrappers_current_accumulators':2000}
 derived=sum(workspace.values());live=2*(n+1000)+derived+12000+4000;assert live<60000 and live<150000
 rawfuture=2_000_000+live*(2*math.ceil(globalbits*math.log10(2))+42);assert rawfuture<32*1024*1024
 data.update({'input_numeric_before_finalmeta':n,'parsed_workspace_book':workspace,'future_numeric_forecast':live,'future_RAW_forecast':rawfuture,
  'guard_extra_NEWstate_copy_reserve':12000,'output_reserve':4000,'metadata_reserve_each_inputcopy':1000,
  'fresh_immutable_runtime_worker_custody_and_independent_SOURCEpeer_and_ROOT_GO_REQUIRED':True})
 ep='native-shared75-formal-to72-actual-to605-source-input.json';dump(ep,data)
 owned=count([br,b,mr,m,meta,selected,ids,restrictions,rawids,intrinsic,p,pruned,basis,idx,table,nus,envelopes,mapping,data,workspace,pins])+10000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_positive_and_linearity_only':True,'native_actual_response_values_UNRUN':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep),'phase-before-positive-coordinate-composition.json':sha(out/'phase-before-positive-coordinate-composition.json')},'source_owned_upper':owned,'source_RAW_upper':raw+sum(p.stat().st_size for p in out.iterdir()),'native_basis_prior':341,'native605_mapping_prior':globalbits,'native605_mapping_temporary':3*globalbits+1,'future_numeric_forecast':live,'future_RAW_forecast':rawfuture,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'source':owned,'forecast':live,'basis_prior':341,'mapping_prior':globalbits}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)
