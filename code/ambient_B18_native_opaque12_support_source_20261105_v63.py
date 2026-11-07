#!/usr/bin/env python3
"""Literal support selection from frozen617->105; no scalar response arithmetic."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE_OPAQUE12_SUPPORT_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):return int(all(x.lstrip('-').isdigit() for x in v.split('/')) and len(v.split('/'))<=2)
 return int(isinstance(v,(int,float)))
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));start=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={};raw=0
 def read(p,h=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024
  if h is not None:assert pins[p]==h
  return json.loads(q.read_text())
 s=read('reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/all617-native-constant-sources-factor-through-fixed105-kernel.json','1b2e135146084f7e99a36a200d0d62ff868bf9711384cf084a40f8c39a8f6774')
 assert read('reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/report.json','faec4ff020e7f1af70cbf58d825ee09c85aeb31000d540b0f261c48ade511754')['all_passed']
 sh=read('reports/AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63/native-shared75-formal-to72-actual-to605-source-input.json','b51fa2bb95a98383987c3d718c6b87ceb25b15a6fdf6edd208d6a3e916d96da5')
 b=read('reports/AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63/constant-kernel105-basis-unit23-minor-and-native-functional-determinacy.json','c1d117db77a0d295652cd53fff3380ab3c48abb7af14bb213010d3eb303e68ce')
 assert read('reports/AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63/report.json','15b6ff370efd02939a1fac290ca6c24c1237ab3046f6a82984d22b1071e8d534')['all_passed']
 assert read('reports/AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW/report.json','0384d2623fab486eb0dbd1a4df239bd04e6a7ff778f7464c454ec324c541683e')['all_passed']
 assert read('reports/AMBIENT_B18_UPPER20_INDEPENDENT_POSTCHECK_V2_20261105_v63/report.json','435cbe8574ef62f3de1e62182cd131ac032944ddc324f4010449344ba73fc5e7')['all_passed']
 assert s['all617_literal_full48S0cN_zero_verified'] and s['all617_full128_reconstructions_verified']
 cases=[579,580,593,594,597,598,605,606,611,613,615,616];sel=[r for r in s['all617_source_maps'] if r['free_case619'] in cases];assert [r['free_case619'] for r in sel]==cases
 basis=b['explicit105_kernel_basis_columns'];used=sorted({i for r in sel for i,q in r['kernel105_coordinates']});old=set(sh['kernel72_original105_indices']);new=sorted(set(used)-old);reused=sorted(set(used)&old)
 rawall=sorted({i for k in used for i,q in basis[k]['column']});rawnew=sorted({i for k in new for i,q in basis[k]['column']});assert len(used)==46 and len(new)==30 and len(reused)==16 and len(rawall)==53 and len(rawnew)==37
 assert sum(len(r['kernel105_coordinates']) for r in sel)==109
 records=[{'kernel105':k,'explicit_original128_source_combination':basis[k]['column'],'kind':basis[k]['kind']} for k in used]
 e={'SOURCE_metadata_only_no_new_response_or_positive_fold':True,'opaque12':cases,'all12_corrected_source_maps':sel,'kernel46_literal_restrictions':records,'kernel46_ids':used,'kernel30_absent_shared72':new,'kernel16_already_shared72':reused,'raw53_all_kernel46':rawall,'raw37_for_new_kernel30':rawnew,'all109_constant_coordinate_terms_kept':True,'scalar_mask0_raws_retained_FORMAL_until_restriction':[i for i in rawall if i%4==0],
 'no_pointzero_inference_from_scalar_mask_or_M0cancellation':True,'all12_fullS0_and128_factor_proof':'Inherited literally fromfaec+1b2 andindependent0384/full105diagonal-complement proof; thispacketmerelyselectsunchangedrecords.',
 'current_comparison_frame_unchanged':'SamecurrentGamma/D/fullPsi/currentM0/nativeNtranspose andvaryingtracebasis; nofield/gaugechange orcanonicalM1zero.',
 'source_correction_beforeinterpretation':'Native12mapsalreadyincludesubtract563609; neveraddbasisnativeagain. CorrespondingOLD/FIRST/outercorrectionremainsseparate.',
 'not_a_response_cache_certificate':'16overlapkernelsareSOURCE/FINE overlap only; notall72valueswerecomputedby20job. Only15actual20kernels+10explicitSOURCEzerosand554previousresponseexist.',
 'future_native_requirements':'Anynewresponsejobmustpruneactualneed,bindall48S0kernelrestrictions,retainformalrawpointtermsandallnu/firstcross,derivepositivewholeheight/ownershipandfreshimmutableworker/custody/ROOTGO.',
 'honest_OLD_frontier_unchanged':'12opaqueOLDfunctionals/24SECOND remainopaque; canonicalnativecachedoesnotassignhonestquotientcoefficients.',
 'inherited_constant_only_prior72':'NoQQarithmeticorcoefficientcomparisonisperformedhere;allselectedQQstringscopiedbyte-for-bytefromacceptedsource.',
 'no_HF2_length_H0_claim':True}
 owned=2*count([s,sh,b])+2*count(e)+4000;assert owned<60000
 ep='twelve-opaque-native46-basis-new30-and-raw53-literal-support.json';(out/ep).write_text(json.dumps(e,separators=(',',':'))+'\n');assert raw+sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 r={'id':STEM,'all_passed':True,'metadata_SOURCE_only':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'source_owned_upper':owned,'source_RAW_charged':raw,'kernel_counts':[46,30,16],'raw_counts':[53,37],'constant_source_only_no_arithmetic':True,'no_signed_native_OLD_or_upper_values':True,'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(r,separators=(',',':'))+'\n');print(json.dumps({'report_sha256':sha(out/'report.json'),'owned':owned,'RAW':raw,'kernel_counts':[46,30,16],'raw_counts':[53,37]}))
if __name__=='__main__':main()
