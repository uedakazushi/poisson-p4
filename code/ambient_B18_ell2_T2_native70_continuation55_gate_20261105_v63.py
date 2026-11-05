#!/usr/bin/env python3
"""Source/count/permanent bounds for 55 uncached images and fixed 4x70 rank."""
import itertools,json,math
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_INPUT_GATE_20261105_v63'
def main():
 pins={}
 def read(p):pins[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 b='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63';r=read(b+'/report.json');d=read(b+'/complete70-cached-input-source-and-positive-whole-image-bounds.json');assert r['all_passed'] and r['future_numeric_bound']==51263
 live=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_LIVE_BOUND_20261105_v63/report.json');assert live['all_passed'] and live['actual_worker_owned_bound']==49734
 p='reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63';pr=read(p+'/report.json');old=read(p+'/actual-native70-lower-prefix-and-witness-if-any.json');assert pr['all_passed'] and pr['computed_prefix_count']==15 and old['nonzero_witness']['case619']==563
 peer=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WITNESS563_20261105_v63_REVIEW/report.json');assert peer['all_passed'] and pins[p+'/report.json']['sha256'] in json.dumps(peer)
 prefix=old['computed_prefix'];assert [v['case619'] for v in prefix]==list(range(549,564)) and len(prefix)==15
 rest=[r['case619'] for r in d['source70_records'][15:]];assert rest==list(range(564,619)) and len(rest)==55
 DU=[tuple(map(int,v)) for v in d['all70_four_images_global_D_U']];w=list(map(Q,old['nonzero_witness']['actual_lower4_values']));W=[]
 for (D,U),x in zip(DU,w,strict=True):assert D%x.denominator==0;z=x.numerator*(D//x.denominator);assert abs(z)<=U;W.append(z)
 per=[]
 for size in range(1,5):
  for rows in itertools.combinations(range(4),size):
   bound=math.factorial(size-1)*sum(abs(W[i])*math.prod(DU[j][1] for j in rows if j!=i) for i in rows)
   per.append({'rows':list(rows),'permanent_bound_bits':bound.bit_length()})
 assert max(x['permanent_bound_bits'] for x in per)==374<385
 # All signed determinant operations are TEMPORARY integers. Only original
 # QQ entries <=128 bits and literal minor/factor recipes are committed.
 cached=count(prefix);assert cached==255;owned=49734+2*cached+540+100;assert owned==50884<51263
 ev={'source_gate_and_count_addendum_pinned':True,'cached_first15_prefix':prefix,'cached_prefix_count':15,'new_case619_indices55':rest,'rank_start_witness_column_index70':14,'rank_start_witness_case619':563,'row_integer_clearing_D_U':[[str(D),str(U)] for D,U in DU],'known_integer_scaled_witness4':list(map(str,W)),'all15_row_subset_permanent_bounds':per,
 'rank_method':'StartwithknownwitnessWandonefixednonzerorow. For eachcandidatecolumn, testborderedminorsagainsttheselectedbasisandallrowsoutsidefixedpivotrows. A nonzero borderedminor enlargesbothsourcebasisandpivotrowset. Every determinant containsW; allpartialsareboundedbythelistedpermanent. Afterrankfreezes, recheckall70 columnsbythefixedborderedminors. SelectedoriginalQQcolumnsandfixedrow/columnindicesarethecertificate. Coefficientsintheirspanaretheliteraladjugate/determinantrecipe; doNOTcommitlargeevaluatedQQfactorsorlargeintegerminors. Independentpostpeerreplaysallminorsfromtheoriginal70QQvectors.',
 'source_committed_entry_bits':128,'rank_temporary_integer_bits':385,'source_permanent_maximum_bits':374,'cached_prefix_two_owned_copies':2*cached,'extra_rank_working_reserve':540,'extra_profile_and_control_reserve':100,'whole_conservative_owned':owned,'unchanged_whole_numeric_cap':51263,'unchanged_RAW_cap':8305349,'new_values_and_rank_NOT_started':True,'scope':'Onlyremaining55LOWERfourimagesandfull70lowerprojectedrank/span. First15cachevaluesarenotrecomputed.No619upper2/HF2equality/total17claim.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'cached15-plus-new55-and-known-witness-bordered-rank-source-input.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'new_remaining55_values_or_rank_started':False,'source_committed_bits':128,'rank_temporary_bits':385,'positive_rank_permanent_bits':374,'whole_numeric_bound':owned,'unchanged_cap':51263,'whole_RAW_cap':8305349,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned,'rank_temp':374}))
if __name__=='__main__':main()
