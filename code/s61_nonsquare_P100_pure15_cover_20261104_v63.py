#!/usr/bin/env python3
"""ONE fixed56 replacement-minor/gcd cover, stopping only at unit gcd."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE15_COVER_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 pf=base/'S61_NONSQUARE_P100_PURE15_COVER_PREFLIGHT_20261104_v63';pr=read(pf/'report.json');pc=read(pf/'certificate.json');assert sha(pf/'report.json')=='26333b08e690c71504f19070bbf9f406c70053782c61bb65c7f4e55a166dab6e'and sha(pf/'certificate.json')==pr['certificate_sha256']
 rb=base/'S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63';rr=read(rb/'report.json');r=read(rb/'native.json');assert rr['generic_rank_saved_pure8x98']==8 and sha(rb/'native.json')==rr['evidence_sha256']['native.json']
 kb=base/'S61_NONSQUARE_P100_K2_PURE7_20261104_v63';kr=read(kb/'report.json');ks=read(kb/'source-binding.json');assert kr['complete'];kpaths=[kb/f'pure98-row{i}.json'for i in range(7)];kd=[read(p)for p in kpaths]
 rp=base/'S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW/report.json';kp=base/'S61_NONSQUARE_P100_K2_PURE7_20261104_v63_REVIEW/report.json'
 assert sha(rp)=='8ef758e0872b7e8b54b478fe82e275daffc24605f1b92a80ff17d472e373ae39'and sha(kp)=='ea4e71d95e47a71069392aff8ac64170ab6e5681e7de2cbc8485b41aa21e77b3'
 for p in kpaths:assert sha(p)==kr['seven_row_sha256'][p.name]==pr['input_sha256'][str(p.relative_to(ROOT))]
 cs=r['selected_columns'];assert cs==pc['fixed_target_local98_columns8'];powers=[v['new_known_unit_LCD_powers5']for v in pc['profiles7']]
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp'
 paths=[Path(__file__),cpp,utility,pf/'report.json',pf/'certificate.json',rb/'report.json',rb/'native.json',rb/'source-labels.json',rb/'selected-source-binding.json',rp,kb/'report.json',kb/'source-binding.json',kp]+kpaths
 # Conservative reserve includes full old18193 and previous K2serialized15464,
 # even though native reads only the selected8 new entries per source row.
 reserve=18193+15464+r['serialized_terms']+256;assert reserve==44870
 text=[str(reserve)]
 def poly(p):text.append(str(len(p)));text.extend(f'{m[0]} {c}'for m,c in p)
 def rf(v):poly(v['numerator']);poly(v['denominator'])
 poly(r['determinant_nonzero'])
 for M in[r['selected_square'],r['final_adjugate']]:
  for row in M:
   for p in row:poly(p)
 for ex in powers:text.append(' '.join(map(str,ex)))
 for d in kd:
  for c in cs:rf(d['values98'][c])
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<3*1024*1024
 binding={'old_selected_source_binding':read(rb/'selected-source-binding.json'),'new_corrected_K2_source_binding':ks,
          'fixed_old_selected_target98_columns':cs,'fixed_original122_target_labels':[kd[0]['target122_columns'][j]for j in cs],
          'new_sources_in_order':[{'kernel':2,'actual_nonpivot135':23,'physical_light':d['physical_light'],'known_unit_rowLCD_powers5':ex}for d,ex in zip(kd,powers)],
          'minor_formula':'Replace old selected rowpositioni (keep allother roworder) by the newLCD-cleared K2row; determinant=sum_k vclear[k]*oldAdj[k,i]. No extra permutation sign.',
          'row_unit_effect':'Old and new LCDs are known wholeP100units. The polynomial minors differ from actual projected source minors only by these rowunits.',
          'no_new_source_coefficient_NF_row_or_rank':True,'literal_reserve':reserve,'allowed_fixed_candidates_only':56,'stop_only_when_exact_gcd_is_unit':True}
 (out/'source-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'cover-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:run=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 nd=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':run.returncode};complete=run.returncode==0 and nd['all_passed']
 unit=complete and nd['stopped_at_knownunit_gcd'];gd=max((m[0]for m,c in nd['final_gcd']),default=-1)if complete else None
 if complete:
  assert nd['literal_sum_Ai_actualMinor_i_equals_E_times_gcd_verified']and nd['all56_selected_LCD_clearings_verified']and nd['peak_live_polynomial_terms_with_literal_reserve']<=150000
  assert nd['computed_fixed_replacements']==56 or unit
 cert={'whole_P100_pointwise_original_quartic_lower149':unit,'chosen_minor_family_gcd_degree':gd,
       'proof_if_unit':'The literal actual-minor polynomial identity sums known projected8minors to E, a product only of five known wholecurveunits. One8minor is nonzero at eachpoint, giving8pure directions plus accepted141block.',
       'if_nonunit_scope':'Unknowngcd only describes the unresolved finite chosen-minor cover; it does not prove actualrankfailure. Higher row replacements or other sources may rescue. No automatic extension.',
       'unknown_factors_inverted':False,'no_fixed149_original_minor_assumed':True,'actual_minor_and_primitive_Bezout_in_native':True,'existing_whole142_and_generic149_remain_valid':True,'no_full239_or_globalFULLU_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 audit={k:v for k,v in nd.items()if k not in['old_minor_strip','candidates','final_gcd','literal_primitive_Bezout_contributors','common_actual_minor_unit_polynomial','actual_minor_polynomial_Bezout_contributors']}
 report={'id':ID,'complete':complete,'whole_P100_quartic_pointwise_lower149':unit,'chosen_cover_unknown_gcd_degree':gd,'native_audit':audit,'no_rank_search_newNF_or_automatic_two_row_replacements':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
