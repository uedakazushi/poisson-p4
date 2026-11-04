#!/usr/bin/env python3
"""ONE compact8x98 generic QQ(z) rank/minor gate; no uniform-rank claim."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';pf=base/'S61_NONSQUARE_P100_PURE8_RANK_SOURCE_PREFLIGHT_20261104_v63';pr=read(pf/'report.json');pc=read(pf/'certificate.json');assert sha(pf/'certificate.json')==pr['certificate_sha256']and pr['known_unit_row_clear_polynomial_terms_upper']==23618
 sp=base/'S61_NONSQUARE_P100_PURE97_SPLIT_20261104_v63/preflight.json';s=read(sp)
 rp=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/pure98.json';g=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63'
 rowpaths=[rp]+[g/f'pure98-row{i}.json'for i in range(7)];rd=[read(p)for p in rowpaths];rows=[d['values98']for d in rd];assert all(d['target122_columns']==rd[0]['target122_columns']for d in rd)
 for p in rowpaths:assert sha(p)==pr['input_sha256'][str(p.relative_to(ROOT))]
 powers=[s['q5_common_denominator_powers5']]+[p['common_denominator_powers5']for p in s['seven_profiles']]
 reserve=pc['RF_input_terms']+10;assert reserve==18203 and sum(sum(v)for v in powers)<400
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp';paths=[Path(__file__),cpp,utility,pf/'report.json',pf/'certificate.json',sp]+rowpaths
 text=[str(reserve)]+[' '.join(map(str,p))for p in powers]
 for row in rows:
  for v in row:
   for k in['numerator','denominator']:text.append(str(len(v[k])));text.extend(f'{m[0]} {c}'for m,c in v[k])
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<3*1024*1024
 labels={'source_order':['actual_q5H_ROW0_minus_full_cubic_H1_lift']+[{'actual135':18,'physical_light':j,'source319_index':174+i}for i,j in enumerate([0,1,2,3,4,6,7])],
         'target122_labels98':rd[0]['target122_columns'],'known_curve_factors_order':['z','z+20','z+60','203z+2700','z+12'],'row_LCD_powers':powers,
         'minor_scaling':'Selected minor of original RF rows equals native determinant divided by product of selected rowLCDs, all whole-curve units.',
         'algorithm':'Increment selected square via D*d-c*adj(S)*b; update adjugate by exact polynomial divisions and verify both products. Never Schur-update full8x98.',
         'degree_bounds':{'selected_determinant':233,'adjugate':212},'retained_polynomial_terms_bound_below150k':True}
 (out/'source-labels.json').write_text(json.dumps(labels,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'rank-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 nd=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode}
 complete=r.returncode==0 and'generic_rank'in nd
 rank=nd.get('generic_rank')if complete else None
 if complete:
  assert nd['all784_known_unit_row_clearings_verified']and nd['all_selected_adjugate_left_right_products_verified']
  assert nd['rank_upper_from_number_of_rows_if8']or nd['all784_span_upper_identities_verified']
  assert nd['peak_live_polynomial_terms_with_literal_input']<=150000
  nd['selected_original122_target_labels']=[rd[0]['target122_columns'][j]for j in nd['selected_columns']]
  # Keep native bytes immutable; labels are recorded in this separate binding.
  (out/'selected-source-binding.json').write_text(json.dumps({'selected_sources':[labels['source_order'][i]for i in nd['selected_rows']],'selected_local98_columns':nd['selected_columns'],'selected_original122_target_labels':nd['selected_original122_target_labels'],'selected_rowLCD_powers':[powers[i]for i in nd['selected_rows']]},indent=2)+'\n')
 report={'id':ID,'complete':complete,'generic_rank_saved_pure8x98':rank,'generic_rank_saved_pure7x97':rank-1 if complete else None,'full_original_quartic_pointwise_lower_on_D_selected_minor':141+rank if complete else None,
         'whole_P100_lower142_unchanged':True,'selected_minor_unknown_nonunit_zeros_NOT_excluded':True,'no_uniform149_no_full239_no_other_Gi_or_automatic_nextjob':True,
         'native_audit':{k:v for k,v in nd.items()if k not in['determinant_nonzero','selected_square','final_adjugate','determinant_steps','row_span_coordinates_numerator8xr']},
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
