#!/usr/bin/env python3
"""Only saved-row univariate gcd/Bézout; no source coefficient evaluation."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import sha,read
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE98_ROW_GCD_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 source=ROOT/'reports/S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63';cpp=Path(__file__).with_suffix('.cpp');paths=[Path(__file__),cpp,source/'report.json',source/'pure98.json',source/'source-binding.json']
 sr=read(paths[2]);data=read(paths[3]);assert sr['complete']and sha(paths[3])==sr['evidence_sha256']['pure98.json'];assert data['all98_denominators_known_whole_curve_units']
 text=[]
 def poly(p):text.append(str(len(p)));text.extend(f'{m[0]} {v}'for m,v in p)
 for r in data['values98']:poly(r['numerator']);poly(r['denominator'])
 factors=[('z',[[[1],'1']]),('z_plus20',[[[1],'1'],[[0],'20']]),('z_plus60',[[[1],'1'],[[0],'60']]),('203z_plus2700',[[[1],'203'],[[0],'2700']]),('z_plus12',[[[1],'1'],[[0],'12']])]
 for name,p in factors:text.append(name);poly(p)
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'gcd-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20);subprocess.run([str(exe),str(ip),str(out/'native.json')],check=True,timeout=max(1,55-(time.monotonic()-start)))
 native=read(out/'native.json');g=native['remaining_unknown'];degree=max(m[0]for m,v in g);unit=native['row_unimodular_on_whole_P100'];assert native['literal_sum_Ai_Ni_equals_gcd']
 cert={'actual_source_row':str(paths[3].relative_to(ROOT)),'actual_target122_labels':data['target122_columns'],'curve_unit_factors5':factors,
       'exact_row_ideal':'Localized row-entry ideal equals the principal ideal of the saved remaining_unknown polynomial; all98dividebygcd andliteralextendedgcdarechecked.',
       'row_unimodular_on_whole_P100':unit,'regular_rowdual_source':'native.json global_regular_row_dual; coefficients Ai*Di/G, literal sum Ai*Ni=G',
       'pointwise_rank_lower142_on_whole_P100':unit,'if_not_unit_scope':'Only this row vanishes on its remaining polynomial finite zero set; other generators can rescue. No actual quartic failure claim.',
       'no_column_GL_or_further_Gaussian_or_new_source_values':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'exact98_numerator_gcd_and_all98_divisions_verified':True,'literal_extended_gcd_verified':True,'contributor_positions':[x['position']for x in native['literal_Bezout_contributors']],
         'row_unimodular_on_whole_P100':unit,'remaining_unknown_degree':degree,'native_stored_terms':native['stored_terms'],
         'pointwise_fullquartic_lower142_on_whole_P100':unit,'no_rank98_full239_or_automatic_nextjob':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'native_input_sha256':sha(ip),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
