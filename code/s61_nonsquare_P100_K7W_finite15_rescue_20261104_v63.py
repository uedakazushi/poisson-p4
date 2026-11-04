#!/usr/bin/env python3
"""ONE NEW C68W border determinant/remainder/wholefiniteH15 gate."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 sb=base/'S61_NONSQUARE_P100_K7W_FINITE15_SOURCE_20261104_v63';sr=read(sb/'report.json');s=read(sb/'source.json');ip=sb/'source-input.txt';assert sr['source_passed']and sr['NEW_source319row226']and sha(ip)==sr['native_input_sha256']and sha(sb/'source.json')==sr['source_sha256']
 pb=base/'S61_NONSQUARE_P100_K7W_FINITE15_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');assert pr['det_degree_upper']==218 and pr['RF_terms']==1057
 hb=base/'S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';hr=read(hb/'report.json');hn=read(hb/'native.json');assert hr['exact_QQ_gcd_degree']==15 and hn['cofactor_QQ_den_and_LC_primeunit_gate']and hn['oldprime_cofactor_gcd_is1']
 nb=base/'S61_NONSQUARE_P100_K7W_PURE98_20261104_v63';nr=read(nb/'report.json');assert nr['complete']and nr['one_literal_K7W_row_only']
 cpp=Path(__file__).with_suffix('.cpp');engine=ROOT/'code/s61_nonsquare_P100_finite15_rescue_20261104_v63.cpp';utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp'
 paths=[Path(__file__),cpp,engine,utility,sb/'report.json',sb/'source.json',ip,pb/'report.json',pb/'certificate.json',hb/'report.json',hb/'native.json',nb/'report.json',nb/'source-binding.json',nb/'pure98-row3.json']
 source={'NEW_bottom_source':s['NEW_bottom_source'],'top8_source_ids':s['case9_top8_source_ids'],'nine_LCD_powers5':s['nine_LCD_powers5'],'nine_oldlocal98_targets':s['nine_oldlocal98_targets'],'nine_original122_targets':s['original122_targets'],
         'not_old_failed_K2W_data':'Thepreviousgenericfixed9engineisreusedasCODEonly. No previousK2Wdet,input,remainderorbottomrowdataisread.',
         'finite_model':'S=Z_(32003)[z]/Hmonic15 finitefree; nofield/reducednessassumption. EntireSbarunitgate; nocurvelocalizerpointdiscard.',
         'allowed_inverses':'OnlyglobalQQpolynomialBareissexactdivisionandmonicHremainder; no oldDiinverse or finite-ringmatrixpivotinverse.',
         'char0_source_units':'Fullknown5linearcurveunitLCDs, originalrow226literalnormalizedC68W; allrawpuretailskept.'}
 (out/'source-binding.json').write_text(json.dumps(source,indent=2)+'\n')
 phase=ROOT/'work/marked-square-v63/K7W-finite15-phase/phase.json';phase.parent.mkdir(parents=True,exist_ok=True);phase.write_text(json.dumps({'status':'source_frozen_before_NEW_math','native_input':str(ip.relative_to(ROOT)),'native_input_sha256':sha(ip),'output':str(out.relative_to(ROOT)),'old_failed_data_NOT_read':True},indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'K7W-finite15-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 d=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode};complete=r.returncode==0 and d['all_passed'];unit=complete and d.get('finite_special_fiber_unit',False)
 degN=max((m[0]for m,c in d['exact_determinant_N']),default=-1)if complete else None
 if complete:
  assert d['all81_source_LCD_identities_verified']and d['literal_N_equals_HQ_plus_R_verified']and degN<=218 and d['peak_live_polynomial_terms_with_literal_reserve']<=150000
  if unit:assert d['all_QQ_denominators_prime_units']and d['literal_modular_Bezout_verified']and d['entire_special_fiber_quotient_dimension']==0
 cert={'whole_P100_original_quartic_pointwise_lower150':unit,'proof_if_unit':['TheNEWNisoneactual9minor: originalsource226C68W, topoldcase9sources, allunitLCDsverified.','ExactN=H*Q+R, withHmonic15andallQQcoefficients32003-integral.','LiteralFp32003Bezout(Hbar,Rbar)=1 meansENTIRESbar/(Rbar)=0,includingnilpotents.','S/(R)isfiniteoverlocalZ_(32003), soNakayama givesS/(R)=0; afterQQ N hasno commonrootwithH.','Theold5actual9minorideal(H), withknowncurveunitnormalizations, plusN coversallT100points, giving141+9=150.'],
       'if_nonunit':'InconclusiveforthisoneNEWsourceborder; noactualrankfailure,newcolumn,source,primeorautomaticchase.',
       'NO_N_LC_condition_needed':True,'no_QQ_Bezout_global9inverse_or_full239_claim':True,'previous_global8_targetsplit_pure90_unchanged':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');audit={k:v for k,v in d.items()if k not in['exact_determinant_N','monic_H15','division_quotient','finite_remainder_R','Bareiss_steps','exact_quotient_DAG','H_reduction','R_reduction','finite_modular_gcd','modular_Bezout_H_cofactor','modular_Bezout_R_cofactor']}
 report={'id':ID,'complete':complete,'NEW_literalC68W_bottom_only':True,'one_actual_det_degree':degN,'whole_P100_original_quartic_pointwise_lower150':unit,'native_audit':audit,'NO_old_K2W_det_rerun_or_new_NF_query':True,'no_auto_nextcolumn_prime_source_or_rankjob':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+sum(p.stat().st_size for p in sb.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 phase.write_text(json.dumps({'status':'NEW_math_complete','complete':complete,'output':str(out.relative_to(ROOT)),'report_sha256':sha(out/'report.json'),'whole150':unit,'no_rerun_authorized':True},indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
