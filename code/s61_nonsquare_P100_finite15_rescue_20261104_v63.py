#!/usr/bin/env python3
"""ONE fixed9det and whole finiteH15 Nakayama gate; no new source query."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_FINITE15_RESCUE_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 sb=base/'S61_NONSQUARE_P100_FINITE15_RESCUE_SOURCE_20261104_v63';sr=read(sb/'report.json');s=read(sb/'source.json');ip=sb/'source-input.txt';assert sr['source_passed']and sha(ip)==sr['native_input_sha256']and sha(sb/'source.json')==sr['source_sha256']
 pb=base/'S61_NONSQUARE_P100_FINITE15_RESCUE_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');assert pr['RF_input_terms']==1043 and pr['det_terms_upper']==215
 hb=base/'S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';hr=read(hb/'report.json');hn=read(hb/'native.json');assert hr['exact_QQ_gcd_degree']==15 and hn['cofactor_QQ_den_and_LC_primeunit_gate']and hn['oldprime_cofactor_gcd_is1']
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp';paths=[Path(__file__),cpp,utility,sb/'report.json',sb/'source.json',ip,pb/'report.json',pb/'certificate.json',hb/'report.json',hb/'native.json']
 binding={'case':s['exact_case'],'finite_H_model':s['finite_H_model'],'actual_source_order':[9,2,3,5,7,6,1,4,11],'target_oldlocal98_order':[0,1,2,3,4,5,6,9,28],'new_original122_target':s['exact_case']['additional_target_original122'],
          'full_curve_poles':'AllRFsourceunitdenominatorsareonlyz,z+20,z+60,203z+2700,z+12. ClearpolynomialrowLCDsBEFOREdet; noneoftheoldDiorfinite-ringmatrixpivotsisinverted.',
          'finite_model_context':'Hmonic15p-integral: S=Z_(32003)[z]/H finitefree15; nofieldorreducednessassumed. TheunitgateusesENTIREspecialfiber, notalocalizedopen.',
          'H_char0_unit_factors':'Hdividesfullyunit-strippedoldprimitivepolynomials,soitschar0rootsexcludethe5knownlinearunitzeros. LCDmultipliersarevalidunitsforactualminorsonT100.',
          'no_new_N_LC_condition':'OnlyHmonic/p-integralanchorrequired; N/Rmusthavep-integralcoefficients. DegreeofRmaydropmodpwithoutinvalidatingfiniteNakayama.'}
 (out/'source-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'finite15-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 d=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode};complete=r.returncode==0 and d['all_passed'];unit=complete and d.get('finite_special_fiber_unit',False)
 degN=max((m[0]for m,c in d['exact_determinant_N']),default=-1)if complete else None
 if complete:
  assert d['all81_source_LCD_identities_verified']and d['literal_N_equals_HQ_plus_R_verified']and degN<=214 and d['peak_live_polynomial_terms_with_literal_reserve']<=150000
  if unit:assert d['all_QQ_denominators_prime_units']and d['literal_modular_Bezout_verified']and d['entire_special_fiber_quotient_dimension']==0
 cert={'whole_P100_quartic_pointwise_lower150':unit,'proof_if_unit':[
       'N isoneactual9minor withknownunitrowLCDs, sourceordersandall81clearingsverified.','N=H*Q+R isliteraloverQQ[z]; allH,Q,Rcoefficientsare32003-integral andHismonic15.',
       'ThepolynomialFpBezoutidentity(Hbar,Rbar)=1 killsENTIRESbar/(Rbar),includingnilpotents.','M=S/(R), S=Z_(32003)[z]/H, isfiniteoverthelocalringZ_(32003). M/pM=0 impliesM=0byNakayama, noflatnessassumptionneeded.',
       'AfterQQbasechange N has nocommonrootwithH. Theoldfive9minorshaveexactgcdH, sooldfiveplusNcoverEVERYP100point. SourceLCDsandQQcontentsarecurveunits, yielding9pure directions+accepted141block.'],
       'if_nonunit_scope':'Onlyinconclusivefinite-special-fiberfrontier; noactualQQrankfailure orautomaticnewcolumn/prime/source.',
       'previous_8target_pure90_constructor_unchanged':True,'no_QQ_Bezout_global9inverse_factorsearch_NF_or_newpoint':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');audit={k:v for k,v in d.items()if k not in['exact_determinant_N','monic_H15','division_quotient','finite_remainder_R','Bareiss_steps','exact_quotient_DAG','H_reduction','R_reduction','finite_modular_gcd','modular_Bezout_H_cofactor','modular_Bezout_R_cofactor']}
 report={'id':ID,'complete':complete,'one_actual_det_degree':degN,'whole_P100_pointwise_quartic_lower150':unit,'native_audit':audit,'no_auto_nextcolumn_prime_source_QQ_Bezout_or_rankjob':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+sum(p.stat().st_size for p in sb.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
