#!/usr/bin/env python3
"""ONE fixed5 Bareiss9/gcd-mod32003 gate from frozen source preflight."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 sb=base/'S61_NONSQUARE_P100_NEXT9_COVER_SOURCE_20261104_v63';sr=read(sb/'report.json');s=read(sb/'source.json');ip=sb/'source-input.txt'
 assert sr['source_passed']and sr['QQ_input_prime_unit_gate']and sha(ip)==sr['native_input_sha256']and sha(sb/'source.json')==sr['source_sha256']
 pb=base/'S61_NONSQUARE_P100_NEXT9_COVER_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');assert pr['five_degree_bounds']==[231,226,214,227,233]
 coverpeer=base/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63_REVIEW/report.json';assert sha(coverpeer)=='0e830d25d8631a7e38da157d15e7be181dac4dc80bc6c282381cc9e64f901f1a'
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp';paths=[Path(__file__),cpp,utility,sb/'report.json',sb/'source.json',ip,pb/'report.json',pb/'certificate.json',coverpeer]
 labels={'new_source':'CorrectedpurekernelK2(actualnonpivot23)*physicalv3; Mrow11, not uncorrectedC23.',
         'ninth_oldlocal98_target':25,'actual_residual90_column':17,'ninth_original122_target':30,'ninth_parent239_target':115,
         'five_cases':s['cases5'],'same8target_cover_preserved':'Allold/newLCDscovertheninthcolumn; extraunitrowpowersallzero. Source405entryidentitiesverifiednative.',
         'determinant_algorithm':'Fixed9x9row-pivotBareiss, exactpolynomialdivisions. Savedpivots/row-swaps/exact-quotientDAGpermitdeterministicreplay; no ranksearch across source/targetfamilies.',
         'primitive_definition':'Afteronly5knownlinearunitfactors, clearQQcoefficientdenominators anddivideintegercontent; choosepositiveleadingprimitiveZZpolynomial.',
         'modular_gate':'Existingprime32003appliedtotheWHOLEcoefficientpolynomial, no numericzetapoint. RequireALLQQdenominatorsp-units,primitiveleadingcoefficientnonzeroanddegreepreserved. Otherwiseinconclusive.',
         'no_expanded_QQ_Bezout_or_nine_target_inverse':True}
 (out/'source-binding.json').write_text(json.dumps(labels,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'next9-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 nd=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode};complete=r.returncode==0 and nd['all_passed'];unit=complete and nd.get('whole_curve_quartic_lower150',False)
 degrees=[]
 if complete:
  assert nd['all405_source_LCD_identities_verified']and nd['peak_live_polynomial_terms_with_literal_reserve']<=150000
  degrees=[max((m[0]for m,c in d['exact_determinant']),default=-1)for d in nd['determinants5']]
  assert all(d<=b for d,b in zip(degrees,[231,226,214,227,233]))
  if unit:assert nd['all_QQ_denominators_prime_units']and nd['all_nonzero_primitive_LC_prime_units']and nd['all_reduced_degrees_preserved']and nd['literal_modular_Bezout_verified']and nd['modular_gcd_is1']
 cert={'pointwise_original_quartic_lower150_on_whole_P100':unit,'Gauss_bridge_if_unit':[
       'Eachnonzero primitiveZi liesinZZ[z]withleadingcoefficienta32003-unit; allcoefficientdenominatorgatespass.',
       'AnycommonmonicfactoroverQQwould, byGauss/monicfactorintegralityoverZ_(32003), havep-integralcoefficientsandretainpositive degree onreduction.',
       'TheliteralpolynomialFpBezoutsum bi*Zi_bar=1 contradictssuchacommonfactor. Thusgcd_QQ(Zi)=1, withno irreducibilityassumption.',
       'TheactualminorsdifferfromZi onlybynonzeroQQcontentandfivecurveunitfactors, sohaveNOcommoncomplexzeroonT100.',
       'Oneactual9minor isnonzero ateachpoint; withtheaccepted141blockthisgivesoriginalweightedquarticrank>=150.'],
       'inconclusive_scope':'FailedLC/p-den/unitgcdgate onlyinconclusive; thechosenfamilymayleavefiniteorwholeunresolvedset. No rankfailure, newprime, extra source/column or two-rowminor follows.',
       'not_constructed':['QQpolynomialBezout','global9targetrightinverse','pure89constructor','anynewgenericclassification'],'previous8targetsplit_remains_valid':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 audit={k:v for k,v in nd.items()if k not in['determinants5','modular_gcd','modular_Bezout5','all5_modular_primitive_polynomials']}
 report={'id':ID,'complete':complete,'five_actual_det_degrees':degrees,'whole_P100_quartic_pointwise_lower150':unit,'native_audit':audit,'no_QQ_Bezout_ranksearch_NF_query_point_or_newprime':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+sum(p.stat().st_size for p in sb.iterdir()if p.is_file())<32*1024*1024
 signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
