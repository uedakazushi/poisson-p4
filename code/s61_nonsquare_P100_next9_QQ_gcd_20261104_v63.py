#!/usr/bin/env python3
"""ONE exact gcd-only on the five saved primitive polynomials; no newminor."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 nb=base/'S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63';nr=read(nb/'report.json');nd=read(nb/'native.json');assert sha(nb/'report.json')=='5670647ac97cb82816852b21b50e70783d17c114e478a7cfdfd95b2553cba63c'and sha(nb/'native.json')==nr['evidence_sha256']['native.json']
 assert nr['complete']and nd['all_QQ_denominators_prime_units']and nd['all_nonzero_primitive_LC_prime_units']and nd['all_reduced_degrees_preserved']
 polys=[d['primitive_ZZ_polynomial']for d in nd['determinants5']];assert sum(map(len,polys))==409 and all('/'not in c for p in polys for m,c in p)
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp';paths=[Path(__file__),cpp,utility,nb/'report.json',nb/'native.json',nb/'source-binding.json',nb/'certificate.json']
 text=[str(10629+409)]
 for p in polys:text.append(str(len(p)));text.extend(f'{m[0]} {c}'for m,c in p)
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n')
 source={'fixed_input':'The5savedprimitiveZZpolynomialsofexact9minors; no source/target/minor recomputation.','cover_basis_order':[d['cover_basis_index']for d in nd['determinants5']],'degrees':[max(m[0]for m,c in p)for p in polys],'integer_coefficient_terms':409,'existing_mod_degree15_upper_not_equality_assumed':True,'whole_curve_units_only':'ActualminorsareQQunitcontent*known5curveunitpowers*thesepolynomials.','no_QQ_Bezout_output':True}
 (out/'source-binding.json').write_text(json.dumps(source,indent=2)+'\n');del nd
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'QQ-gcd-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:run=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 d=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':run.returncode};complete=run.returncode==0 and d['all_passed']
 degree=d.get('gcd_degree')if complete else None
 if complete:assert d['all5_literal_gcd_times_quotient_equals_input_verified']and d['exact_gcd_justified_by_divisibility_and_coprime_cofactors']and(d['oldprime_cofactor_gcd_is1']or d['primitive_PRS_fallback_used'])and d['peak_live_polynomial_terms_with_literal_input']<=150000
 cert={'exact_chosen_minor_ideal_generator_monic_degree':degree,'whole_P100_quartic_lower150':degree==0 if complete else False,
       'exact_gcd_proof':'Hdividesall5literally; cofactorsarecoprimeoverQQ bydegree-preservedoldprimeBezoutorcompleteprimitivePRS. HenceH isexactmonicgcd, notmerelyamodularupper.',
       'finite_cover_if_nonconstant':'Thechosen9minoridealinT100equals(H). Hhasnoknownlinearunitfactor sinceitdividesalreadyfullystrippedinput. QQ[z]/H isafinitealgebra ofdimensiondegreeH, notassumedfield/reduced. ActualquarticrankfailureisNOTinferred.',
       'proper_integral_model_if_gate':'IfHcoefficientsare32003-integral(aschecked), V[z]/H withV=Z_(32003)isfinitefree. A futureextra-polynomialunitproofmustchecktheENTIREspecialfiber, neverjustanopenorapoint.',
       'no_QQ_Bezout_factorsearch_newminor_or_automatic_rescue':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 audit={k:v for k,v in d.items()if k not in['exact_gcd_monic','primitive_ZZ_gcd','all5_quotients','modular_cofactor_gcd','cofactor_modular_Bezout5','cofactor_reductions5','primitive_PRS_transcript']}
 report={'id':ID,'complete':complete,'exact_QQ_gcd_degree':degree,'whole_P100_pointwise_lower150':degree==0 if complete else False,'native_audit':audit,'no_QQ_Bezout_factorsearch_newminor_point_or_prime':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
