#!/usr/bin/env python3
"""Metadata-only explicit reading of prior literal-copy coefficient bound."""
from pathlib import Path
from fractions import Fraction as Q
import json
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_HEIGHT_20261105_v63'
def main():
 base='reports/AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_SOURCE_20261105_v63/';p=ROOT/base/'report.json';r=json.loads(p.read_text());assert r['all_passed'] and sha(p)=='d490c810468a9b5fc3011e762406e7342fbc8dbdc1978f9837e235664eb00c3b'
 sp='reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json';B=fields(ROOT/sp,['YB_in_YA_linear'])['YB_in_YA_linear'];height=max(max(abs(Q(x).numerator).bit_length(),Q(x).denominator.bit_length()) for row in B[:7] for x in row[:7]);assert height==7
 priorphase=json.loads((ROOT/base/'phase-positive-before-source-checks.json').read_text());assert 'copiesofliteralBlinearcoefficients' in priorphase['positive_literal_product_bound']
 first='reports/AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_SOURCE_20261105_v63_FIRST_MONOMIAL_READER/';old=ROOT/first/'failed-source-reader.py';failure=ROOT/first/'failure.json';f=json.loads(failure.read_text());assert sha(old)==f['failed_code_sha256'] and not f['native_or_functionresponse_values_evaluated']
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 report={'id':STEM,'all_passed':True,'metadata_only_no_product_or_NF_rerun':True,'primary_unchanged_sha256':sha(p),'primary_positive_prior1_scope':'lambda10timesNF660dotsonly','prior_literal_copy_bound_expanded':'CAlinearformsareunitcoordinateorzero,soeach64productcoefficientiscopyofBlinearcoefficientornothing. Thepriorphasealreadyrecordedthisboundbeforeproducts.','literal_product_copy_numerator_denominator_maxbits':height,'all_new_source_check_arithmetic_prior_bits':max(height,r['positive_prior_bits']),'source_owned_prior_unchanged':r['source_owned_prior'],'whole_RAW_prior_unchanged':r['whole_RAW_prior'],'no_functionresponse_or_new_budget':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p),sp:sha(ROOT/sp),str((ROOT/base/'phase-positive-before-source-checks.json').relative_to(ROOT)):sha(ROOT/base/'phase-positive-before-source-checks.json'),str(old.relative_to(ROOT)):sha(old),str(failure.relative_to(ROOT)):sha(failure)},'code_sha256':sha(Path(__file__))};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'sourceprior':height}))
if __name__=='__main__':main()
