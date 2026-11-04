#!/usr/bin/env python3
"""Only source/degree bounds for the saved8x98 pure block; NO field rank."""
import argparse,json,time
from pathlib import Path
from fractions import Fraction as Q
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE8_RANK_SOURCE_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 split=base/'S61_NONSQUARE_P100_PURE97_SPLIT_20261104_v63';pre=read(split/'preflight.json');sr=read(split/'report.json');assert sr['complete']
 paths=[Path(__file__),split/'report.json',split/'preflight.json',base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/pure98.json']+[base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63'/f'pure98-row{i}.json'for i in range(7)]
 rows=[read(p)['values98']for p in paths[3:]];assert len(rows)==8 and all(len(r)==98 for r in rows)
 for p in paths[3:]:assert sha(p)==sr['input_sha256'][str(p.relative_to(ROOT))]
 deg=lambda p:max((m[0]for m,c in p),default=-1)
 lcd=[pre['q5_common_denominator_degree']]+[p['source_common_denominator_degree']for p in pre['seven_profiles']]
 profiles=[]
 for row,d in zip(rows,lcd):
  nd=max([deg(v['numerator'])for v in row]);dd=max(deg(v['denominator'])for v in row)
  # Multiplication by the existing five-factor row LCD is a unit scaling.
  clear_degree=max([deg(v['numerator'])+d-deg(v['denominator'])for v in row if v['numerator']])
  profiles.append({'nonzero_entries':sum(bool(v['numerator'])for v in row),'RF_terms':sum(len(v['numerator'])+len(v['denominator'])for v in row),'row_LCD_known_curve_unit_degree':d,'cleared_polynomial_max_degree_upper':clear_degree,'cleared_polynomial_terms_upper':98*(clear_degree+1)})
 maxbits=max(max(abs(Q(c).numerator).bit_length(),Q(c).denominator.bit_length())for row in rows for v in row for k in['numerator','denominator']for m,c in v[k])
 cert={'matrix_shape':[8,98],'source_order':['accepted_q5H_after_H1_clear']+[f'G18*v{j}'for j in[0,1,2,3,4,6,7]],'row_profiles':profiles,'RF_input_terms':sum(p['RF_terms']for p in profiles),'max_input_coefficient_bits':maxbits,'cleared_polynomial_terms_upper':sum(p['cleared_polynomial_terms_upper']for p in profiles),'generic_rank_bridge':'rank_QQ(z)(original8x98)=1+rank_QQ(z)(saved7x97) by the accepted explicit invertible target/source transport. No rank value is assumed.','possible_next_algorithm':'ONE bounded exact QQ(z) row-rank/selected-minor transcript on this compact8x98, after known-unit rowLCD clearing; arbitrary polynomial pivots are allowed ONLY for generic rank, never a whole-curve unit claim. Not executed here.','whole_curve_gate':'Generic full row rank would leave finite unknown zeros of selected determinant. A Smith/determinantal-ideal unit proof is separate.','rank_or_minor_or_Gaussian_executed':False,'no_new_source_value_or_other_Gi':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 report={'id':ID,'source_only':True,'shape':[8,98],'input_terms':cert['RF_input_terms'],'max_input_coefficient_bits':maxbits,'known_unit_row_clear_polynomial_terms_upper':cert['cleared_polynomial_terms_upper'],'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'no_fieldjob_launched':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
