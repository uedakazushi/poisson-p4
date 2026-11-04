#!/usr/bin/env python3
"""Freeze ONE newborder source/H15 input before determinant arithmetic."""
import argparse,json,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_FINITE15_RESCUE_SOURCE_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports';pb=base/'S61_NONSQUARE_P100_FINITE15_RESCUE_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');pc=read(pb/'certificate.json');assert sha(pb/'certificate.json')==pr['certificate_sha256']
 hb=base/'S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';hr=read(hb/'report.json');hn=read(hb/'native.json');assert hr['exact_QQ_gcd_degree']==15 and sha(hb/'native.json')==hr['evidence_sha256']['native.json']
 ob=base/'S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63';op=read(ob/'operator.json');refs=op['nodes'][op['exports']['M15x98']]['row_refs'];M=[];paths=[]
 for r in refs:
  p=ROOT/r['path'];assert sha(p)==r['sha256'];M.append(read(p)['values98']);paths.append(p)
 ca=pc['single_fixed_border'];rows=ca['nine_source_rows_in_order'];cols=ca['nine_oldlocal98_target_columns'];ex=ca['known_unit_LCD_powers5'];assert rows==[9,2,3,5,7,6,1,4,11]and cols[-1]==28
 H=hn['exact_gcd_monic'];assert H[-1]==[[15],'1'];text=['51034','32003'];punit=all(Q(c).denominator%32003!=0 for m,c in H)
 def poly(p):text.append(str(len(p)));text.extend(f'{m[0]} {c}'for m,c in p)
 def rf(v):
  nonlocal punit
  for k in['numerator','denominator']:punit&=all(Q(c).denominator%32003!=0 for m,c in v[k]);poly(v[k])
 poly(H);text.append('9');text.extend(' '.join(map(str,e))for e in ex)
 for r in rows:
  for c in cols:rf(M[r][c])
 ip=out/'source-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<1024*1024
 source={'exact_case':ca,'finite_H_model':pc['finite_algebra'],'all_source_QQ_denominators_prime_units':punit,'known_LCD_source_proof':'Prescribedwhole-rowfive-unitLCDprofilesfromacceptedold8/K2coefficientcachescoverALL98targets; nativewillverifyeach81polynomialdivisionliterally. No unknownDiinverse.',
         'source_RF_terms':1043,'cleared_entry_terms_upper':1500,'det_degree_upper':214,'det_terms_upper':215,'ring_operations':'GlobalQQ[z]BareissN, thenmonicH-polynomialremainder; NO finite-ringmatrixpivotdivision.','source_frozen_BEFORE_determinant':True}
 (out/'source.json').write_text(json.dumps(source,indent=2)+'\n');allpaths=[Path(__file__),pb/'report.json',pb/'certificate.json',hb/'report.json',hb/'native.json',ob/'operator.json']+paths
 report={'id':ID,'source_passed':punit,'source_only':True,'all_source_QQ_punit':punit,'H_monic_degree15':True,'source_RF_terms1043':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in allpaths},'source_sha256':sha(out/'source.json'),'native_input_sha256':sha(ip),'seconds':time.monotonic()-start,'frozen':True,'math_not_started':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
