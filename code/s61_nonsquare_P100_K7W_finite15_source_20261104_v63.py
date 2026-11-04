#!/usr/bin/env python3
"""Freeze the NEW K7W fixedborder input, not the previous K2W border."""
import argparse,json,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_K7W_FINITE15_SOURCE_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports';pb=base/'S61_NONSQUARE_P100_K7W_FINITE15_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');pc=read(pb/'certificate.json');assert sha(pb/'certificate.json')==pr['certificate_sha256']and pr['support_matching']
 ob=base/'S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63';op=read(ob/'operator.json');refs=op['nodes'][op['exports']['M15x98']]['row_refs'];paths=[];rows=[]
 for i in pc['old_source_rows_in_order8']:
  r=refs[i];p=ROOT/r['path'];assert sha(p)==r['sha256'];rows.append(read(p)['values98']);paths.append(p)
 vp=ROOT/pc['new_bottom_source']['file'];assert sha(vp)==pc['new_bottom_source']['sha256'];rows.append(read(vp)['values98']);assert pc['new_bottom_source']['original319_source_index']==226
 hb=base/'S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';hr=read(hb/'report.json');hn=read(hb/'native.json');assert hr['exact_QQ_gcd_degree']==15 and sha(hb/'native.json')==hr['evidence_sha256']['native.json'];H=hn['exact_gcd_monic'];assert H[-1]==[[15],'1']
 text=['51034','32003'];punit=all(Q(c).denominator%32003!=0 for m,c in H);terms=0
 def poly(p):text.append(str(len(p)));text.extend(f'{m[0]} {c}'for m,c in p)
 poly(H);text.append('9');text.extend(' '.join(map(str,e))for e in pc['nine_LCD_powers5'])
 for row in rows:
  for c in pc['nine_oldlocal98_target_columns']:
   v=row[c]
   for k in['numerator','denominator']:punit&=all(Q(q).denominator%32003!=0 for m,q in v[k]);terms+=len(v[k]);poly(v[k])
 assert terms==1057;ip=out/'source-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<1024*1024
 source={'case9_top8_source_ids':pc['old_source_rows_in_order8'],'NEW_bottom_source':pc['new_bottom_source'],'nine_LCD_powers5':pc['nine_LCD_powers5'],'nine_oldlocal98_targets':pc['nine_oldlocal98_target_columns'],'original122_targets':pc['original122_target_columns'],'monicH15_ref':pc['finite_H15_external'],'all_selected_QQ_den_punit':punit,'same8cover_weights_retained':True,'NO_old_failed_K2W_bottom_source':True,'math_not_started':True}
 (out/'source.json').write_text(json.dumps(source,indent=2)+'\n');allpaths=[Path(__file__),pb/'report.json',pb/'certificate.json',ob/'operator.json',hb/'report.json',hb/'native.json',vp]+paths
 report={'id':ID,'source_passed':punit,'source_only':True,'NEW_source319row226':True,'RF_input_terms':terms,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in allpaths},'source_sha256':sha(out/'source.json'),'native_input_sha256':sha(ip),'seconds':time.monotonic()-start,'frozen':True,'math_not_started':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
