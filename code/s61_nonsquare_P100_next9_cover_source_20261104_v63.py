#!/usr/bin/env python3
"""Freeze fixed5 source/LCD/old-prime gates BEFORE any determinant arithmetic."""
import argparse,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_NEXT9_COVER_SOURCE_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 pb=base/'S61_NONSQUARE_P100_NEXT9_COVER_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');pc=read(pb/'certificate.json');assert sha(pb/'certificate.json')==pr['certificate_sha256']
 ob=base/'S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63';op=read(ob/'operator.json');refs=op['nodes'][op['exports']['M15x98']]['row_refs'];rows=[];paths=[]
 for v in refs:
  p=ROOT/v['path'];assert sha(p)==v['sha256'];rows.append(read(p)['values98']);paths.append(p)
 factors=[(Q(1),Q(0)),(Q(1),Q(20)),(Q(1),Q(60)),(Q(203),Q(2700)),(Q(1),Q(12))];cache={};punit=True
 def dv(poly):
  key=tuple((m[0],c)for m,c in poly)
  if key in cache:return cache[key]
  r={m[0]:Q(c)for m,c in poly};es=[]
  for a,b in factors:
   n=0
   while r and max(r)>=1:
    q={};v=r.copy()
    while v and max(v)>=1:
     j=max(v);c=v[j]/a;q[j-1]=q.get(j-1,Q())+c
     for k,t in[(j,c*a),(j-1,c*b)]:
      v[k]=v.get(k,Q())-t
      if not v[k]:del v[k]
    if v:break
    r=q;n+=1
   es.append(n)
  assert r and max(r)==0;cache[key]=es;return es
 text=['51034','32003','5'];cases=[];inputterms=0;matrixterms=0
 def rf(v):
  nonlocal punit,inputterms
  for k in['numerator','denominator']:
   punit&=all(Q(c).denominator%32003!=0 for m,c in v[k]);inputterms+=len(v[k]);text.append(str(len(v[k])));text.extend(f'{m[0]} {c}'for m,c in v[k])
 for case in pc['five_fixed_cases']:
  rr=case['nine_source_rows_in_order'];cc=case['nine_local98_target_columns'];assert len(rr)==len(set(rr))==9 and cc==pc['first8_target_local98_columns']+[25]
  old=case['nine_LCD_powers5'];new=[];extra=[]
  for row,e in zip(rr,old):
   needed=[max(dv(rows[row][col]['denominator'])[j]for col in cc)for j in range(5)];nn=[max(a,b)for a,b in zip(e,needed)];new.append(nn);extra.append([a-b for a,b in zip(nn,e)])
  text.append(str(case['old_cover_basis_index']));text.extend(' '.join(map(str,e))for e in new)
  for row in rr:
   for col in cc:rf(rows[row][col])
  degs=[]
  for row,e in zip(rr,new):
   ds=[max(m[0]for m,c in rows[row][col]['numerator'])+sum(e)-max(m[0]for m,c in rows[row][col]['denominator'])for col in cc if rows[row][col]['numerator']];degs.append(max(ds,default=0));matrixterms+=sum(d+1 for d in ds)
  cases.append({**case,'effective_LCD_powers5':new,'extra_unit_scalings_over_accepted8cover':extra,'effective_det_degree_upper':sum(degs)})
 assert all(not any(e)for c in cases for e in c['extra_unit_scalings_over_accepted8cover'])
 assert [c['effective_det_degree_upper']for c in cases]==[231,226,214,227,233]and matrixterms==8206 and inputterms==6164
 ip=out/'source-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<2*1024*1024
 source={'curve_ring':'T100','prime_for_coefficient_polynomials_ONLY':32003,'cases5':cases,'input_RF_terms_with_duplicate_uses':inputterms,'cleared_entry_terms_upper':matrixterms,'all405_denominators_divide_known5factor_LCD':True,'extra_scalings_not_needed':True,'all_input_QQ_coefficient_denominators_32003_units':punit,'old_parent_literal_reserve_plus_current_input':51034,'future_native_one_matrix_streamed':True,'no_determinant_or_series_point_rank_gcd_computed':True}
 (out/'source.json').write_text(json.dumps(source,indent=2)+'\n')
 allpaths=[Path(__file__),pb/'report.json',pb/'certificate.json',ob/'report.json',ob/'operator.json']+paths
 report={'id':ID,'source_passed':punit,'source_only':True,'fixed_cases5':True,'QQ_input_prime_unit_gate':punit,'degree_bounds':[231,226,214,227,233],'input_terms':inputterms,'cleared_entry_terms_upper':8206,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in allpaths},'source_sha256':sha(out/'source.json'),'native_input_sha256':sha(ip),'seconds':time.monotonic()-start,'frozen':True,'determinant_math_not_started':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
