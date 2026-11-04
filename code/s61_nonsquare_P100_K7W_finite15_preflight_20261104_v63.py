#!/usr/bin/env python3
"""Source-only SAME finiteH15 border with NEW K7W bottom; no determinant."""
import argparse,json,math,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_K7W_FINITE15_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports';nb=base/'S61_NONSQUARE_P100_K7W_PURE98_20261104_v63';nr=read(nb/'report.json');vp=nb/'pure98-row3.json';v=read(vp);assert nr['complete']and sha(vp)==nr['evidence_sha256'][vp.name]
 ob=base/'S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63';op=read(ob/'operator.json');case=next(c for c in op['five_minors']if c['basis_index']==9);refs=op['nodes'][op['exports']['M15x98']]['row_refs'];top=case['source_rows_in_order'];rows=[];paths=[]
 for i in top:
  p=ROOT/refs[i]['path'];assert sha(p)==refs[i]['sha256'];d=read(p);assert d['target122_columns']==v['target122_columns'];rows.append(d['values98']);paths.append(p)
 rows.append(v['values98']);C=op['context']['selected_local98_columns8'];cols=C+[28];assert top==[9,2,3,5,7,6,1,4]
 factors=[(1,0),(1,20),(1,60),(203,2700),(1,12)];cache={}
 def denvec(p):
  key=tuple((m[0],c)for m,c in p)
  if key in cache:return cache[key]
  r={m[0]:Q(c)for m,c in p};es=[]
  for a,b in factors:
   n=0
   while r and max(r)>=1:
    t=r.copy();q={}
    while t and max(t)>=1:
     j=max(t);c=t[j]/a;q[j-1]=q.get(j-1,Q())+c
     for k,x in[(j,c*a),(j-1,c*b)]:
      t[k]=t.get(k,Q())-x
      if not t[k]:del t[k]
    if t:break
    r=q;n+=1
   es.append(n)
  assert r and max(r)==0;cache[key]=es;return es
 ex=[max(denvec(x['denominator'])[j]for x in rows[-1])for j in range(5)];es=case['LCDpowers_in_order']+[ex];deg=lambda p:max((m[0]for m,c in p),default=-1)
 def height(p):
  qq=[Q(c)for m,c in p];d=1
  for q in qq:d=math.lcm(d,q.denominator)
  return d,max([abs(q.numerator*(d//q.denominator)).bit_length()for q in qq]+[0])
 score=0;entryterms=0;degrees=[];heights=[];dens=[];punit=True
 for row,e in zip(rows,es):
  vals=[row[c]for c in cols];dd=[height(x['numerator'])[0]for x in vals];common=math.lcm(*dd);qh=sum(k*((a+b).bit_length()+a.bit_length())for k,(a,b)in zip(e,factors));hh=[];ds=[]
  for x,d in zip(vals,dd):
   assert all(a<=b for a,b in zip(denvec(x['denominator']),e));score+=len(x['numerator'])+len(x['denominator']);punit&=all(Q(c).denominator%32003!=0 for part in['numerator','denominator']for m,c in x[part])
   if x['numerator']:
    k=deg(x['numerator'])+sum(e)-deg(x['denominator']);ds.append(k);entryterms+=k+1;hh.append(height(x['numerator'])[1]+qh+(common//d).bit_length())
  degrees.append(max(ds,default=0));heights.append(max(hh,default=0));dens.append(common)
 pattern=[[bool(row[c]['numerator'])for c in cols]for row in rows];used={}
 def dfs(r,seen):
  for c in range(9):
   if pattern[r][c]and c not in seen:
    seen.add(c)
    if c not in used or dfs(used[c],seen):used[c]=r;return True
  return False
 matching=all(dfs(r,set())for r in range(9));db=sum(degrees);nbits=sum(heights)+sum((x+1).bit_length()for x in degrees)+math.factorial(9).bit_length();dbytes=(db+1)*(32+math.ceil((nbits+math.prod(dens).bit_length())*math.log10(2)))
 hb=base/'S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';hn=read(hb/'native.json');assert hn['gcd_degree']==15 and hn['cofactor_QQ_den_and_LC_primeunit_gate']
 cert={'old8cover_basis_index':9,'old_source_rows_in_order8':top,'new_bottom_source':{'file':str(vp.relative_to(ROOT)),'sha256':sha(vp),'original319_source_index':226,'literal_normalizedC68_times_physicalv3':True},'nine_LCD_powers5':es,'nine_oldlocal98_target_columns':cols,'original122_target_columns':[v['target122_columns'][c]for c in cols],
       'structural_perfect_matching_support':matching,'RF_input_terms':score,'cleared_entry_poly_terms_upper':entryterms,'det_degree_upper':db,'det_terms_upper':db+1,'raw_det_JSON_bytes_upper':dbytes,'all_source_QQ_den_punit':punit,
       'finite_H15_external':{'file':str((hb/'native.json').relative_to(ROOT)),'sha256':sha(hb/'native.json'),'pointer':'/exact_gcd_monic'},
       'future_recipe':'ONEnew9detQQbyglobalpolynomialBareiss; literalremmonicH15; wholeFp32003gcd/Bezout(Hbar,Rbar); finiteV-Nakayamaifunit. NO Diinverse, nofinite-ringmatrixpivot, noNLCcondition.',
       'source_scope':'UsesNEWliteralC68W row, notoldfailedK2Wminor. Accepted8coverweightsandcolumnsareunchanged. AllfullpuretailsremainfromactualNF117source.',
       'if_unit':'Old5minorideal(H)plusnewminorunitcoversT100=>originalquartic>=150; noQQ9inverseconstructed.','if_nonunit':'Onlyinconclusiveforfixedsource/column; noactualrankfailure orautomaticothercase/column/prime.',
       'no_determinant_remainder_gcd_or_rank_computed':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');allpaths=[Path(__file__),nb/'report.json',vp,nb/'source-binding.json',ob/'operator.json',hb/'report.json',hb/'native.json']+paths
 report={'id':ID,'source_only':True,'same_case9_and_column28':True,'new_source319row226':True,'support_matching':matching,'RF_terms':score,'cleared_entry_terms_upper':entryterms,'det_degree_upper':db,'det_terms_upper':db+1,'raw_det_JSON_bytes_upper':dbytes,'input_punit':punit,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in allpaths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'no_new9minor_or_finite_unit_job':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
