#!/usr/bin/env python3
"""Source-only one additional9border over the exact finiteH15 frontier."""
import argparse,json,math,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_FINITE15_RESCUE_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports';hb=base/'S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63';hr=read(hb/'report.json');hn=read(hb/'native.json');assert hr['complete']and hr['exact_QQ_gcd_degree']==15 and hn['cofactor_QQ_den_and_LC_primeunit_gate']and hn['oldprime_cofactor_gcd_is1']
 ob=base/'S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63';op=read(ob/'operator.json');M=[];rowpaths=[]
 for r in op['nodes'][op['exports']['M15x98']]['row_refs']:
  p=ROOT/r['path'];assert sha(p)==r['sha256'];M.append(read(p)['values98']);rowpaths.append(p)
 coverp=base/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63/source-binding.json';cv=read(coverp);w=11;we=cv['new_sources_in_order'][3]['known_unit_rowLCD_powers5'];C=op['context']['selected_local98_columns8'];T=op['context']['remaining_local98_columns90']
 def match(mat):
  used={}
  def dfs(r,seen):
   for c in range(9):
    if mat[r][c]and c not in seen:
     seen.add(c)
     if c not in used or dfs(used[c],seen):used[c]=r;return True
   return False
  return all(dfs(r,set())for r in range(9))
 degree=lambda p:max((m[0]for m,c in p),default=-1)
 candidates=[]
 for case in op['five_minors']:
  rows=case['source_rows_in_order']+[w];es=case['LCDpowers_in_order']+[we]
  for col in T:
   if col==25:continue
   cols=C+[col]
   if not match([[bool(M[r][c]['numerator'])for c in cols]for r in rows]):continue
   terms=sum(len(M[r][c]['numerator'])+len(M[r][c]['denominator'])for r in rows for c in cols)
   degs=[max((degree(M[r][c]['numerator'])+sum(e)-degree(M[r][c]['denominator'])for c in cols if M[r][c]['numerator']),default=0)for r,e in zip(rows,es)]
   candidates.append((terms,sum(degs),case['basis_index'],col,rows,es))
 score,detdegree,basis,col,rows,es=min(candidates);assert(score,detdegree,basis,col)==(1043,214,9,28);cols=C+[col]
 # Height and regularity forecasts only; no global polynomial or finite residue.
 def height(p):
  qq=[Q(c)for m,c in p];d=1
  for q in qq:d=math.lcm(d,q.denominator)
  return d,max([abs(q.numerator*(d//q.denominator)).bit_length()for q in qq]+[0])
 factors=[(1,0),(1,20),(1,60),(203,2700),(1,12)];rowheight=[];rowden=[];rowdegrees=[];entries=0;punit=True
 for r,e in zip(rows,es):
  vv=[M[r][c]for c in cols];ds=[height(v['numerator'])[0]for v in vv];common=math.lcm(*ds);qb=sum(x*((a+b).bit_length()+a.bit_length())for x,(a,b)in zip(e,factors));hh=[];degs=[]
  for v,d in zip(vv,ds):
   punit&=all(Q(c).denominator%32003!=0 for k in['numerator','denominator']for m,c in v[k])
   if not v['numerator']:continue
   nd=degree(v['numerator'])+sum(e)-degree(v['denominator']);entries+=nd+1;degs.append(nd);hh.append(height(v['numerator'])[1]+qb+(common//d).bit_length())
  rowheight.append(max(hh,default=0));rowden.append(common);rowdegrees.append(max(degs,default=0))
 nbits=sum(rowheight)+sum((d+1).bit_length()for d in rowdegrees)+math.factorial(9).bit_length();dbits=math.prod(rowden).bit_length();byte=(detdegree+1)*(32+math.ceil((nbits+dbits)*math.log10(2)))
 H=hn['exact_gcd_monic'];assert degree(H)==15 and H[-1][1]=='1';h_punit=all(Q(c).denominator%32003!=0 for m,c in H)
 hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';h=read(hp);t122=op['context']['pure98_target122_labels'][col];t239=h['parent_target122_monomial_indices'][t122]
 cert={'finite_algebra':{'QQ':'QQ[z]/(Hmonic15)','integral':'V[z]/H,V=Z_(32003)','dimension':15,'monic_H_external_ref':{'path':str((hb/'native.json').relative_to(ROOT)),'sha256':sha(hb/'native.json'),'json_pointer':'/exact_gcd_monic'},'H_coefficients_32003_integral':h_punit,'not_assumed_field_or_reduced':True},
       'single_fixed_border':{'old8cover_basis_index':basis,'nine_source_rows_in_order':rows,'known_unit_LCD_powers5':es,'nine_oldlocal98_target_columns':cols,'additional_target_local98':col,'additional_target_original122':t122,'additional_parent239':t239,'additional_joint_monomial':h['parent_joint_target239_monomials'][t239]},
       'cost':{'support_scan_only_no_minor_values':True,'RF_input_terms':score,'cleared_entry_polynomial_terms_upper':entries,'determinant_degree_upper':detdegree,'determinant_terms_upper':detdegree+1,'raw_det_JSON_bytes_upper':byte,'all_selected_input_QQ_denominators_32003_units':punit,'future_polynomial_remainder_degree_bound':14},
       'safe_next_recipe':['ClearrowLCDsusingonlythe5knowncurvefactors. ComputeONEfixed9detN globallyoverQQ[z]withBareissexactdivision, notinfinitequotientring.','ReduceNmodthemonicH bypolynomialdivision. ThisinvertsnoDi,nofinite-ringmatrixpivot,nounknownfactor.','CheckallnewQQdenominators32003-units; thencomputeoldprimegcd/BezoutofHbarandtheENTIREresidueNbar.','Ifunit, M=V[z]/(H,N) isfiniteoverVandM/pM=0,soNakayamaM=0. Char0chosenexceptiondisappears; withold5minorcoveroffH thisproveswholeP100quartic>=150.'],
       'important_conditions':['NoassumptionthatanyoldDiisunitmodH. Thefullknown8coverisretainedandonlyitsweightedselectorisused.','DoNOTdiscardanyspecialfiberpointwhereacurveunitreduceszero; positiveunitgatemustholdontheENTIREfiniteSbar,includingnilpotents.','NoLCconditiononnewNisneededforthefiniteNakayamagate; monicHanchorsfiniteness. Allcoefficientsmustremainp-integral.','Hhasnoexcludedchar0linearfactorbecauseitdividesfullyunit-strippedoldprimitiveinputs; originalrowLCDsarethereforeunitsonS_Q.'],
       'if_nonunit':'Inconclusive; no actualrankfailure, alternativeprime, nextcolumn/sourceoradditionalborder isautomatic.',
       'no_determinant_remainder_finite_unit_rank_or_source_query_executed':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),hb/'report.json',hb/'native.json',ob/'report.json',ob/'operator.json',coverp,hp]+rowpaths
 report={'id':ID,'source_only':True,'one_fixed9_border':True,'basis_index':basis,'additional_target28':True,'finite_dimension15':True,'RF_input_terms':score,'cleared_entry_terms_upper':entries,'det_terms_upper':215,'raw_det_JSON_bytes_upper':byte,'input_QQ_punit':punit,'H_QQ_punit':h_punit,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'no_determinant_or_finite_residue_job':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
