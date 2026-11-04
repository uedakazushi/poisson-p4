#!/usr/bin/env python3
"""Source-only five9minor cover forecast. No determinant, gcd or residual value."""
import argparse,json,math,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_NEXT9_COVER_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports';ob=base/'S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63';op=read(ob/'operator.json');orr=read(ob/'report.json');assert sha(ob/'operator.json')==orr['operator_sha256']
 refs=op['nodes'][op['exports']['M15x98']]['row_refs'];paths=[ROOT/v['path']for v in refs];M=[]
 for p,v in zip(paths,refs):assert sha(p)==v['sha256'];M.append(read(p)['values98'])
 cb=base/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63/source-binding.json';cover=read(cb);C=op['context']['selected_local98_columns8'];T=op['context']['remaining_local98_columns90'];cases=op['five_minors'];w=11;we=cover['new_sources_in_order'][3]['known_unit_rowLCD_powers5'];assert cover['new_sources_in_order'][3]['physical_light']==3
 def match(mat):
  used={}
  def dfs(r,seen):
   for c in range(9):
    if mat[r][c]and c not in seen:
     seen.add(c)
     if c not in used or dfs(used[c],seen):used[c]=r;return True
   return False
  return all(dfs(r,set())for r in range(9))
 def degree(p):return max((m[0]for m,c in p),default=-1)
 candidates=[]
 for col in T:
  score=0;degrees=[];ok=True
  for ca in cases:
   rows=ca['source_rows_in_order']+[w];es=ca['LCDpowers_in_order']+[we];cols=C+[col]
   assert len(set(rows))==9
   if not match([[bool(M[r][c]['numerator'])for c in cols]for r in rows]):ok=False;break
   score+=sum(len(M[r][c]['numerator'])+len(M[r][c]['denominator'])for r in rows for c in cols)
   degrees.append(sum(max((degree(M[r][c]['numerator'])+sum(e)-degree(M[r][c]['denominator'])for c in cols if M[r][c]['numerator']),default=0)for r,e in zip(rows,es)))
  if ok:candidates.append((score,col,degrees))
 assert candidates;score,col,degree_bounds=min(candidates);assert col==25
 # Integer coefficient-height bounds only, not QQ[z] products or determinants.
 def height(p):
  qq=[Q(c)for m,c in p];d=1
  for q in qq:d=math.lcm(d,q.denominator)
  return d,max([abs(q.numerator*(d//q.denominator)).bit_length()for q in qq]+[0])
 factors=[(1,0),(1,20),(1,60),(203,2700),(1,12)];profiles=[];totalbytes=0;totalentries=0
 for ca,db in zip(cases,degree_bounds):
  rows=ca['source_rows_in_order']+[w];es=ca['LCDpowers_in_order']+[we];cols=C+[col];rowheights=[];rowdens=[];rowdegrees=[];entryterms=0
  for r,e in zip(rows,es):
   vals=[M[r][c]for c in cols];dens=[height(v['numerator'])[0]for v in vals];common=1
   for d in dens:common=math.lcm(common,d)
   qbound=sum(ex*((abs(a)+abs(b)).bit_length()+a.bit_length())for ex,(a,b)in zip(e,factors));hh=[];degrees=[]
   for v,d in zip(vals,dens):
    if not v['numerator']:continue
    assert Q(v['denominator'][-1][1])==1
    nd=degree(v['numerator'])+sum(e)-degree(v['denominator']);degrees.append(nd);entryterms+=nd+1
    hh.append(height(v['numerator'])[1]+(common//d).bit_length()+qbound)
   rowheights.append(max(hh,default=0));rowdens.append(common);rowdegrees.append(max(degrees,default=0))
  den=math.prod(rowdens);nb=sum(rowheights)+sum((d+1).bit_length()for d in rowdegrees)+math.factorial(9).bit_length();byte=(db+1)*(32+math.ceil((nb+den.bit_length())*math.log10(2)));totalbytes+=byte;totalentries+=entryterms
  profiles.append({'old_cover_basis_index':ca['basis_index'],'nine_source_rows_in_order':rows,'nine_LCD_powers5':es,'nine_local98_target_columns':cols,'determinant_degree_upper':db,'determinant_terms_upper':db+1,'all81_cleared_entry_terms_upper':entryterms,'det_common_QQ_den_bit_upper':den.bit_length(),'det_integer_num_bit_upper':nb,'raw_det_JSON_bytes_upper':byte})
 hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';h=read(hp);lab122=op['context']['pure98_target122_labels'][col];lab239=h['parent_target122_monomial_indices'][lab122]
 cert={'fixed_new_source':{'M_row':w,'kernel':2,'actual_nonpivot135':23,'physical_light':3,'whole_curve_kernel_source':True},'selected_ninth_local98_target':col,'original122_target':lab122,'parent239_target':lab239,'joint_weight4_monomial':h['parent_joint_target239_monomials'][lab239],
       'first8_target_local98_columns':C,'support_only_column_scan':{'remaining_columns':90,'structurally_admissible_for_all5_cases':len(candidates),'least_RF_input_score':score,'no_polynomial_value_or_rank_used':True},'five_fixed_cases':profiles,
       'all5_det_terms_upper':sum(d+1 for d in degree_bounds),'all5_cleared_matrix_entries_terms_upper':totalentries,'all5_raw_det_JSON_bytes_upper':totalbytes,
       'literal_minor_definition':'For each acceptedweighted8selectorP_i use rowsP_i plus LCDw*e_w; columnsC+[25]. Source order8old thennewW; no rowpermutation. These are actual9minors withknownunitrowLCDs.',
       'positive_uniform_gate':'Exactgcd ofthese5 fixed9minors knownunit afteronly5curvelocalizers =>wholeP100originalweightedquartic>=150. Ifnonunit, retainfiniteunresolvedcover; no actualrankfailure orautomaticnewcolumn/source.',
       'global_residual_relation':'WithN_i thedefined9minor, E*LCDw*(R_w,25)=sum_i Ai*N_i. Prefer5smallminorpolynomials toexpandingAi-weightedglobalRentry.',
       'future_algorithm_bound':'Streamone9x9 polynomialmatrix atatime. Fixedexactdeterminant/source-LCDidentities andliteralgcdcertificate only; no full15x98Gaussian. Determinantandgcdoutputstillneedlive150k/32MiBguards.',
       'no_minor_gcd_residual90_value_or_rank_computed':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 allpaths=[Path(__file__),ob/'report.json',ob/'operator.json',cb,hp]+paths
 report={'id':ID,'source_only':True,'fixed_five9minor_candidate':True,'ninth_local98_target':col,'five_degree_bounds':degree_bounds,'all5_det_terms_upper':cert['all5_det_terms_upper'],'all5_raw_det_JSON_bytes_upper':totalbytes,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in allpaths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'no_determinant_gcd_or_rank_job':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
