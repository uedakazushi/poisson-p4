#!/usr/bin/env python3
"""Postverify114 provided own-source sixth coefficients; never solve them."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import numeric_eval,series_mul
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63_REVIEW'
A='AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63'
PATHS={
 'report':f'reports/{A}/report.json','CA':f'reports/{A}/CA-source6-and-regular-quotient-third.json','CB':f'reports/{A}/CB-source6-and-regular-quotient-third.json',
 'CAcheck':f'reports/{A}/CA-new-sixth-coefficient-checkpoint.json','CBcheck':f'reports/{A}/CB-new-sixth-coefficient-checkpoint.json','values':f'reports/{A}/actual-eight-g-third-values.json',
 'input':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/exact-two57-univariate5-input.json',
 'input_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW/report.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'CA5':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/CA-source-and-regular-quotient-second.json',
 'CB5':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/CB-source-and-regular-quotient-second.json',
 'peer5':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 'own3':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/independent-complete619-symbolic-source-postverification.json',
 'own4':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()};dump(out,'phase.json',{'phase':'provided114_order6_source_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['report','input_peer','peer5'])
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 nodes=D['source']['DAG_scalar_nodes'];yb=[list(map(Q,p)) for i,p in D['own3']['actual_YB_last2_through3']]
 assert [p[1] for p in yb]==[Q(17,12),Q(7,12)] and all(p[0]==p[2]==0 for p in yb)
 cbfirst=[list(map(Q,p)) for p in D['own4']['families']['CB']['quotient_firsts4_by2']]
 actual={};proof=[]
 for f in D['input']['families']:
  name=f['family'];P=D[name];old=D[name+'5'];check=D[name+'check'];u=f['original_unknown57_order'];E=f['original_equations57'];g=f['actual_g_original_node'];bn=f['actual_bbar4_original_nodes']
  fixed=D['FJ']['both57_source2'][name];assert fixed['unknown57']==u and fixed['equations57']==E
  J=[list(map(Q,r)) for r in fixed['fixed_J57']];Ji=[list(map(Q,r)) for r in fixed['fixed_J57_inverse']]
  assert all(sum(J[i][k]*Ji[k][j] for k in range(57))==sum(Ji[i][k]*J[k][j] for k in range(57))==Q(i==j) for i in range(57) for j in range(57))
  lower={i:tuple(map(Q,p)) for i,p in old['all57_source_jets_through5']};candidate={i:tuple(map(Q,p)) for i,p in P['all57_source_jets_through6']}
  assert list(candidate)==u and all(candidate[i][:6]==lower[i] for i in u)
  roots=E+[g]+bn;values,units=numeric_eval(nodes,roots,candidate,7,f['source_needed_original_DAG_nodes'])
  zero6={i:p+(Q(0),) for i,p in lower.items()};before,_=numeric_eval(nodes,E,zero6,7,f['source_needed_original_DAG_nodes'])
  slopes=list(map(Q,f['last2_target_slopes']));expect=[(Q(0),slopes[i-55] if i>=55 else Q(0),0,0,0,0,0) for i in range(57)]
  assert all(values[i]==p for i,p in zip(E,expect))
  rhs=[before[i][6] for i in E];new=[candidate[i][6] for i in u]
  assert list(map(str,rhs))==check['degree6_residual57_before'] and list(map(str,new))==check['new_degree6_corrections57']
  assert all(sum(J[i][j]*new[j] for j in range(57))==-rhs[i] for i in range(57))
  assert [[i,list(map(str,values[i]))] for i in E]==P['all57_source_equations_through6']
  assert [[i,list(map(str,candidate[i]))] for i in range(2,11)]==P['physical_x9_readbacks_through6']
  gp=values[g];assert gp[:3]==(0,0,0) and gp[3] and gp[4]==0 and list(map(str,gp))==P['g_through6']
  computed=[]
  for j,node in enumerate(bn):
   h1=Q(old['regular_quotient_origin0_and_first'][j][1]);h2=Q(old['new_four_regular_quotient_second_line'][j]);num=values[node]
   assert num[:4]==(0,0,0,0) and num[4]==gp[3]*h1 and num[5]==gp[3]*h2
   h3=(num[6]-gp[5]*h1)/gp[3];h=(Q(0),h1,h2,h3,0,0,0)
   assert series_mul(gp,h,7)==num
   bend=sum(cbfirst[j][k]*yb[k][3] for k in range(2)) if name=='CB' else Q(0)
   record=P['four_quotient_records'][j]
   assert record['new_own_h3']==str(h3) and record['known_common_bending']==str(bend) and record['actual_common_h3']==str(h3+bend)
   assert record['bbar_node']==node and record['bbar_through6']==list(map(str,num))
   computed.append(h3+bend)
  ids=[566400,566401,566402,566403] if name=='CA' else [330164,330289,330524,330594]
  actual.update(zip(ids,computed));proof.append({'family':name,'all57_equations_all7_coefficients_zero_or_target':True,
   'unitJ0_inverse_and_provided_correction_equation_checked':True,'physical9_readbacks_checked':True,'quotient_product_all7_coefficients_checked':True,
   'source_g_through6':list(map(str,gp)),'actual_four_common_h3':list(map(str,computed)),'unit_inverse_checks':units})
 assert [[i,str(c)] for i,c in sorted(actual.items())]==D['values']['actual_eight_common_regular_g_third_values']
 diffs=[actual[a]-actual[b] for a,b in zip([330164,330289,330524,330594],[566400,566401,566402,566403])]
 assert list(map(str,diffs))==D['values']['actual_four_oldCB_minus_newCA_third_differences']
 assert diffs==[Q(68165,3456),0,0,Q(1435,256)]
 ev={'source_manifest':inputs,'kind':'NEW higher physical order6 primary, independently verified providedcoefficients without resolving',
  'families':proof,'all114_provided_source6_coefficients_unique':'all0..5 coefficients fixed, exactsame source equations, invertibleJ0',
  'actual_eight_common_g3':[[i,str(c)] for i,c in sorted(actual.items())],
  'actual_oldCB_minus_newCA_g3':list(map(str,diffs)),
  'common_CB_bending_proof':'YB2=0 andglobalh0=0, so own-slopeh3+linear_h1 dotYB3; no quadraticH2 bending term',
  'full_ratio_times_f_source_ancestry':'literal original sparse-index sourceDAG and fixed signed regular numerator/denominator nodes retained',
  'no_new_unknown_solve_common105sixth_mixedthird_gauge3_or_HF4':True}
 ep=out/'independent-provided114-source6-and-eight-g3.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'provided114_source6_and_all57eq_through6_independently_accepted':True,
  'actual_eight_common_g3_and_bending_accepted':True,'actual_old_minus_new_g3':list(map(str,diffs)),
  'new_unknown_source_solve_or_HF4_by_reviewer':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
