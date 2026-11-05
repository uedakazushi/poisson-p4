#!/usr/bin/env python3
"""ONE 46x41 integer Bareiss certificate; known linear strand only."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};cursor={};steps=[];selectedR=[];selectedC=[];maxcommit=0;maxtemp=0
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63/';g=read(base+'report.json');assert g['all_passed'] and pins[base+'report.json']=='9f1627a8f293869e96ab9854befc9e5e9c320f1a822103416d37fd63cbda2772'
 ip=base+'pruned46-by41-integer-rank-source-and-live-bound.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and pins[base+'report.json'] in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert sha(ROOT/d['source_profile_path'])==d['source_profile_sha256'] and d['all238_A_times_G3_zero_and_G3rank62_certified'] and count([g,p,c,pins])<1000
 assert count(d)<3500 and d['whole_numeric_cap']==16784 and d['all_minors_order_le29_bits306']==306 and d['integer_product_difference_bits613']==613 and d['hard_committed384']==384 and d['hard_temporary1153']==1153
 M=[list(row) for row in d['residual_integer_matrix46_by41']];assert len(M)==46 and all(len(row)==41 for row in M);R=list(range(46));C=list(range(41));previous=1
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir())<33554432
 def stop(reason,x):
  state={'reason':reason,'cursor':cursor,'uncommitted_integer':str(x),'current_integer_matrix':M,'completed_steps':steps,'selected_rows':selectedR,'selected_cols':selectedC,'remaining_rows':R,'remaining_cols':C,'source_input_sha256':pins[ip],'no_retry_or_cap_raise':True};assert 2*count(d)+count([state,M,steps,g,p,c,pins])+1000<16784;dump('guard-frontier.json',state);raise RuntimeError(reason)
 def commit(x):
  nonlocal maxcommit
  b=abs(x).bit_length();maxcommit=max(maxcommit,b)
  if b>306:stop('source306committedbound',x)
  return x
 def temp(x):
  nonlocal maxtemp
  b=abs(x).bit_length();maxtemp=max(maxtemp,b)
  if b>613:stop('source613temporarybound',x)
  return x
 dump('phase.json',{'phase':'authorized_fixed46x41_known_integer_test_before_rank','source_gate_sha256':pins[base+'report.json'],'hard_guards':[384,1153],'derived_guards':[306,613],'no_physical_or_jet_query':True})
 while len(steps)<28:
  found=next(((i,j) for i in R for j in C if M[i][j]),None)
  if found is None:break
  r,s=found;pv=commit(M[r][s]);newR=[i for i in R if i!=r];newC=[j for j in C if j!=s]
  for i in newR:
   for j in newC:
    cursor={'stage':'fixed_integer_Schur','step':len(steps),'pivot_row46':r,'pivot_col41':s,'target_row46':i,'target_col41':j};num=temp(temp(pv*M[i][j])-temp(M[i][s]*M[r][j]));q,rem=divmod(num,previous)
    if rem:stop('source_exact_integer_division_failed',num)
    M[i][j]=commit(q)
  steps.append({'pivot_row46':r,'pivot_col41':s,'pivot_integer':str(pv),'exact_previous_divisor':str(previous)});selectedR.append(r);selectedC.append(s)
  for i in newR:M[i][s]=0
  for j in newC:M[r][j]=0
  R=newR;C=newC;previous=pv
  dump('phase.json',{'phase':'fixed_integer_Schur_prefix','completed_pivot_steps':steps,'remaining_shape':[len(R),len(C)],'maximum_committed_bits':maxcommit,'maximum_temporary_bits':maxtemp})
 cursor={'stage':'all_remaining_Schur_entries_zero','completed_steps':len(steps)}
 for i in R:
  for j in C:
   if M[i][j]:stop('known_rank_upper_or_frame_inconsistent',M[i][j])
 rank=len(steps);assert rank<=28 and (rank==0 or previous)
 full=162+rank;nullity=252-full;minor=previous if rank else 1
 kernel={'residual_free_columns41':C,'residual_generator_recipe':'forfreej: e_j-E_selected*adj(U)*originalR[selected_rows,j]/det(U), U=originalR[selected_rows,selected_cols]','selected_rows46':selectedR,'selected_cols41':selectedC,'determinant_literal':str(minor),'additional_zero_columns49_original252':d['additional_zero_source_columns49'],'reverse162_unit_source_lift_recipe':'reverseprofilefixedpivots:row-singletonsetsremovedsourcecoordinate0;column-singletonsetsitminus(sumremainingrowcoeff*currentcoordinates)/literalpivot. Coefficientsaretheunchangedliteralprofilematrix;noimplicitfunctionquery.','complete_active_kernel_dimension':nullity,'extra_dimension_over_known_G3_62':nullity-62}
 certificate={'residual_rank':rank,'active_G2_linear_strand_rank':full,'active_kernel_dimension':nullity,'fixed_integer_pivot_steps':steps,'selected_residual_rows46':selectedR,'selected_residual_columns41':selectedC,'selected_original340_rows':[d['original340_row_indices46'][i] for i in selectedR],'selected_original252_columns':[d['original252_column_indices41'][j] for j in selectedC],'fixed_residual_integer_minor_nonzero':True,'fixed_residual_integer_minor':str(minor),'all_remaining_Schur_entries_exact_zero':True,'all41_column_span_recipe':'selectedoriginalresidualcolumns timesadj(U)*originalR[selected_rows,j]/det(U);zeroSchurverifieseveryrow','full_nonzero_minor_factored_recipe':'product162literalprofilepivots *residualinteger_minor/productselectedrowclearingD;orderedrow/columnLaplacecertificate,noevaluatedlargeQQdeterminant','complete_kernel_recipes':kernel,'active_coefficient_degree1_kernel_equals_G3_image62':rank==28,'full918_strand_kernel_equals_G3linear_image242_if_tensor_source_premise':rank==28,'not_general_associated_graded_exactness':True,'b1_functional_annihilation_still_separate':True,'no_source_physical_gauge_or_higher_jet_query':True,'maximum_committed_integer_bits':maxcommit,'maximum_temporary_integer_bits':maxtemp}
 owned=2*count(d)+count([M,steps,selectedR,selectedC,R,C,kernel,certificate,g,p,c,pins])+2000;assert owned<16784
 ep='fixed46-by41-integer-Schur-rank-and-full-strand-kernel-recipes.json';dump(ep,certificate)
 report={'id':STEM,'all_passed':True,'known_fixed_integer_matrix_only':True,'residual_rank':rank,'active_linear_strand_rank':full,'active_kernel_dimension':nullity,'active_kernel_equals_known_G3':rank==28,'conditional_full918_exactness':rank==28,'actual_T2_higher_relation_HF_or_length_computed':False,'committed_integer_bits':maxcommit,'temporary_integer_bits':maxtemp,'owned_numeric_bound':owned,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({'all_passed':True,'residual_rank':rank,'active_rank':full,'kernel_dim':nullity,'sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
