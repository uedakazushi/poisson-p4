#!/usr/bin/env python3
"""STOP-only four-row source/code/whole custody peer; no signed row arithmetic."""
from pathlib import Path
from fractions import Fraction as Q
import ast,hashlib,json,resource,signal,subprocess,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_INPUT_GATE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_INPUT_GATE_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='6354180e36b18d7d2f14ee45fbfacd10edd48a3efbb258eafeddc7db7068ae06'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'pruned-four-point-dualrow-STOP-only-input.json')
 old='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_INPUT_GATE_20261105_v63/report.json';oldr=rd(old);assert oldr['all_passed'] and oldr['raw_row_construction_prior_bits']==424 and not oldr['row_stage_ready_with_unchanged384']
 kp='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63/actual396-Kfull-P5-and-native23-existing-cache-values.json';K=fld(kp,['K99_two_channels_rr_rell'])['K99_two_channels_rr_rell'];assert d['K99_two_channels_rr_rell']==K
 kr='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63_REVIEW/report.json';assert rd(kr)['all_passed'] and pins[kr]=='251c318156c9f6868f3c5ae68c24bcd40b7c9261cf65757e3d1f88b7dde33386'
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';z=fld(mu,['normalized_basis_rows','normalized_basis_original_selected_source_representations','actual_selected_Y0_columns102'])
 assert d['normalized99_basis']==z['normalized_basis_rows'] and d['normalized99_selected_representations']==z['normalized_basis_original_selected_source_representations'] and d['selected99_Y0_columns102']==z['actual_selected_Y0_columns102']
 mr='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json';assert rd(mr)['all_passed'] and pins[mr]=='6428e2178c41cfd7cfcf271402c4d13898879414a0b5129fd409e3c97300d460'
 norm={i:sp(v) for i,v in d['normalized99_basis']};rep={i:sp(v) for i,v in d['normalized99_selected_representations']};Y=list(map(sp,d['selected99_Y0_columns102']));KF=[[[Q(x) for x in pair] for pair in c] for c in K]
 assert len(norm)==len(rep)==len(Y)==len(KF)==99 and set(norm)==set(rep) and set(range(102))-set(norm)==set(d['free_Y3_point_axes'])=={1,36,38}
 assert all(c.get(i)==1 and all(j>=i for j in c) for i,c in norm.items()) and all(0<=j<99 for c in rep.values() for j in c)
 initial=[x for c in KF for pair in c for x in pair]+[x for A in [norm,rep] for c in A.values() for x in c.values()]+[x for c in Y for x in c.values()];height=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in initial);assert height==60==r['initial_cached_operand_maximum_bits'];del initial
 assert d['four_row_order']==[['channel0','rr'],['channel0','rell'],['channel1','rr'],['channel1','rell']]
 assert d['STOP_only_no_completion_bound'] and d['prior_coarse424_is_not_execution_growth_promise'] and d['hard_committed_bits']==384 and d['hard_temporary_bits']==1153
 assert d['derived_binary_mul_or_div_bits']==768 and d['derived_binary_add_or_sub_bits']==769 and d['single_uncommitted_offender_maxbits']==769
 assert count(d)==3176 and count([KF,norm,rep,Y])==3164 and 2*3176+3164+6000==15516==r['explicit_future_live_bound']<r['future_numeric_cap']==20406
 cb='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_PREPARED_CUSTODY_V2_20261105_v63/';cu=rd(cb+'report.json');assert cu['all_passed'] and pins[cb+'report.json']=='d62bf895364ac55280f70e9c5d6f9116c389833190822d56adc12b879862cc80'
 for p,h in cu['literal_input_sha256'].items():assert sha(ROOT/p)==h
 closure=rd(cb+'guard-state-and-all-private-live-closure.json');ledger=rd(cb+'STOP-only-preoperation-and-ownership-ledger.json')
 assert pins[cb+'guard-state-and-all-private-live-closure.json']==cu['guard_private_live_closure_sha256'] and pins[cb+'STOP-only-preoperation-and-ownership-ledger.json']==cu['ownership_ledger_sha256']
 assert closure['single_full_state_numeric_upper']==603+198+198+99+10+4+1+1==1114
 assert closure['three_state_copies_for_live_state_string_conversion_and_RAW']==3*1114==3342
 assert closure['guard_path_total_reserve']==3342+1000+198+20+200==4760<6000
 assert closure['guard_whole_upper']==9516+4760==14276<15516<20406 and cu['guard_path_numeric_upper']==14276
 worker=ROOT/cu['worker'];assert sha(worker)==cu['worker_sha256']=='5956d12678224828b59d3f37ba5a0c944f447c6e3932692f70acb9337369a96f'
 code=worker.read_text();tree=ast.parse(code);func={n.name:n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef)};binary=ast.get_source_segment(code,func['binary']);stop=ast.get_source_segment(code,func['stop'])
 assert "if bits>384:stop('384bit_STOP_before_commit',q,op,a,b)" in binary and binary.index('if bits>384')<binary.index('return q')
 assert all(s in binary for s in ['q=a*b','q=a+b','q=a-b','<=384']) and 'else:raise AssertionError(op)' in binary
 assert all(s in stop for s in ["'operation_cursor':dict(cursor)","'operand_a':str(a)","'operand_b':str(b)","'uncommitted_fraction':str(result)","'current_accumulator_before_operation':str(acc)","'completed_point_rows':completed","'all_current_rhs_pivots':rec(rhs)","'current_backsubstitution_row':rec(current)","'fixed_pivot_order':sorted(norm)","'input_sha256':pins","'no_retry_or_cap_raise':True"])
 # Source-AST control proof, not execution of any response recurrence.
 assert "term=binary(K[j][t][pair],x,'multiply');cursor['pending']='add_rhs';nxt=binary(acc,term,'add');acc=nxt" in code
 assert "term=binary(x,current.get(j,F(0)),'multiply');cursor['pending']='subtract';nxt=binary(acc,term,'subtract');acc=nxt" in code
 assert "term=binary(current.get(j,F(0)),x,'multiply');cursor['pending']='add';nxt=binary(acc,term,'add');acc=nxt" in code
 assert "for pivot in reversed(sorted(norm))" in code and "if j==pivot:continue" in code and "current[pivot]=acc" in code
 assert "expected=rhs[col] if kind=='normalized_basis' else K[col][t][pair]" in code
 helpresult=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helpresult.returncode==0 and hashlib.sha256(helpresult.stdout).hexdigest()==cu['safe_help_stdout_sha256']
 assert cu['helpers']=='NONE:standardlibraryonly' and cu['hard_guards']==[384,1153] and cu['binary_temporary_bound']==769 and cu['future_numeric_cap']==20406
 rawset=set(cu['literal_input_sha256'])|{cu['worker'],cb+'report.json',cb+'STOP-only-preoperation-and-ownership-ledger.json',cb+'guard-state-and-all-private-live-closure.json'};rawbytes=sum((ROOT/p).stat().st_size for p in rawset);rawbound=rawbytes+1000000+2000000+100000;assert rawbound<33554432 and rawbytes<150804
 proof={'all_K396_literal_values_and_point99_basis_representations_bound_to_originals':True,'source_F2_pointorder102_and_selected99_frame_unchanged':True,'free_axes1_36_38_fixed_not_unknown_values':True,'normalized_pivot_diagonal1_and_upper_support':True,'fourrow_equation_proof':'Eachnormcolumn=Yselected*rep, so rhs_p=K*rep_p. Descendingnormalizedpivotsolvesrow*norm=rhs withfree1,36,38=0. All99sourceand99normmembershipsperrowarestillrequiredonsuccess.','initial_operand_height60':True,'STOP_before_use_proof':'Everymultiplyreturnsfrombinaryonlyafter<=384commitcheck;onlythenmayitbeadded/subtracted. Everynewaccumulatorisassignedonlyafteritsbinarysumcheck. Binaryinputsnumerators/denominators<=384 implyproductintegers<=768,sumcrossproducts<=769;normalizedoffendercannotgrow. Hard1153unchanged. No424completionclaim.','exact_frontier_scope':'AtanyheightSTOP max3completedrows;all99rhs,current99,pivotorder99,cursor10,operands/result/acc4 plusflags. State1114;3copies3342+wrappers/meta1418=4760<6000. Completedlastrowmembershipusesnewcurrentand3completed,not4.','guard_path_whole14276_and_prearith15516':True,'success_path_fourrows_only':True,'future_numeric_cap20406':True,'literal_RAW_bytes':rawbytes,'whole_RAW_upper_with_metadata_output':rawbound,'custody_whole_RAW3150804':cu['whole_RAW_upper'],'immutable_worker_and_safehelp_literalPASS':True,'old99f3_coarse424_and_oldcustody60c95_preserved':True,'no_signed_pointrow_alpha_G3_RHS_b1_upper_or_higher_values_evaluated':True}
 ep='independent-STOP-fourrow-source-code-private-state-and-whole-bound.json';put(ep,proof)
 owned=count([r,d,oldr,K,z,norm,rep,Y,KF,cu,closure,ledger,proof,pins])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'SOURCE_code_custody_STOP_only_no_signedrows':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(out/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);assert (out/'report.json').stat().st_size<1000000;print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned}))
if __name__=='__main__':main()
