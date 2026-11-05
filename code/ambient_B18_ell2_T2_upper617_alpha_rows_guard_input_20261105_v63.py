#!/usr/bin/env python3
"""Prune a STOP-only four-point-dual-row job; no signed row arithmetic."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 old='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_INPUT_GATE_20261105_v63/';r=read(old+'report.json');assert r['all_passed'] and r['raw_row_construction_prior_bits']==424 and not r['row_stage_ready_with_unchanged384']
 ip=old+'Kfull396-to-full918-alpha-factor-input-and-positive-stage-bounds.json'
 raw={k:field(ip,k) for k in ['K99_two_channels_rr_rell','selected99_Y0_columns102','normalized99_basis','normalized99_selected_representations','free_Y3_point_axes']}
 for p in ['reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json']:assert read(p)['all_passed']
 K=[[[F(x) for x in q] for q in ch] for ch in raw['K99_two_channels_rr_rell']];norm={i:H.parse(v) for i,v in raw['normalized99_basis']};rep={i:H.parse(v) for i,v in raw['normalized99_selected_representations']};Y=list(map(H.parse,raw['selected99_Y0_columns102']))
 assert len(K)==len(norm)==len(rep)==len(Y)==99 and set(norm)==set(rep) and set(range(102))-set(norm)==set(raw['free_Y3_point_axes'])=={1,36,38}
 for i,v in norm.items():assert v.get(i)==1 and all(j>=i for j in v)
 operands=[x for ch in K for q in ch for x in q]+[x for m in [norm,rep] for v in m.values() for x in v.values()]+[x for v in Y for x in v.values()]
 initial=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in operands);assert initial<=384
 payload={**raw,'four_row_order':[['channel0','rr'],['channel0','rell'],['channel1','rr'],['channel1','rell']],'equations':'s_t,pair*Ystar0=K[:,t,pair];s_free1,36,38=0;rhs_p=sum_j K[j,t,pair]*repr[p,j];s_p=rhs_p-sum_(j>p)norm[p,j]*s_j.','STOP_only_no_completion_bound':True,'prior_coarse424_is_not_execution_growth_promise':True,'hard_committed_bits':384,'hard_temporary_bits':1153,'derived_binary_mul_or_div_bits':768,'derived_binary_add_or_sub_bits':769,'each_product_commit_guard_before_subsequent_add':True,'single_uncommitted_offender_maxbits':769,'no_alpha918_clearing_G3dots_RHSZ2_or_response_values':True}
 n=count(payload);parsed=count([K,norm,rep,Y]);whole=2*n+parsed+6000;assert whole<20406
 ep=out/'pruned-four-point-dualrow-STOP-only-input.json';H.dump(ep,payload)
 report={'id':STEM,'all_passed':True,'source_only_no_signed_dualrow_operation':True,'STOP_only_execution_method_not_promised_completion':True,'coarse424_blocked_history_preserved':True,'initial_cached_operand_maximum_bits':initial,'hard_committed_bits':384,'hard_temporary_bits':1153,'binary_temporary_bound':769,'uncommitted_offender_bits_bound':769,'pruned_numeric_records':n,'parsed_source_numeric_records':parsed,'explicit_future_live_bound':whole,'future_numeric_cap':20406,'future_RAW_output_cap':2000000,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','initial_cached_operand_maximum_bits','pruned_numeric_records','explicit_future_live_bound','future_numeric_cap','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
