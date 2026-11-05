#!/usr/bin/env python3
"""Known A*G3 zero/source/minor bound gate; never query residual rank."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);hs={};maximum=0
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 pb='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_SOURCE_PROFILE_20261105_v63/';p=read(pb+'report.json');assert p['all_passed'] and hs[pb+'report.json']=='effebaaa9b673f6a19bdde737ea9c06819b96db755f9cbfac136787a22ec9c6a'
 d=read(pb+'literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json');assert len(d['fixed_singleton_unit_pivots'])==162 and d['remaining_rank_upper']==28
 A=[H.parse(v) for v in d['literal_matrix_columns']];gp='reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_INPUT_GATE_20261105_v63/actual-full-G3-linear-zero3-split-and-active-kernel-input.json';G=[H.parse(v) for v in field(gp,'actual_active_G3linear252_by238_columns')];assert len(A)==252 and len(G)==238
 zero=(1,0)
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());assert maximum<=384;return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else zero
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else zero
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 for col in G:
  q={}
  for j,c in col.items():
   for i,x in A[j].items():q[i]=plus(q.get(i,zero),mul(ex(c),ex(x)))
 H.dump(out/'phase.json',{'phase':'positive_bound_before_all238_known_A_G3_zero_checks','positive_QQ_prior':maximum,'no_residual_rank':True})
 actual=0
 def ck(x):
  nonlocal actual
  actual=max(actual,abs(x.numerator).bit_length(),x.denominator.bit_length());assert actual<=maximum;return x
 for col in G:
  q={}
  for j,c in col.items():
   for i,x in A[j].items():
    v=ck(q.get(i,F(0))+ck(c*x))
    if v:q[i]=v
    else:q.pop(i,None)
  assert not q
 rankp='reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63/actual-full-G3linear-active-fixed-rank-and-kernel.json';assert field(rankp,'active_rank')==62
 assert read('reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63_REVIEW/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json')['all_passed']
 R=d['remaining_nonzero_row_indices'];C=d['remaining_nonzero_column_indices'];assert [len(R),len(C)]==[46,41];ri={r:i for i,r in enumerate(R)};ci={c:i for i,c in enumerate(C)};clear={r:(int(D),int(U),bits) for r,D,U,bits in d['remaining_row_integer_D_L1_bounds']};M=[[0]*41 for _ in range(46)]
 for j,col in d['remaining_literal_matrix_columns']:
  for i,s in col:
   x=F(s);D,U,bits=clear[i];assert D%x.denominator==0;n=x.numerator*(D//x.denominator);assert abs(n)<=U and abs(n).bit_length()<=17;M[ri[i]][ci[j]]=n
 assert sum(bool(x) for row in M for x in row)==102
 for i,r in enumerate(R):assert sum(map(abs,M[i]))==clear[r][1]
 top=sorted((bits for D,U,bits in clear.values()),reverse=True);assert sum(top[:28])==304 and sum(top[:29])==306
 zeroC=sorted(set(range(252))-{v['column252'] for v in d['fixed_singleton_unit_pivots']}-set(C));assert len(zeroC)==49
 payload={'residual_integer_matrix46_by41':M,'original340_row_indices46':R,'original252_column_indices41':C,'row_integer_scaling_D46':[str(clear[r][0]) for r in R],'fixed_singleton_unit_pivots162':d['fixed_singleton_unit_pivots'],'additional_zero_source_columns49':zeroC,'known_active_G3_image_rank62':62,'known_full_G3_image_rank242':242,'known_total_active_rank_upper190':190,'residual_rank_upper28':28,'integer_entry_bits17':17,'all_minors_order_le29_bits306':306,'integer_product_difference_bits613':613,'hard_committed384':384,'hard_temporary1153':1153,'time_seconds':60,'memory_bytes':2000000000,'whole_numeric_cap':16784,'whole_RAW_cap':33554432,
 'all238_A_times_G3_zero_and_G3rank62_certified':True,'source_profile_path':pb+'literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json','source_profile_sha256':hs[pb+'literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json'],
 'tensor_source_linkage':'peeredexactactiveunitresolution tensoredK(zA2,zA3,zA6);linearpointunitsretainthissourceframe.ConditionalactiveH2internal3=0 =>fullH2internal3=0,notallordergradedexactness.',
 'failure_or_success_live_bound':{'payload_double_including_meta_upper':7000,'working_integer_matrix':1886,'pivot_vectors_cursor_and_steps_upper':500,'output_or_STOP_serialization_upper':4500,'gate_peer_custody_metadata_upper':1000,'sum':14886},'rank_test_NOT_run':True}
 assert count(payload)<3500 and sum(payload['failure_or_success_live_bound'][k] for k in payload['failure_or_success_live_bound'] if k!='sum')==14886<16784
 ep=out/'pruned46-by41-integer-rank-source-and-live-bound.json';H.dump(ep,payload)
 addendum={'id':'AMBIENT_B18_ELL2_T2_G2_STRAND_SOURCE_RESOURCE_ADDENDUM_20261105_v63','all_passed':True,'old_profile_sha256':hs[pb+'report.json'],'old_source_count':16784,'omitted_old_G3_certificate':1664,'graph_sets_conservative_bound':2008,'corrected_source_phase_upper':20456,'hard_source_cap_unchanged':150000,'rank_worker_cap_unchanged':16784,'old_bytes_unchanged':True};H.dump(out/'source-profile-count-addendum.json',addendum)
 owned=count([p,d,A,G,R,C,clear,M,payload,addendum])+18000;assert owned<150000
 report={'id':STEM,'all_passed':True,'known_source_product_checks_only_no_residual_rank':True,'all238_A_G3_zero':True,'positive_product_prior_bits':maximum,'actual_known_product_bits':actual,'source_owned_bound':owned,'rank_worker_numeric_cap':16784,'rank_worker_live_upper':14886,'derived_integer_committed_bound':306,'derived_integer_temporary_bound':613,'hard_integer_guards':[384,1153],'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep),'source-profile-count-addendum.json':H.sha(out/'source-profile-count-addendum.json')},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_product_prior_bits','actual_known_product_bits','source_owned_bound','rank_worker_live_upper']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
