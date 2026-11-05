#!/usr/bin/env python3
"""Prepared tiny actual G3-plane quadratic values/rank; root GO required."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_PLANE_QUADRATIC_20261105_v63'
def main(out,gate,sourcepeer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_signed_known_resolution_quad'})
 gp=ROOT/gate;gr=json.loads(gp.read_text());assert gr['all_passed'] and not gr['signed_actual_quadratic_values_evaluated']
 ip=gp.parent/'actual-enlarged-G3quadratic-factored-source-and-known-NF-input.json';assert H.sha(ip)==gr['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 pp=ROOT/sourcepeer;assert json.loads(pp.read_text())['all_passed']
 A=list(map(H.parse,d['actual_enlarged_selected126_correction_linear_columns']));P=list(map(H.parse,d['actual_selected126_rawD31_Pi42_linear_columns']));NF=list(map(H.parse,d['PiY_L48_degree2_NF42_times20_columns']));pairs={tuple(v):i for i,v in enumerate(d['degree2_active_pair_order20'])}
 # Positive input envelopes only; no signed values enter this height proof.
 Dall=math.lcm(*(int(q[0]) for q in d['positive_bound_per_enlarged_common_D_U']));U=max(int(q[1])*(Dall//int(q[0])) for q in d['positive_bound_per_enlarged_common_D_U'])
 minorbits=6*U.bit_length()+10;committed_bound=minorbits+Dall.bit_length();temporary_bound=3*committed_bound+1
 assert committed_bound<=535 and temporary_bound<2048
 cols=[];norm={};reps={};original=[];selected=[];factors=[];receipts=[];maxbits=0;temporary_peak=0;cursor={}
 class Stop(Exception):pass
 def guard(x,commit=True):
  nonlocal maxbits,temporary_peak
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());temporary_peak=max(temporary_peak,b)
  if b>2048:
   cursor['offending_uncommitted_value']=str(x);cursor['height_bits']=b;raise Stop('2048_beforecommit')
  if commit:maxbits=max(maxbits,b)
  return x
 def axpy(a,b,c):
  v=dict(a)
  for i,x in b.items():
   cursor['term_index']=i;cursor['uncommitted_accumulator']=v
   y=guard(v.get(i,F(0))+guard(x*c,False))
   if y:v[i]=y
   else:v.pop(i,None)
  return v
 def image(basis,w):
  v={}
  for j,c in w.items():v=axpy(v,basis[j],c)
  return v
 try:
  for j,a in enumerate(A):
   cursor={'stage':'signed_quad','kernel_source_index':j};q={}
   for k,x in a.items():
    source,u=divmod(k,6)
    for i,y in P[source].items():
     row,v=divmod(i,6);pair=tuple(sorted((u,v)))
     if pair==(5,5):continue
     q=axpy(q,NF[row*20+pairs[pair]],guard(-x*y,False))
   D,Uj=map(int,d['positive_bound_per_enlarged_common_D_U'][j])
   assert all(D%c.denominator==0 and abs(c.numerator)*(D//c.denominator)<=Uj for c in q.values())
   cols.append(q)
  for j,col in enumerate(cols):
   cursor={'stage':'rank','relation_column':j};v=dict(col);a={};steps=[]
   while v and min(v) in norm:
    p=min(v);c=v[p];v=axpy(v,norm[p],-c);a=axpy(a,reps[p],c);steps.append([p,str(c)])
   if v:
    p=min(v);c=v[p];idx=len(original)
    n={i:guard(x/c) for i,x in v.items()};w=axpy({idx:F(1)},a,-1);w={i:guard(x/c) for i,x in w.items()}
    assert image(original+[col],w)==n
    norm[p]=n;reps[p]=w;original.append(col);selected.append(j);a={idx:F(1)}
    receipts.append({'source_kernel_index':j,'pivot':p,'divide_by':str(c),'subtract':steps})
   assert image(original,a)==col;factors.append(a)
 except Stop:
  acc=cursor.pop('uncommitted_accumulator',{})
  H.dump(out/'guard-frontier.json',{'cursor':cursor,'partial_accumulator':srec(acc),'completed_relation_columns':[srec(v) for v in cols],
   'committed_basis':[[p,srec(v),srec(reps[p])] for p,v in norm.items()],'selected':selected,'factors':[srec(v) for v in factors],
   'no_retry_or_cap_raise':True})
  print(json.dumps({'all_passed':False,'guard_stop':cursor}));return
 rank=len(norm);assert rank<=6 and maxbits<=committed_bound
 free=[j for j in range(len(cols)) if j not in selected]
 for j in free:
  w={j:F(1)}
  for k,c in factors[j].items():w[selected[k]]=w.get(selected[k],F(0))-c
  assert not image(cols,w)
 slots=sum(map(len,A+P+NF+cols+factors+original))+sum(map(len,norm.values()))+sum(map(len,reps.values()))+18000
 assert slots<=gr['next_whole_slot_bound']<60000
 ep=out/'actual218-G3-quadratic-initial-classes-and-fixed-rank.json';H.dump(ep,{
  'actual218_projected_quadratic_initial_relation_columns':[srec(v) for v in cols],
  'rank':rank,'selected_original_kernel_relation_indices':selected,'fixed_original_minor_rows':[q['pivot'] for q in receipts],
  'fixed_minor_determinant_product':[q['divide_by'] for q in receipts],'fixed_rank_receipts':receipts,
  'normalized_basis':[[p,srec(v)] for p,v in sorted(norm.items())],
  'basis_source_representations':[[p,srec(reps[p])] for p in sorted(norm)],
  'all218_factorizations_in_selected_original_relation_basis':[srec(v) for v in factors],
  'kernel_free_source_indices':free,'kernel_factored_recipe':'for eachfreej, e_j-sum_k factor[j][k]*e_selected[k]; allzeroimages checked andfree-index blockI',
  'all_products_and_positive_bounds_checked':True,'known_quadratic_cover_standard_dimension':6-rank,
  'actual_HF2_or_total_length_not_computed':True,'new_physical_or_implicit_coefficients':False,
  'source_scope':'true quadraticinitial relations from enlarged218 activeG3constant family only; all378 unknownunitderivative products vanishNF, actualsourcefirstfullztermsretained',
  'global_positive_output_D_U':[str(Dall),str(U)],'integer_minor_height_bound':minorbits,'source_normalized_height_bound':committed_bound})
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'known_quadratic_initial_relation_rank':rank,'known_quadratic_cover_standard_dimension':6-rank,
  'actual_HF2_total_length_or_H0_not_computed':True,'new_physical_or_implicit_coefficients':False,
  'input_sha256':{gate:H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip),sourcepeer:H.sha(pp)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'owned_QQ_slots':slots,'max_committed_height_bits':maxbits,
  'source_derived_committed_bound_bits':committed_bound,'source_derived_uncommitted_bound_bits':temporary_bound,
  'actual_temporary_peak_bits':temporary_peak,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','known_quadratic_initial_relation_rank','known_quadratic_cover_standard_dimension','owned_QQ_slots','max_committed_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
