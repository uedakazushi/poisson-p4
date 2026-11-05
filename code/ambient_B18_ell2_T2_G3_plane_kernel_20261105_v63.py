#!/usr/bin/env python3
"""Prepared bounded plane-kernel worker; run only under root GO."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63'
def main(out,gate,sourcepeer):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_distinct_G3_plane_kernel'})
 gp=ROOT/gate;gr=json.loads(gp.read_text());assert gr['all_passed'] and gr['rank_not_computed']
 ip=gp.parent/'actual378-annihilation-checks-and105by238-condition.json';assert H.sha(ip)==gr['evidence_sha256'][ip.name]
 pp=ROOT/sourcepeer;pr=json.loads(pp.read_text());assert pr['all_passed']
 A=list(map(H.parse,H.selected(str(ip.relative_to(ROOT)),'condition105_by238_columns')));assert len(A)==238
 norm={};reps={};kernels=[];selected=[];receipts=[];peak=0;maxbits=0
 cursor={};temporary_peak=0
 class GuardStop(Exception):pass
 def temporary(x):
  nonlocal temporary_peak
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());temporary_peak=max(temporary_peak,b)
  # Committed inputs are384 bits. A product has<=768 bits; adding it
  # to one committed384-bit value has<=1153 uncancelled numerator bits.
  # Division of two committed values has<=768. One offender is retained.
  assert b<=1153,('uncancelled_temporary_bound',b)
  return x
 def scalar(x):
  nonlocal maxbits
  temporary(x);b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maxbits=max(maxbits,b)
  if b>384:
   cursor['offending_uncommitted_value']=str(x);cursor['offending_height_bits']=b
   raise GuardStop('dynamic384_before_commit')
  return x
 def guard(v):
  for x in v.values():scalar(x)
  return v
 def axpy(a,b,c):
  v=dict(a)
  for i,x in b.items():
   cursor.update({'arithmetic':'accumulator_plus_source_times_scalar','term_index':i,'scalar':str(c),'uncommitted_accumulator':v})
   z=temporary(x*c);y=scalar(v.get(i,F(0))+z)
   if y:v[i]=y
   else:v.pop(i,None)
  return v
 def normalize(a,c):
  v={}
  for i,x in a.items():
   cursor.update({'arithmetic':'normalization_division','term_index':i,'scalar':str(c),'uncommitted_accumulator':v})
   v[i]=scalar(x/c)
  return v
 def image(w):
  v={}
  for j,c in w.items():v=axpy(v,A[j],c)
  return v
 try:
  for j,col in enumerate(A):
   cursor={'source_column':j,'stage':'begin'};v=dict(col);w={j:F(1)}
   while v and min(v) in norm:
    p=min(v);c=v[p];cursor.update({'stage':'reduce_target','pivot':p});v=axpy(v,norm[p],-c)
    cursor['stage']='reduce_source_representation';w=axpy(w,reps[p],-c)
   if v:
    p=min(v);c=v[p];cursor['stage']='normalize_target';v=normalize(v,c)
    cursor['stage']='normalize_source_representation';w=normalize(w,c)
    cursor['stage']='verify_new_membership';assert image(w)==v
    norm[p]=v;reps[p]=w;selected.append(j);receipts.append({'source_column':j,'pivot_row':p,'divide_by':str(c)})
   else:
    cursor['stage']='verify_new_kernel';assert w.get(j)==1 and all(k<=j for k in w) and not image(w)
    kernels.append((j,w))
   live=sum(map(len,A))+sum(map(len,norm.values()))+sum(map(len,reps.values()))+sum(len(w) for _,w in kernels)+4*(105+238)+4*238+5000
   peak=max(peak,live);assert peak<=gr['next_rank_whole_QQ_bound']<150000,('slot_guard',j,peak)
   if (j+1)%16==0:H.dump(out/'phase.json',{'phase':'active_rank_and_kernel','completed':j+1,'active_rank':len(norm),'kernel_count':len(kernels),'peak_QQ_slots':peak,'height_bits':maxbits})
 except GuardStop as e:
  # Preserve one frontier INSTEAD of final output; do not duplicate basis data.
  acc=cursor.pop('uncommitted_accumulator',{})
  with (out/'guard-frontier.json').open('w') as f:
   f.write(json.dumps({'stopped':str(e),'cursor':cursor,'work_target':srec(v),'work_source':srec(w),'partial_uncommitted_operation_accumulator':srec(acc),'source_gate_sha':H.sha(gp)},separators=(',',':'))[:-1])
   f.write(',"committed_basis_and_source_representations":[')
   for k,q in enumerate(receipts):
    if k:f.write(',')
    p=q['pivot_row'];f.write(json.dumps({'receipt':q,'basis':srec(norm[p]),'source':srec(reps[p])},separators=(',',':')))
   f.write('],"committed_kernel":[')
   for k,(free,z) in enumerate(kernels):
    if k:f.write(',')
    f.write(json.dumps([free,srec(z)],separators=(',',':')))
   f.write(']}\n')
  H.dump(out/'phase.json',{'phase':'GUARD_STOP_before_commit','cursor':cursor,'active_committed_rank':len(norm),'kernel_count':len(kernels),'temporary_peak_bits':temporary_peak,'no_retry_or_cap_raise_authorized':True})
  print(json.dumps({'all_passed':False,'guard_stop':cursor,'committed_rank':len(norm)}));return
 r=len(norm);assert r+len(kernels)==238 and r<=62
 # At each new source its coefficient on its own source index is nonzero;
 # kernel vectors have distinct free index1. Thus the combined source change
 # is triangular and invertible, certifying the entire original source space.
 assert [q['source_column'] for q in receipts]==selected
 assert set(selected).isdisjoint(j for j,_ in kernels)
 ep=out/'actual-G3-plane-condition-fixed-rank-and-kernel.json'
 fixed={'condition_shape':[105,238],'condition_rank':r,'enlarged_constant_kernel_dimension':238-r,
  'selected_original_source_columns':selected,'fixed_original_minor_rows':[q['pivot_row'] for q in receipts],
  'fixed_minor_determinant_product':[q['divide_by'] for q in receipts],
  'fixed_minor_determinant_proof':'basisrepresentations are triangular in selectedsource order withdiagonal1/pivot; normalizedleadingbasis minor determinant1, so originalminor productpivotsnonzero',
  'all_basis_membership_and_all_kernel_zero_products_verified':True,
  'source_family':'first238 activeG3sources; fullplane andprojectedlinear zero; z-linear terms retained under378NF annihilation gate',
  'actual_T2_HF2_or_total_length_computed':False,'Taylor2_coefficients_queried':False}
 # Stream certificate records, avoiding a second coefficient-array copy.
 with ep.open('w') as f:
  f.write('{')
  for i,(k,v) in enumerate(fixed.items()):
   if i:f.write(',')
   f.write(json.dumps(k)+':'+json.dumps(v,separators=(',',':')))
  f.write(',"normalized_basis_and_original_source_representations":[')
  for i,q in enumerate(receipts):
   if i:f.write(',')
   p=q['pivot_row'];f.write(json.dumps({'receipt':q,'normalized':srec(norm[p]),'original238_source_representation':srec(reps[p])},separators=(',',':')))
  f.write('],"full_constant_kernel_columns":[')
  for i,(j,w) in enumerate(kernels):
   if i:f.write(',')
   f.write(json.dumps({'free_source_index':j,'original238_source_combination':srec(w)},separators=(',',':')))
  f.write(']}\n')
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'condition_rank':r,'enlarged_constant_kernel_dimension':238-r,
  'actual_T2_HF2_or_total_length_computed':False,'new_physical_or_Taylor2_coefficients':False,
  'input_sha256':{gate:H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip),sourcepeer:H.sha(pp)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'peak_QQ_slots':peak,'max_height_bits':maxbits,'uncancelled_scalar_temporary_bound_bits':1153,'actual_temporary_peak_bits':temporary_peak,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','condition_rank','enlarged_constant_kernel_dimension','peak_QQ_slots','max_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
