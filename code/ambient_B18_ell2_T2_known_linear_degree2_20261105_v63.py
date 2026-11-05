#!/usr/bin/env python3
"""Complete degree2 homogeneous closure of actual knownL48, upper only."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
from ambient_B18_ell2_T2_first_Hilbert_source_20261105_v63 import quotient
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63'
def main(out,gate,sourcepeer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/gate;g=json.loads(gp.read_text());assert g['all_passed'] and not g['degree2_rank_started'];ip=gp.parent/'actual-complete-L48-and-degree2-rewrite-source-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text());pp=ROOT/sourcepeer;peer=json.loads(pp.read_text());assert peer['all_passed'];maxbits=0
 def guard(v):
  nonlocal maxbits
  maxbits=max([maxbits]+[max(abs(c.numerator).bit_length(),c.denominator.bit_length()) for c in v.values()]);assert maxbits<=2048;return v
 nodes=d['degree2_monomial_rewrite_nodes'];memo={}
 def nf(i):
  if i in memo:return memo[i]
  n=nodes[i]
  if n[0]=='zero_ell2':v={}
  elif n[0]=='standard':v={n[1]:F(1)}
  else:
   v={}
   for j,c in n[2]:v=guard(add(v,nf(j),F(c)))
  memo[i]=v;return v
 for roots in d['all264_NF_root_nodes']:
  for i in roots:nf(i)
 assert maxbits<=d['degree2_NF_positive_bound_bits']
 mons=[tuple(m) for m in d['degree2_monomial_order44']];mi={m:i for i,m in enumerate(mons)};rows={p:H.parse(v) for p,v in d['normalized_linear_rewrite_rows']};cols=[];labels=[]
 for p,row in rows.items():
  for var in range(9):
   q={}
   for k,c in row.items():
    target,v=divmod(k,9);m=[0]*9;m[var]+=1;m[v]+=1
    if m[8]<2:q=guard(add(q,memo[d['all264_NF_root_nodes'][target][mi[tuple(m)]]],c))
   cols.append(q);labels.append([p,var])
 assert len(cols)==432;norm={};reps={};selected=[];original=[];factors=[];receipts=[]
 for j,orig in enumerate(cols):
  v=dict(orig);a={};steps=[]
  while v and min(v) in norm:
   r=min(v);c=v[r];v=guard(add(v,norm[r],-c));a=guard(add(a,reps[r],c));steps.append([r,str(c)])
  if v:
   r=min(v);c=v[r];k=len(original);original.append(orig);selected.append(j);norm[r]=guard(scale(v,1/c));reps[r]=guard(scale(add({k:F(1)},a,-1),1/c));a={k:F(1)};receipts.append({'source_column':j,'pivot':r,'divide_by':str(c),'subtract':steps})
  assert mv(original,a)==orig;factors.append(a)
 q=len(d['degree2_linear_standard_module_monomials']);free,Pi=quotient(q,norm);assert all(not mv(Pi,v) for v in cols);h2=len(free);prefix=12+h2
 ep=out/'actual432-overlap-rank-and-known-linear-cover-degree2.json';H.dump(ep,{'all264_linear_NF_columns':[[i,srec(v)] for i,v in sorted(memo.items())],
  'all432_true_source_multiple_labels':labels,'all432_remainder_columns':[srec(v) for v in cols],
  'rank':len(norm),'fixed_original_minor_source_columns':selected,'fixed_rank_receipts':receipts,
  'normalized_degree2_relation_rows':[[i,srec(v)] for i,v in sorted(norm.items())],
  'all432_original_selected_source_factorizations':[srec(v) for v in factors],
  'known_degree2_quotient_projection_columns':[srec(v) for v in Pi],'degree2_quotient_free_standard_indices':free,
  'known_homogeneous_linear_cover_HF_through2':[6,6,h2],'known_cover_prefix_length':prefix,
  'actual_T2_HF2_not_claimed':True,'no_new_actual_quadratic_source_relations_or_coefficients':True,
  'stop_rule':'Ifprefix>20 knowncovercannotimproveactualupper20; noactualT2lowerinference. Otherwisefreshdegree3source/slotgate required.'})
 slots=sum(map(len,cols+factors+list(memo.values())+list(norm.values())+list(reps.values())+Pi))+10000;assert slots<d['whole_degree2_closure_slot_bound']<60000
 report={'id':STEM,'all_passed':True,'all432_source_multiples_and_overlaps_retained':True,'known_cover_HF_through2':[6,6,h2],
  'known_cover_prefix_length':prefix,'degree2_relation_rank':len(norm),'actual_T2_HF2_total_length_or_H0_computed':False,
  'input_sha256':{gate:H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip),sourcepeer:H.sha(pp)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},
  'owned_QQ_slot_bound':slots,'actual_max_height_bits':maxbits,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','known_cover_HF_through2','known_cover_prefix_length','degree2_relation_rank','actual_max_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
