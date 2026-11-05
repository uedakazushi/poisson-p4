#!/usr/bin/env python3
"""Prepared first-cache lower4-by9 values; requires parent GO and source peer."""
import argparse,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_LOWER_VALUES_20261105_v63'
def main(out,gate,peer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/gate;g=json.loads(gp.read_text());pp=ROOT/peer;p=json.loads(pp.read_text());assert g['all_passed'] and p['all_passed'];assert H.sha(gp) in json.dumps(p),'source peer must literally pin this gate'
 assert not g['new_rank_or_actual_lower_value_job']
 ip=gp.parent/'exact171-decomposition-162-zero-9-lower-first-source-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 assert d['all_eight_dual_Q_annihilations'] and not d['source_actual_lower_values_evaluated']
 hs={gate:H.sha(gp),peer:H.sha(pp),str(ip.relative_to(ROOT)):H.sha(ip)}
 for path,sha in g['input_sha256'].items():assert H.sha(ROOT/path)==sha
 rank_bound=d['integer4by9_rank_bound_bits'];Dall=int(d['common_integer_lower_D']);Uall=int(d['common_integer_lower_abs_U']);assert Dall.bit_length()<=rank_bound<=384
 # Input is already source-bound; all dot-term denominators divide Dall and
 # their absolute cleared sum is bounded by the saved positive envelope.
 committed=max(rank_bound,g['cached_first_prior_bits']);temporary=3*committed+1;assert committed<=384 and temporary<=1153
 # Source's10000 margin also covers its three appended bound metadata fields.
 # Keep9000 for this precheck; the explicit live-output check below uses5000.
 assert count(d)+9000<=g['future_lower_owned_bound']<60000
 H.dump(out/'phase.json',{'phase':'authorized_cache_lower_only_before_signed_dots','source_gate_sha256':H.sha(gp),'source_peer_sha256':H.sha(pp),'committed_prior_bits':committed,'temporary_prior_bits':temporary,'whole_numeric':g['future_lower_owned_bound']})
 Hr=list(map(H.parse,d['Hstar1_r_selected99_columns']));He=list(map(H.parse,d['Hstar1_ell_selected99_columns']));cols=[];integers=[];maxbits=0
 def guard(x):
  nonlocal maxbits
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maxbits=max(maxbits,b);assert b<=committed,('source_height_guard',b,committed);return x
 for j,r in enumerate(d['remaining9_sources']):
  br=H.parse(r['b1_r_selected99']);be=H.parse(r['b1_ell_selected99']);col={}
  for a,ell in enumerate([[0,1,4],[11,0,270]]):
   for mixed in [False,True]:
    s=F(0)
    for i in sorted(set(br)|set(be)):
     for k,l in enumerate(ell):
      if not l:continue
      s=guard(s-guard(F(l)*Hr[i].get(k,F(0))*(be if mixed else br).get(i,F(0))))
      if mixed:s=guard(s-guard(F(l)*He[i].get(k,F(0))*br.get(i,F(0))))
    row=2*a+int(mixed);Dj,Uj=map(int,d['all36_lower_positive_D_U'][4*j+row]);assert Dj%s.denominator==0 and abs(s.numerator)*(Dj//s.denominator)<=Uj
    if s:col[row]=s
  z={i:guard(x*Dall) for i,x in col.items()};assert all(x.denominator==1 and abs(x)<=Uall for x in z.values());cols.append(col);integers.append(z)
 # Rank over QQ of the integer-cleared original columns. The usual normalized
 # elimination entries/representations are ratios of minors of order<=4.
 norm={};reps={};original=[];selected=[];factors=[];receipts=[]
 def axpy(v,w,c):
  q=dict(v)
  for i,x in w.items():
   y=guard(q.get(i,F(0))+guard(x*c))
   if y:q[i]=y
   else:q.pop(i,None)
  return q
 def image(b,w):
  q={}
  for j,c in w.items():q=axpy(q,b[j],c)
  return q
 for j,v0 in enumerate(integers):
  v=dict(v0);a={};steps=[]
  while v and min(v) in norm:
   k=min(v);c=v[k];v=axpy(v,norm[k],-c);a=axpy(a,reps[k],c);steps.append([k,str(c)])
  if v:
   k=min(v);c=v[k];idx=len(original);original.append(v0);selected.append(j);norm[k]={i:guard(x/c) for i,x in v.items()};reps[k]={i:guard(x/c) for i,x in axpy({idx:F(1)},a,-1).items()};a={idx:F(1)}
   assert image(original,reps[k])==norm[k];receipts.append({'source9_index':j,'original418_source':d['remaining9_sources'][j]['source_index'],'pivot':k,'divide_by':str(c),'subtract':steps})
  assert image(original,a)==v0;factors.append(a)
 rank=len(original);rows=[r['pivot'] for r in receipts];det=F(1) if not rank else F(0)
 if rank:
  for perm in itertools.permutations(range(rank)):
   sign=(-1)**sum(perm[i]>perm[j] for i in range(rank) for j in range(i+1,rank));term=F(sign)
   for i in range(rank):term=guard(term*integers[selected[perm[i]]].get(rows[i],F(0)))
   det=guard(det+term)
  assert det
 kernels=[]
 for j in range(9):
  if j in selected:continue
  w={j:F(1)}
  for k,x in factors[j].items():w[selected[k]]=w.get(selected[k],F(0))-x
  assert not image(integers,w);kernels.append(srec(w))
 assert maxbits<=committed and count([d,cols,integers,norm,reps,original,factors,receipts,kernels])+5000<g['future_lower_owned_bound']
 ep=out/'actual-four-lower-by-nine-source-values-and-fixed-rank.json';H.dump(ep,{'original418_sources':[r['source_index'] for r in d['remaining9_sources']],
  'actual_lower4_by9_columns':[srec(v) for v in cols],'common_integer_denominator':str(Dall),'integer_columns':[srec(v) for v in integers],
  'known_lower_quadratic_relation_rank':rank,'selected_original9_columns':selected,'fixed_original_minor_rows':rows,'fixed_integer_minor_determinant':str(det),
  'fixed_rank_receipts':receipts,'all9_original_basis_factorizations':[srec(v) for v in factors],'normalized_basis':[[i,srec(v),srec(reps[i])] for i,v in norm.items()],
  'kernel_source9_combinations':kernels,'all_factor_and_kernel_identities_verified':True,
  'scope':'Only four lower-dual images. Rank is a lower bound on actual quadratic relation rank, giving actual HF2<=6-rank. Upper components remain unknown; no actual HF2 equality/total/T2 torsion inference.',
  'upper_components_not_zero_filled':True,'new_source_physical_gauge_or_implicit_coefficients':False})
 assert sum(p.stat().st_size for p in out.iterdir())<g['future_lower_RAW_bound']<32*1024*1024
 report={'id':STEM,'all_passed':True,'known_lower_quadratic_relation_rank':rank,'actual_HF2_upper_bound':6-rank,'actual_HF2_equality_or_total_not_claimed':True,
  'source_gate_sha256':H.sha(gp),'source_peer_sha256':H.sha(pp),'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'source_proven_committed_bit_bound':committed,'source_proven_temporary_bit_bound':temporary,'actual_maximum_committed_bits':maxbits,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','known_lower_quadratic_relation_rank','actual_HF2_upper_bound','actual_maximum_committed_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
