#!/usr/bin/env python3
"""Prepared first-only witness scan, no new second/native/gauge coefficient."""
import argparse,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_20261105_v63'
def main(out,source_peer,formal_peer):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_INPUT_20261105_v63/report.json';g=json.loads(gp.read_text());pp=ROOT/source_peer;peer=json.loads(pp.read_text());assert g['all_passed'] and peer['all_passed'] and sha(gp) in json.dumps(peer)
 fp=ROOT/formal_peer;fpeer=json.loads(fp.read_text());assert fpeer['all_passed'] and '83f5ab73f6021da914ce14206897e09e7e40178425a009cf4025229f6ccbd761' in json.dumps(fpeer)
 ip=gp.parent/'pruned548-first-only-cache-source-and-positive-bound.json';assert sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 assert d['no_actual548_lower_values'] and d['source_prior_ready'] and d['positive_committed_bit_bound']==60 and d['positive_temporary_bit_bound']==181
 for h in g['input_sha256'].values():assert sha(ROOT/h['path'])==h['sha256']
 assert count([g,peer,fpeer])<1000 and 2*count(d)+18*len(d['typed_closure_nodes'])+8*548+6000<g['future_numeric_bound']==40368
 (out/'phase.json').write_text(json.dumps({'phase':'AUTHORIZED_first_only548_before_any_images','source_gate_sha256':sha(gp),'source_peer_sha256':sha(pp),'formal619_postpeer_sha256':sha(fp),'no_fresh_function_coefficients':True,'arbitrary_second_extension_only_after_complete_symbolic_cancellation':True},separators=(',',':'))+'\n')
 cap=60;maximum=0
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=cap
  return x
 pairs=d['ordinary_pair_order8'];first={i:sp(v) for i,v in d['actual_C0_leaf_first10']};ops=dict(d['typed_closure_nodes']);jets={}
 def jm(a,b):
  v=[ck(a[0]*b[0])]
  for k in range(1,10):v.append(ck(ck(a[0]*b[k])+ck(b[0]*a[k])))
  for j,(u,w) in enumerate(pairs):
   q=ck(ck(a[0]*b[10+j])+ck(b[0]*a[10+j]));q=ck(q+ck(a[u+1]*b[w+1]))
   if u!=w:q=ck(q+ck(a[w+1]*b[u+1]))
   v.append(q)
  return v
 for i,n in sorted(ops.items()):
  op=n[0]
  if op=='actual_C0_scalar_function':v=[first[n[1]].get(j,Q(0)) for j in range(10)]+[Q(0)]*8
  elif op=='c':v=[Q(n[1],n[2])]+[Q(0)]*17
  elif op=='sum':
   v=[Q(0)]*18
   for j,a,b in n[1]:v=[ck(x+ck(Q(a,b)*y)) for x,y in zip(v,jets[j],strict=True)]
  else:
   assert op=='product';v=[Q(n[1],n[2])]+[Q(0)]*17
   for child in n[3]:v=jm(v,jets[child])
  jets[i]=v
 G=list(map(sp,d['Gamma4_source56_times8']));B={j:dict(v) for j,v in d['oldB_source_rows_pruned']};R={p:sp(v) for p,v in d['actual_pivot_r_first_pruned']};E={p:sp(v) for p,v in d['actual_old_pivot_ell_first_pruned']}
 L0=list(map(sp,d['native_L0_rows_pruned']));Le=list(map(sp,d['native_L1ell_rows_pruned']));S0={p:sp(v) for p,v in d['native_pivot_S0_columns']};Se={p:sp(v) for p,v in d['native_pivot_S1ell_columns']}
 for p,j in d['native_pivot_local128_ids']:
  v={}
  for A,b in [(L0,Se[p]),(Le,S0[p])]:
   for col,c in b.items():
    for row,x in A[col].items():
     y=ck(v.get(row,Q(0))+ck(c*x))
     if y:v[row]=y
     else:v.pop(row,None)
  E[p]=v
 computed=[];witness=None
 for rec in d['selected548_source_records']:
  source=rec['original_old640_source'];a1=sp(rec['actual_a1_27_times9']);assert all(k%9<7 for k in a1)
  values=[]
  for target,row in enumerate(G):
   x=Q(0)
   for flat,c in row.items():
    source_row,pair=divmod(flat,8);root=B[source].get(source_row)
    if root is not None:x=ck(x+ck(c*jets[root][10+pair]))
    u,v=pairs[pair]
    if u<7:
     for p in d['needed_pivot_indices']:
      a=a1.get(9*p+u,Q(0))
      if a:
       firstv=R[p].get(source_row,Q(0)) if v==7 else E[p].get(source_row,Q(0))
       x=ck(x-ck(c*ck(firstv*a)))
   D,U=map(int,d['all548_four_images_global_D_U'][target]);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
   values.append(x)
  computed.append({'case619':rec['case619'],'original_old640_source':source,'actual_lower4_values':list(map(str,values))})
  if any(values):witness={**computed[-1],'actual_a1_27_times9':rec['actual_a1_27_times9'],'first_nonzero_lower_channel':next(i for i,x in enumerate(values) if x)};break
 assert count([d,jets,G,B,R,E,L0,Le,S0,Se,computed,witness])+5000<g['future_numeric_bound']
 path=out/'actual548-first-only-lower-prefix-and-witness-if-any.json';path.write_text(json.dumps({'computed_source_prefix':computed,'selected_family_count':548,'computed_prefix_count':len(computed),'nonzero_witness':witness,
 'all548_zero':witness is None and len(computed)==548,'actual_function_second_parts_not_set_zero':True,'proof':'Thecompleteformaloldsecondweightvanishesforeachselectedcase/all4channels. Hencearbitrary0secondextensionoftheexactcachedfirstscomputesthetrueweightedX2image. FirstcrossproductsandactualindividualnativepivotLeS0areincluded.',
 'no_nativeX2_or_higher_function_response':True,'scope':'OnenonzerolowervectorwouldgiveatruequadraticinitialrelationandHF2<=5.Allzeroonlycloses548lower4subfamily;upperchannelsandremaining71cyclesstayopen.Norankorlengthclaimed.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'computed_prefix_count':len(computed),'nonzero_lower_witness_found':witness is not None,'all548_lower4_zero':witness is None and len(computed)==548,'no_new_physical_gauge_or_second_response':True,
 'source_prior_bits':60,'actual_bits':maximum,'input_sha256':{str(gp.relative_to(ROOT)):sha(gp),str(pp.relative_to(ROOT)):sha(pp),str(fp.relative_to(ROOT)):sha(fp),str(ip.relative_to(ROOT)):sha(ip)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<g['future_RAW_bound'];print(json.dumps({'all_passed':True,'prefix':len(computed),'witness':witness,'actual_bits':maximum,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--source-peer',required=True);p.add_argument('--formal-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.source_peer,a.formal_peer)
