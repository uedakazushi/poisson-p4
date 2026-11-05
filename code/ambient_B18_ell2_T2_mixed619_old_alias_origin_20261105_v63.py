#!/usr/bin/env python3
"""Stream selected arithmetic aliases to cached boundaries, not fresh jets."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_OLD_ALIAS_ORIGIN_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 prof=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_CALLBACK_PROFILE_20261105_v63/literal-oldB-native-X2-callback-and-firstcarry-profile.json')
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 known=read(mug,['all_exact_cached_leaf_bindings_origin_and_ds'])['all_exact_cached_leaf_bindings_origin_and_ds'];K={r['original_node']:r for r in known}
 leaves=prof['reachable_actual_C0_scalar_labels'];missing=sorted(set(leaves)-set(K));assert len(missing)==12
 cp='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';inputs[cp]={'path':cp,'sha256':sha(ROOT/cp)}
 text=(ROOT/cp).read_text();assert len(text.encode())<12*1024*1024
 dec=json.JSONDecoder();field=text.index('"DAG_scalar_nodes"');begin=text.index('[',field)+1
 nodes={};front=set(missing);passes=0
 while front:
  wanted=front-set(K)-set(nodes);front=set()
  if not wanted:break
  passes+=1;assert passes<=24
  pos=begin
  for i in range(max(wanted)+1):
   while text[pos].isspace() or text[pos]==',':pos+=1
   n,pos=dec.raw_decode(text,pos)
   if i in wanted:nodes[i]=n
   del n
  assert count(nodes)<4000 and len(nodes)<1000
  for i in wanted:
   n=nodes[i];op=n[0]
   if op in ('a','m'):front.update(n[1:3])
   elif op in ('n','i'):front.add(n[1])
   else:assert op=='c',('unbound_nonarithmetic_boundary',i,n)
 # Only a small selected closure is decoded/retained, never the whole DAG.
 usedboundaries=sorted({j for n in nodes.values() for j in (n[1:3] if n[0] in ('a','m') else n[1:2] if n[0] in ('n','i') else []) if j in K})
 bounds={i:max(abs(sp(K[i]['origin_and_ds']).get(0,Q(0)).numerator).bit_length(),sp(K[i]['origin_and_ds']).get(0,Q(0)).denominator.bit_length()) for i in K}
 zero={i:not sp(K[i]['origin_and_ds']).get(0,Q(0)) for i in K};prior=max(bounds.values())
 for i,n in sorted(nodes.items()):
  op=n[0]
  if op=='c':b=max(abs(n[1]).bit_length(),abs(n[2]).bit_length());z=n[1]==0
  elif op=='n':b=bounds[n[1]];z=zero[n[1]]
  elif op=='i':assert not zero[n[1]];b=bounds[n[1]];z=False
  else:
   a,b0=n[1:3]
   if op=='m':z=zero[a] or zero[b0];b=1 if z else bounds[a]+bounds[b0]
   else:z=zero[a] and zero[b0];b=bounds[b0] if zero[a] else bounds[a] if zero[b0] else bounds[a]+bounds[b0]+1
  bounds[i]=b;zero[i]=z;prior=max(prior,b)
 assert prior<=384
 vals={i:sp(rec['origin_and_ds']).get(0,Q(0)) for i,rec in K.items()};actual=0
 for i,n in sorted(nodes.items()):
  op=n[0]
  if op=='c':v=Q(n[1],n[2])
  elif op=='n':v=-vals[n[1]]
  elif op=='i':assert vals[n[1]];v=1/vals[n[1]]
  elif op=='a':v=vals[n[1]]+vals[n[2]]
  else:v=vals[n[1]]*vals[n[2]]
  b=max(abs(v.numerator).bit_length(),v.denominator.bit_length());assert b<=bounds[i];actual=max(actual,b);vals[i]=v
 assert all(vals[i]==0 for i in missing)
 assert all(sp(K[i]['origin_and_ds']).get(0,Q(0))==0 for i in leaves if i in K)
 owned=count([prof,known,nodes,bounds,zero,vals])+10000;assert owned<60000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'twelve-arithmetic-alias-closure-and-zero-origin-bindings.json'
 ev={'missing12_literal_nodes':missing,'selected_arithmetic_nodes':[[i,n] for i,n in sorted(nodes.items())],
 'cached_boundary_records_used':[K[i] for i in usedboundaries],'all62_literal_leaf_origins_zero':True,'all12_alias_origins':[[i,str(vals[i])] for i in missing],
 'stream_passes':passes,'selected_numeric_closure':count(nodes),'a_priori_scalar_origin_bits':prior,'actual_scalar_origin_bits':actual,'whole_owned_bound':owned,
 'no_fresh_physical_or_implicit_function_coefficient':True,'no_regular_nonunit_quotient_redivided':True,
 'scope':'Only cachedboundaryconstantsandliteralunitarithmeticat0. Truehonestregularfunctionsremainopaqueatcachedboundaries. Aliasfirsts/secondsandall619valuesnotcomputed.',
 'next_AD_economy':'Everyoneof62typedleaforiginsiszero. Aproductoftwoormorepositiveorderchildrenhaszerolinearsensitivitytohonestsecondcoefficients;firstcrossproductsstillremain. Exactnewweightedgradientfoldneedsitsownsource/resourcegate.'}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all62_leaf_origins_zero':True,'selected_arithmetic_node_count':len(nodes),'source_only_cached_origin_arithmetic':True,'prior_bits':prior,'actual_bits':actual,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'nodes':len(nodes),'passes':passes,'prior':prior,'actual':actual,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
