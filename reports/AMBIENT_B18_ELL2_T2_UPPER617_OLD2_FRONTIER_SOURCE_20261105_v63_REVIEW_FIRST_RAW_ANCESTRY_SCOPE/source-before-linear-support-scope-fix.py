#!/usr/bin/env python3
"""Metadata-only original C0 constructor/cache-label audit; no coefficients."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63';STEM=A+'_REVIEW'
def sha_bytes(x):return hashlib.sha256(x).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple,set)):return sum(map(count,x))
 return int(isinstance(x,(int,float)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;cache={}
 def read(p):
  nonlocal raw
  if p in cache:return cache[p]
  b=(ROOT/p).read_bytes();pins[p]=sha_bytes(b);raw+=2*len(b);assert raw<32*1024*1024
  z=json.loads(b);cache[p]=z;return z
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 r=read(f'reports/{A}/report.json');assert pins[f'reports/{A}/report.json']=='c158d5ae6733823b25dc828c528ad00bc61ba5276ae0cb02e39099ece4e59731'
 ep='literal12-opaque-function-and24-plane-second-frontier-with-existing-rr-labels.json';d=read(f'reports/{A}/{ep}');assert pins[f'reports/{A}/{ep}']==r['evidence_sha256'][ep]
 primary=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json')
 peer=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json');assert peer['all_passed']
 inv=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json')
 outer=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_PIVOT4_SECOND_SOURCE_20261105_v63/literal-four-OLDpivot-Gamma-plane2-callbacks-and-frontier.json')
 supports={}
 for case,chs in primary['all617_combined_OLD2_independent_symbol_coefficients']:
  for ch,v in enumerate(chs):
   for pair,f,q in v:assert q!='0';supports.setdefault(f,set()).add(pair)
 assert [[f,sorted(v)] for f,v in sorted(supports.items())]==d['all42_surviving_formal_symbol_pair_support']
 oldg=[330164,330289,330524,330594];oldh=[330184,330309,330544,330614];newg=list(range(566400,566404));symbols=sorted(oldg+oldh+newg)
 assert symbols==d['opaque_remaining_functions12'] and all(supports[i]=={8,9} for i in symbols)
 assert inv['outer_four_pivot_source']==outer and outer['original4_old_pivots']==[604,613,622,631]
 assert all(j not in outer['original4_old_pivots'] for v in inv['original619_c0_source_combinations'] for j,x in v)
 # Independent metadata closure of D roots: the12 surviving opaque labels have
 # no D-side function occurrence; they come from the outer pivot carry only.
 ops=dict(inv['typed_program']);todo=[node for j,col in inv['D_old_root_columns'] for row,node in col];seen=set();Dleaves=set()
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i);n=ops[i]
  if n[0]=='actual_C0_scalar_function':Dleaves.add(n[1])
  elif n[0]=='sum':todo.extend(v[0] for v in n[1])
  else:assert n[0]=='product';todo.extend(n[3])
 assert not(Dleaves&set(symbols))
 # Stream literal source once, stopping all genuine r/g quotient constructors.
 # This is an ancestry check, never an Artin division or a response evaluation.
 cp='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';q=ROOT/cp;b=q.read_bytes();pins[cp]=sha_bytes(b);raw+=len(b);assert pins[cp]==r['input_sha256'][cp];del b
 text=q.read_text();raw+=q.stat().st_size;assert raw+1_000_000<32*1024*1024
 wanted={i for i,n in d['literal_original_C0_arithmetic_and_regular_quotient_nodes']}|{327554}|{i for i,n in outer['literal_C0_program']}
 dec=json.JSONDecoder();pos=text.index('[',text.index('"DAG_scalar_nodes"'))+1;found={}
 dump('phase.json',{'phase':'metadata_owner_and_wholeRAW_before_original_C0_stream_PASS','RAW_upper_before_output':raw+1_000_000,'no_function_response_values':True})
 for i in range(max(wanted)+1):
  while text[pos].isspace() or text[pos]==',':pos+=1
  n,pos=dec.raw_decode(text,pos)
  if i in wanted:found[i]=n
 del text
 assert set(found)==wanted
 for i,n in d['literal_original_C0_arithmetic_and_regular_quotient_nodes']+d['literal_higher_arithmetic_boundaries']+outer['literal_C0_program']:assert found[i]==n,(i,n,found[i])
 boundary=327554;assert found[boundary]==['r',327553,27948,'f63:z1']
 closure=set();todo=list(symbols)
 while todo:
  i=todo.pop()
  if i==boundary or i in closure:continue
  closure.add(i);n=found[i]
  if n[0] in ['r','g']:continue
  assert n[0] in ['a','m','n','c'];todo.extend(n[1:3] if n[0] in ['a','m'] else n[1:2] if n[0]=='n' else [])
 assert closure=={i for i,n in d['literal_original_C0_arithmetic_and_regular_quotient_nodes']} and len(closure)==24
 rawr=[330168,330293,330528,330598]
 for k,(g,h,n,r0) in enumerate(zip(oldg,oldh,newg,rawr,strict=True)):
  assert found[g][0]=='g' and found[g][2]==299650 and found[n][0]=='g' and found[n][2]==505938
  assert found[r0][0]=='r' and found[r0][2]==27948
  a0,a1=found[h][1:];assert found[h][0]=='a' and a0==r0 and found[a1][0]=='n'
  prod=found[a1][1];assert found[prod]==['m',327554,g]
 # Current negative-boundary source uses exactly this C0 file. Verify the
 # original OLD-minus-NEW and typed h/g routes, without expanding bound aliases.
 negp='reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/report.json';neg=read(negp)
 assert neg['all_passed'] and neg['actual_chain_identity_UP_QB'] and neg['input_hashes']['C0']=={'path':cp,'sha256':pins[cp]}
 typed=dict(outer['typed_program']);literal=dict(outer['literal_C0_program'])
 for k,(source,col) in enumerate(outer['four_source_root_columns']):
  assert [v[0] for v in col]==[4,28,55]
  assert typed[col[0][1][1]]==['actual_C0_scalar_function',oldh[k]]
  assert typed[col[1][1][1]]==['actual_C0_scalar_function',oldg[k]]
  root=col[2][1][1];a,b=literal[root][1:];assert literal[root][0]=='a' and a==oldg[k] and literal[b]==['n',newg[k]]
 # Existing rr cache presence only: take labels, never copy/contract coefficients.
 up='reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json';u=read(up)
 labels=[v[0] for v in u['actual_eight_regular_g_second_values']];assert labels==oldg+newg;del u;cache.pop(up)
 for p0 in ['reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json']:assert read(p0)['all_passed']
 remaining=[[i,p0] for i in symbols for p0 in [8,9] if not(p0==8 and i in labels)]
 assert remaining==d['remaining_individual_scalar_frontier_after_rr_cache_presence'] and len(remaining)==16
 # The source report's literal file hashes are checked through this read cache;
 # any file not needed for parsing is hashed exactly once (no source replay).
 for path,h in r['input_sha256'].items():
  if path not in pins:
   bb=(ROOT/path).read_bytes();raw+=len(bb);pins[path]=sha_bytes(bb)
  assert pins[path]==h,path
 owned=count([cache,d,primary,inv,outer,supports,found,Dleaves,seen,closure,remaining,labels,pins])+10000
 assert owned<60000 and raw+1_000_000<32*1024*1024
 ep2='independent24-C0-constructors-and-eight-rr-cache-labels-only.json';e={'all617_formal_pair_support_and12_opaqueIDs_verified':True,
  'all24_original_C0_constructors_and_boundary327554_verified':True,'exact_hrest_rawr_minus_hf63_oldg_identity':True,
  'rawr_denominator27948_oldg299650_newg505938_preserved':True,'current_negative_boundary_same_C0_and_signed_oldminusnew_verified':True,
  'all12_opaque_labels_absent_D_function_ancestry_only_outercarry':True,'all8_rr_cache_labels_present_no_values_copied':True,
  'remaining16_is_inventory_not_newjet_necessity':True,'no_gradient_or_coefficient_response_evaluated':True,
  'whole_owned_numeric_upper':owned,'whole_RAW_upper':raw+1_000_000}
 dump(ep2,e);rep={'id':STEM,'all_passed':True,'metadata_source_only':True,'input_sha256':pins,'evidence_sha256':{ep2:sha_bytes((out/ep2).read_bytes())},
  'code_sha256':sha_bytes(Path(__file__).read_bytes()),'note_sha256':sha_bytes((ROOT/'research/notes'/f'{STEM}.md').read_bytes()),'whole_numeric_upper':owned,'whole_RAW_upper':raw+1_000_000,'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha_bytes((out/'report.json').read_bytes()),'owned':owned,'RAW':raw+1_000_000,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
