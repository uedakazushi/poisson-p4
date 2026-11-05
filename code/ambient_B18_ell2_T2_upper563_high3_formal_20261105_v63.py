#!/usr/bin/env python3
"""FORMAL origin dual only; never evaluate physical third/upper responses."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_FORMAL_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def enc(a):return [[i,str(x)] for i,x in sorted(a.items())]
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);pins={};bits=0;cursor={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_INPUT_GATE_20261105_v63/';r=read(base+'report.json');assert r['all_passed'] and pins[base+'report.json']=='4fdfae0b28a4a16a3153d4ad27728aa7cb612fb371a1fed32d896153c6293018'
 ip=base+'direct-old-highestthird-origin-dual-source-and-positive-input.json';d=read(ip);assert pins[ip]==r['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and pins[base+'report.json'] in json.dumps(p) and c['all_passed'] and c['worker_sha256']==sha(Path(__file__)) and c['safe_help_literal_PASS']
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([r,p,c,pins])<600 and d['positive_prior_committed_bits']==57 and d['temporary_bound']==172
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir())<d['future_RAW_bound']
 dump('phase.json',{'phase':'authorized_FORMAL_highestthird_origin_dual_before_fold','gate_sha256':pins[base+'report.json'],'no_source_values':True})
 def ck(x):
  nonlocal bits
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());bits=max(bits,b)
  if b>57:
   dump('guard-frontier.json',{'cursor':cursor,'offending_uncommitted_known_QQ':str(x),'no_retry_or_cap_raise':True});raise RuntimeError('source57guard')
  return x
 def add(a,b,c=Q(1)):
  z=dict(a)
  for i,x in b.items():
   q=ck(z.get(i,Q(0))+ck(c*x))
   if q:z[i]=q
   else:z.pop(i,None)
  return z
 def scale(a,c):return {i:ck(x*c) for i,x in a.items() if x*c}
 def within(x,b):D,U=map(int,b);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
 cb={i:v for i,v in d['positive_node_origin_envelopes']};gb={i:{k:q for k,q in v} for i,v in d['positive_node_high3_sensitivity_envelopes']};origin={};dual={};leaf0={i:Q(x) for i,x in d['actual_scalar_leaf_origins']}
 for i,n in d['typed_source_program']:
  cursor={'stage':'formal_origin_dual','node':i}
  if n[0]=='actual_C0_scalar_function':o=leaf0[n[1]];v={n[1]:Q(1)}
  elif n[0]=='c':o=Q(n[1],n[2]);v={}
  elif n[0]=='sum':
   o=Q(0);v={}
   for j,a,b in n[1]:s=Q(a,b);o=ck(o+ck(s*origin[j]));v=add(v,dual[j],s)
  else:
   assert n[0]=='product';o=Q(n[1],n[2]);v={}
   for j in n[3]:v=add(scale(v,origin[j]),dual[j],o);o=ck(o*origin[j])
  origin[i]=o;dual[i]=v;within(o,cb[i]);assert set(v)<=set(gb[i])
  for k,x in v.items():within(x,gb[i][k])
 targets=[{},{}]
 for t in d['original_factor_terms']:
  cursor={'stage':'formal_weighted_highestthird','channel':t['channel'],'source':t['original_source'],'row':t['source_row56']};w=ck(Q(t['coefficient_factors'][0])*Q(t['coefficient_factors'][1]));v=dual[t['root']] if t['root_kind']=='typed' else {t['root']:Q(1)};targets[t['channel']]=add(targets[t['channel']],v,w)
 for row,b in zip(targets,d['positive_two_weighted_high3_envelopes'],strict=True):
  b=dict(b);assert set(row)<=set(b)
  for i,x in row.items():within(x,b[i])
 data={'case619':563,'two_direct_OLD_highestthird_symbol_coefficients':list(map(enc,targets)),'all_symbols_independent_and_unassigned':True,'formal_ordinary_coefficient_r2ell_only':True,'all_node_origin_and_dual_positive_bounds_verified':True,'whole_upper_source_carry_not_evaluated':['native','knownfirstsecondcross','b1','A2','M1'],'no_new_physical_gauge_implicit_or_scalar_jet_values':True,'no_new_higher_source_order_necessity_claim':True}
 owned=2*count(d)+count([origin,dual,targets,data,r,p,c,pins])+1000;assert owned<=d['future_whole_numeric_bound']==7618
 ep='two-formal-direct-old-highestthird-sensitivity-rows.json';dump(ep,data)
 report={'id':STEM,'all_passed':True,'FORMAL_ONLY_no_scalar_or_upper_response_values':True,'nonzero_formal_entries':sum(map(len,targets)),'actual_known_coefficient_bits':bits,'owned_numeric_bound':owned,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({'all_passed':True,'nonzero_formal_entries':report['nonzero_formal_entries'],'sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
