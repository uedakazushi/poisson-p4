#!/usr/bin/env python3
"""FORMAL617 directOLD origin-dual fold ONLY; no scalar/upper responses."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_FORMAL_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def enc(a):return [[i,str(x)] for i,x in sorted(a.items())]
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};bits=0;cursor={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63/';r=read(base+'report.json');assert r['all_passed'] and pins[base+'report.json']=='2c00ddfe0a10dad52925b490f16787c31db03c8d4b032a443c43492469841ef0'
 ip=base+'617-direct-old-highestthird-coherent-source-positive-input.json';d=read(ip);assert pins[ip]==r['evidence_sha256'][Path(ip).name] and d['all_pruned_roots_explicit_type_bound_no_missing_defaults']
 p=read(peer);c=read(custody);assert p['all_passed'] and pins[base+'report.json'] in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([r,p,c,pins])<1000 and d['positive_prior_committed_bits']==76 and d['temporary_prior_bits']==229
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir())<d['future_formal_RAW_bound']
 dump('phase.json',{'phase':'authorizedFORMAL617beforeorigin_dual','gate_sha256':pins[base+'report.json'],'all_scalar_higher_values_unassigned':True})
 def ck(x):
  nonlocal bits
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());bits=max(bits,b)
  if b>76:
   dump('guard-frontier.json',{'cursor':cursor,'offending_uncommitted_known_QQ':str(x),'no_cap_raise_or_retry':True});raise RuntimeError('source76guard')
  return x
 def add(a,b,c=Q(1)):
  z=dict(a)
  for i,x in b.items():
   q=ck(z.get(i,Q(0))+ck(c*x))
   if q:z[i]=q
   else:z.pop(i,None)
  return z
 def scale(a,c):
  if not c:return {}
  return {i:ck(x*c) for i,x in a.items()}
 origins={i:Q(x) for i,kind,x in d['literal_boundary_origins_and_membership']};membership={i:kind for i,kind,x in d['literal_boundary_origins_and_membership']};co={};cg={}
 def C(i):
  if i in co:return co[i],cg[i]
  assert i in origins
  return origins[i],{i:Q(1)}
 for i,n in d['literal_C0_program']:
  cursor={'stage':'formal_C0_dual','node':i};op=n[0]
  if op=='c':o=ck(Q(n[1],n[2]));g={}
  elif op=='n':a,u=C(n[1]);o=ck(-a);g=scale(u,Q(-1))
  elif op=='a':a,u=C(n[1]);b,v=C(n[2]);o=ck(a+b);g=add(u,v)
  else:
   assert op=='m';a,u=C(n[1]);b,v=C(n[2]);o=ck(a*b);g=add(scale(u,b),v,a)
  co[i]=o;cg[i]=g
 to={};tg={}
 for i,n in d['typed_source_program']:
  cursor={'stage':'formal_typed_dual','node':i};op=n[0]
  if op=='actual_C0_scalar_function':o,g=C(n[1]);g=dict(g)
  elif op=='c':o=ck(Q(n[1],n[2]));g={}
  elif op=='sum':
   o=Q(0);g={}
   for j,a,b in n[1]:s=Q(a,b);o=ck(o+ck(s*to[j]));g=add(g,tg[j],s)
  else:
   assert op=='product';o=ck(Q(n[1],n[2]));g={}
   for j in n[3]:g=add(scale(g,to[j]),tg[j],o);o=ck(o*to[j])
  to[i]=o;tg[i]=g
 gamma=[{i:Q(x) for i,x in row} for row in d['Gamma_two']];columns={}
 for j,rr in d['pruned_original_old_columns_roots']:
  cursor={'stage':'formal_original_old_column_high3','source':j};rows=dict(rr);gg=[]
  for ch,row in enumerate(gamma):
   g={}
   for i,x in row.items():
    if i not in rows:continue
    tag=rows[i]
    if tag[0]=='constant':v={}
    elif tag[0]=='C0':v=C(tag[1])[1]
    else:assert tag[0]=='typed' and tag[1] in tg;v=tg[tag[1]]
    g=add(g,v,x)
   gg.append(g)
  columns[j]=gg
 sourceC=d['original619_source768_combinations']
 def case(j,ch):
  z={}
  for k,x in sourceC[j]:
   if k<640:
    assert k in columns
    z=add(z,columns[k][ch],Q(x))
  return z
 basis={j:[case(j,ch) for ch in range(2)] for j in [563,609]};answers=[];union=set();nonzero=0
 for rec in d['case_source_recipes617']:
  cursor={'stage':'formal617_actual_source_residual','free_case619':rec['free_case619']};pair=[]
  for ch in range(2):
   z={}
   for f in rec['constant_case619_source_factors']:
    j=f['case619'];op,x=f['scale_operator'];s=Q(x) if op=='literal' else -Q(x);assert op in ['literal','negate_literal'];z=add(z,basis[j][ch] if j in basis else case(j,ch),s)
   D,U=map(int,d['positive_common_two_row_envelopes'][ch])
   for i,x in z.items():assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
   union.update(z);pair.append(z)
  nonzero+=bool(pair[0] or pair[1]);answers.append({'free_case619':rec['free_case619'],'two_formal_directOLD_high3_symbol_coefficients':list(map(enc,pair))})
 assert len(answers)==617
 data={'all617_directOLD_high3_formal_rows':answers,'nonzero_formal_residual_cases':nonzero,'surviving_literal_boundary_symbols':[[i,membership[i]] for i in sorted(union)],'all_scalar_third_values_unassigned':True,'basis563609_subtracted_before_fold':True,'no_basis_upper_zero_fill':True,'not_evaluated_actual_upper_terms':['native','knowncrossparts','b1','A2','M1'],'no_new_physical_gauge_implicit_or_actual_scalar_jet_value':True,'no_higher_jet_necessity_claim':True}
 owned=2*count(d)+count([co,cg,to,tg,columns,basis,answers,data,r,p,c,pins])+3000;assert owned<=d['future_formal_numeric_bound']==40920
 ep='617-formal-direct-old-highestthird-symbol-residuals.json';dump(ep,data)
 report={'id':STEM,'all_passed':True,'FORMAL_ONLY_no_actual_scalar_thirds_or_upper_values':True,'formal_residual_cases_nonzero':nonzero,'surviving_literal_symbol_count':len(union),'actual_known_QQ_bits':bits,'owned_numeric_bound':owned,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({'all_passed':True,'nonzero_formal_cases':nonzero,'symbols':len(union),'sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
