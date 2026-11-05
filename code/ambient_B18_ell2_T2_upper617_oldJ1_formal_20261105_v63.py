#!/usr/bin/env python3
"""FORMAL J0/Jr/Jell of directOLD only, no higher-coefficient responses."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};bits=0;cursor={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/';r=read(base+'report.json');assert r['all_passed'] and pins[base+'report.json']=='ff4e80e5034da41bb93fe4afd9dfd83970cc7edf0d912afbb45dfeb433fb9110'
 ip=base+'617-direct-old-known-first-formal-Jacobian-positive-input.json';d=read(ip);assert pins[ip]==r['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and pins[base+'report.json'] in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([r,p,c,pins])<1000 and d['source_committed_prior_bits']==82 and d['temporary_prior_bits']==247
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(x.stat().st_size for x in out.iterdir())<d['future_formal_RAW_bound']
 dump('phase.json',{'phase':'authorizedFORMAL_OLD_J1_before_signed_fold','higher_source_values_all_unassigned':True,'gate_sha256':pins[base+'report.json']})
 def ck(x):
  nonlocal bits
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());bits=max(bits,b)
  if b>82:
   dump('guard-frontier.json',{'cursor':cursor,'offending_uncommitted_known_QQ':str(x),'no_cap_raise_or_retry':True});raise RuntimeError('source82guard')
  return x
 Z=(Q(0),)*3;one=(Q(1),Q(0),Q(0))
 def va(a,b):return tuple(ck(x+y) for x,y in zip(a,b,strict=True))
 def vs(a,c):return tuple(ck(x*c) for x in a)
 def vm(a,b):return (ck(a[0]*b[0]),ck(ck(a[0]*b[1])+ck(a[1]*b[0])),ck(ck(a[0]*b[2])+ck(a[2]*b[0])))
 def add(a,b,c=one):
  z=dict(a)
  for i,x in b.items():
   q=va(z.get(i,Z),vm(c,x))
   if any(q):z[i]=q
   else:z.pop(i,None)
  return z
 def scale(a,c):
  z={}
  for i,x in a.items():
   q=vm(c,x)
   if any(q):z[i]=q
  return z
 fv={i:tuple(map(Q,v)) for i,kind,v in d['literal_boundary_firsts_point_r_ell_and_membership']};membership={i:kind for i,kind,v in d['literal_boundary_firsts_point_r_ell_and_membership']};co={};cg={}
 def C(i):
  if i in co:return co[i],cg[i]
  assert i in fv
  return fv[i],{i:one}
 for i,n in d['literal_C0_program']:
  cursor={'stage':'formal_C0_Jfirst','node':i};op=n[0]
  if op=='c':o=(ck(Q(n[1],n[2])),Q(0),Q(0));g={}
  elif op=='n':a,u=C(n[1]);o=vs(a,Q(-1));g=scale(u,(Q(-1),Q(0),Q(0)))
  elif op=='a':a,u=C(n[1]);b,v=C(n[2]);o=va(a,b);g=add(u,v)
  else:assert op=='m';a,u=C(n[1]);b,v=C(n[2]);o=vm(a,b);g=add(scale(u,b),v,a)
  co[i]=o;cg[i]=g
 to={};tg={}
 for i,n in d['typed_source_program']:
  cursor={'stage':'formal_typed_Jfirst','node':i};op=n[0]
  if op=='actual_C0_scalar_function':o,g=C(n[1]);g=dict(g)
  elif op=='c':o=(ck(Q(n[1],n[2])),Q(0),Q(0));g={}
  elif op=='sum':
   o=Z;g={}
   for j,a,b in n[1]:s=Q(a,b);o=va(o,vs(to[j],s));g=add(g,tg[j],(s,Q(0),Q(0)))
  else:
   assert op=='product';o=(ck(Q(n[1],n[2])),Q(0),Q(0));g={}
   for j in n[3]:g=add(scale(g,to[j]),tg[j],o);o=vm(o,to[j])
  to[i]=o;tg[i]=g
 gamma=[{i:Q(x) for i,x in row} for row in d['Gamma_two']];columns={}
 for j,rr in d['pruned_original_old_columns_roots']:
  cursor={'stage':'formal_OLD_Jfirst_column','column':j};roots=dict(rr);pair=[]
  for row in gamma:
   z={}
   for i,x in row.items():
    if i not in roots:continue
    tag=roots[i]
    if tag[0]=='constant':v={}
    elif tag[0]=='C0':v=C(tag[1])[1]
    else:assert tag[0]=='typed' and tag[1] in tg;v=tg[tag[1]]
    z=add(z,v,(x,Q(0),Q(0)))
   pair.append(z)
  columns[j]=pair
 def case(j,ch):
  z={}
  for k,x in d['original619_source768_combinations'][j]:
   if k<640:z=add(z,columns[k][ch],(Q(x),Q(0),Q(0)))
  return z
 basis={j:[case(j,ch) for ch in range(2)] for j in [563,609]};answers=[];union=[set(),set(),set()];nonzero=[0,0,0]
 for rec in d['case_source_recipes617']:
  cursor={'stage':'formal617_Jfirst_source_subtraction','free_case619':rec['free_case619']};pair=[];flags=[False]*3
  for ch in range(2):
   z={}
   for f in rec['constant_case619_source_factors']:
    j=f['case619'];op,x=f['scale_operator'];s=Q(x) if op=='literal' else -Q(x);assert op in ['literal','negate_literal'];z=add(z,basis[j][ch] if j in basis else case(j,ch),(s,Q(0),Q(0)))
   for i,v in z.items():
    for k,x in enumerate(v):
     D,U=map(int,d['positive_common_two_by3_origin_Jr_Jell_envelopes'][ch][k]);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
     if x:union[k].add(i);flags[k]=True
   pair.append([[i,list(map(str,v))] for i,v in sorted(z.items())])
  for k,flag in enumerate(flags):nonzero[k]+=flag
  answers.append({'free_case619':rec['free_case619'],'two_formal_J0_Jr_Jell_symbol_rows':pair})
 data={'all617_formal_directOLD_J0_Jr_Jell':answers,'ordinary_sensitivity_order':['third_r2ell','second_rell','second_rr'],'all_higher_coefficients_unassigned':True,'nonzero_formal_cases_by_component':nonzero,'literal_symbols_by_component':[[[i,membership[i]] for i in sorted(v)] for v in union],'not_evaluated_actual_upper_terms':['first_only_cubic','boundary_second_values','boundary_third_values','native','b1','A2','M1'],'no_new_physical_gauge_implicit_or_actual_scalar_jet_value':True,'no_higher_jet_necessity_claim':True}
 owned=2*count(d)+count([fv,co,cg,to,tg,columns,basis,answers,data,r,p,c,pins])+4000;assert owned<=d['future_formal_numeric_bound']==49214
 ep='617-formal-direct-old-third-and-second-symbol-sensitivities.json';dump(ep,data)
 report={'id':STEM,'all_passed':True,'FORMAL_ONLY_no_actual_second_third_or_upper_values':True,'nonzero_formal_cases_by_component':nonzero,'surviving_symbol_counts_by_component':list(map(len,union)),'actual_known_QQ_bits':bits,'owned_numeric_bound':owned,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({'all_passed':True,'nonzero_cases':nonzero,'symbols':list(map(len,union)),'sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
