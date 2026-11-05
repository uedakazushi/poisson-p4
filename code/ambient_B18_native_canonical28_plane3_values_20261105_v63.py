#!/usr/bin/env python3
"""Prepared ONE C28 cache-only34 coefficient worker; requires source peer.

UNRUN until root authorizes immutable code/custody. No new physical source or
gauge solve, N adapter, native pivot coefficient or comparison response.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_VALUES_20261105_v63'
GATE='reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_SOURCE_V2_20261105_v63/'
GATE_SHA='50565ea2320e80f30e1416b9978594449d985ef51aba58849eca65dc201cdc3a'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,str):
  try:F(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(x,(int,float,F)))
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 pins={};cursor={};maximum=0;ops=0;values=[];powers={};routes={};kval=[];nu=[];current={};accumulator={}
 def read(p):q=ROOT/p;pins[p]=sha(q);return json.loads(q.read_text())
 def dump(n,x):
  (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
  assert sum(p.stat().st_size for p in out.iterdir())<1_000_000
 r=read(GATE+'report.json');assert r['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'native-P7-C28-canonical-third-source-and-positive34-response-input.json';d=read(ip);assert pins[ip]==r['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody)
 assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS']
 assert c['worker_sha256']==sha(Path(__file__)) and c['gate_sha256']==GATE_SHA
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert r['positive_committed_bits']==d['positive_committed_bits']==37 and d['positive_temporary_bits']==112
 assert d['prospective_value_whole_numeric_bound']==r['future_value_numeric_upper']==14403
 cap=14403;prior=37
 assert count(d)==1899 and count([r,p,c,pins])<1000
 assert c['explicit_success_and_guard_numeric_upper']<cap and c['whole_RAW_upper']<32*1024*1024
 dump('phase-before-values.json',{'phase':'literal_source_peer_custody_owner_and_prior_PASS','gate_sha256':GATE_SHA,'source37_temporary112':True,'whole_numeric_cap':cap,'no_new_physical_or_gauge_solve':True})
 def owned(extra=()):return 2*count(d)+count([r,p,c,pins,kval,nu,powers,routes,values,current,accumulator,extra])+1000
 def stop(q,why):
  frontier={'reason':why,'cursor':cursor,'uncommitted_fraction':str(q),'completed32_response_prefix':values,
   'current_term':[[i,str(v)] for i,v in current.items()],'current_accumulator':[[i,str(v)] for i,v in accumulator.items()],
   'power_cache':[[list(k),[[i,str(v)] for i,v in x.items()]] for k,x in powers.items()],
   'input_sha256':pins,'no_retry_or_cap_raise':True}
  assert owned(frontier)<cap;dump('guard-frontier.json',frontier);raise RuntimeError(why)
 def ck(q):
  nonlocal maximum
  maximum=max(maximum,abs(q.numerator).bit_length(),q.denominator.bit_length())
  if maximum>prior:stop(q,'source37 STOP before commit')
  return q
 def plus(a,b):return ck(a+b)
 def times(a,b):
  nonlocal ops
  ops+=1
  if ops>d['positive_propagation_operations_upper']:stop(F(ops),'source operation bound STOP')
  return ck(a*b)
 def product(a,b):
  z={}
  for i,u in a.items():
   for j,v in b.items():
    for k,w in routes.get((i,j),[]):
     q=plus(z.get(k,F(0)),times(times(u,v),w))
     if q:z[k]=q
     else:z.pop(k,None)
  return z
 basis=[tuple(v) for v in d['P7_basis_r_ell']];index={m:i for i,m in enumerate(basis)}
 assert basis==[(0,0),(1,0),(0,1),(2,0),(1,1),(3,0),(2,1)]
 nu=[[ck(F(q)) for q in v] for v in d['actual_nu_two_P7_cached_inputs']]
 kval=[{i:ck(F(q)) for i,q in enumerate(v) if F(q)} for v in d['actual_kappa9_C28_cached_inputs_UNEVALUATED_AS_G8']]
 assert all(0 not in v for v in kval) and all(v[:3]==[F(0)]*3 for v in nu)
 for x in range(28):
  bx,mx=divmod(x,4)
  for y in range(28):
   by,my=divmod(y,4);a,e=basis[bx];b,f=basis[by];collision=mx&my;om=mx^my;route=[]
   if collision==0:
    if (a+b,e+f) in index:route=[(4*index[a+b,e+f]+om,F(1))]
   elif collision in [1,2]:
    for k,(u,v) in enumerate(basis):
     if nu[collision-1][k] and (a+b+u,e+f+v) in index:route.append((4*index[a+b+u,e+f+v]+om,nu[collision-1][k]))
   else:assert collision==3 # νAνB order>=4, zero in P7.
   if route:routes[x,y]=route
 powers.update({(j,1):v for j,v in enumerate(kval)})
 for j,v in enumerate(kval):
  for n in range(2,6):cursor={'phase':'cached_kappa_power','kappa':j,'power':n};powers[j,n]=product(powers[j,n-1],v)
 assert owned()<cap
 for g,polynomial in enumerate(d['G8_formal_kappa9_polynomials_mod_radical6']):
  accumulator={}
  for term,(m,q) in enumerate(polynomial):
   cursor={'phase':'original_exact_G8_polynomial','generator':g,'term':term};current={0:ck(F(q))}
   for j,n in enumerate(m):
    if n:current=product(current,powers[j,n])
   for k,v in current.items():
    q=plus(accumulator.get(k,F(0)),v)
    if q:accumulator[k]=q
    else:accumulator.pop(k,None)
  # Only the requested r²ell four masks are response outputs. Other components
  # of this current algebra element are internal polynomial arithmetic only.
  v=[accumulator.get(24+m,F(0)) for m in range(4)]
  for mask,q in enumerate(v):
   den,u=map(int,d['positive_G8_requested32_D_U'][g][mask]);assert den%q.denominator==0 and abs(q.numerator)*den<=u*q.denominator
  values.append(list(map(str,v)));assert owned()<cap
 monic=[str(v[6]) for v in nu]
 for i,q in enumerate(map(F,monic)):
  den,u=map(int,d['positive_monic_requested2_D_U'][i]);assert den%q.denominator==0 and abs(q.numerator)*den<=u*q.denominator
 ep='actual32-G8-r2ell-and-two-monic-cached-coefficients.json'
 evidence={'source_G8_order':d['source_G8_order'],'ordinary_monomial':'rA^2*ell','free4_basis_mask_order':d['free4_basis_mask_order'],
  'actual_G8_r2ell_eight_by4':values,'actual_nu_A_B_r2ell':monic,
  'all34_positive_source_envelopes_verified':True,'same_original_canonical_unitIFT_branch_and_physical_R0':True,
  'no_new_physical_source_or_gauge_solve':True,'no_N3_Xp2_Xc3_comparison_upper_HF_or_length_query':True}
 owner=owned(evidence);assert owner<cap;dump(ep,evidence)
 report={'id':STEM,'all_passed':True,'only32_G8_r2ell_and2_monic_response_coefficients':True,'actual_maximum_bits':maximum,
  'source_committed37_temporary112':True,'hard_guards':[384,1153],'whole_numeric_cap':cap,'actual_owned_upper':owner,
  'operations':ops,'source_operation_bound':d['positive_propagation_operations_upper'],
  'no_new_source_or_gauge_solve':True,'no_native_Xp_or_upper_HF_length_response':True,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'seconds':time.monotonic()-start}
 dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_maximum_bits','actual_owned_upper','operations','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
