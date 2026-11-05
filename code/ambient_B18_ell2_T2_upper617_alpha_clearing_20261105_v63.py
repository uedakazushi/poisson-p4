#!/usr/bin/env python3
"""Prepared full918 alpha and836 G3 dots only; no rho/RHS/source data."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_INPUT_GATE_20261105_v63/'
GATE_SHA='6e68e307280750ea713beebb73fa1c3979cc7d122b6389c899e11484c9c804c9'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {int(i):F(x) for i,x in v if F(x)}
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def rec(v):return [[i,str(x)] for i,x in sorted(v.items()) if x]
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};cursor={};pending={};accumulator=F(0);maximum=0;alpha=[];ljall=[];Gzero=[]
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<2_000_000
 g=read(GATE+'report.json');assert g['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'actual-four-row-to-full918-alpha-clearing-pruned-input.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([g,p,c,pins])<1000 and g['ready_with_unchanged384'] and g['all836_G3_check_positive_prior_bits']==137 and g['temporary_prior_bits']==412
 assert g['future_value_numeric_bound']==24230 and count(d)==2911 and d['full836_G3_annihilation_checks_required']
 Pi=list(map(parse,d['PiY3_from102_columns']));J=list(map(parse,d['J25_full918_columns']));FF=list(map(parse,d['F25_columns27']));G3=list(map(parse,d['all418_G3_full918_columns']));pre=list(map(parse,d['pre_two_full918_rows']))
 assert [len(Pi),len(J),len(FF),len(G3),len(pre)]==[102,25,27,418,2]
 env=[{i:(int(D),int(U)) for i,D,U in row} for row in d['positive_two_full918_alpha_envelopes']]
 parsed=count([Pi,J,FF,G3,pre,env]);cap=24230
 #3672 records for two dense sparse-format rows(index,value), twice for
 # serialization;836 literal zero checks twice;all source arrays plus2k misc.
 before=2*count(d)+parsed+2*3672+2*836+50+2000;assert before<cap
 dump('phase.json',{'phase':'source_peer_custody_and_prior_PASS_before_signed918clearing','source_gate_sha256':GATE_SHA,'positive_committed_prior137':137,'temporary_prior412':412,'hard_guards':[384,1153],'whole_numeric_cap':cap,'explicit_prearith_live_upper':before,'no_rho_RHSZ2_b1_or_upper_query':True})
 def stop(reason,x):
  state={'reason':reason,'cursor':dict(cursor),'pending_operation':dict(pending),'partial_accumulator_before_operation':str(accumulator),'uncommitted_fraction':str(x),'completed_alpha_prefix':[rec(v) for v in alpha],'completed_LJ25':[[str(x) for x in v] for v in ljall],'completed_G3_literal_dots':Gzero,'input_sha256':pins,'no_retry_or_cap_raise':True};assert 2*count(d)+parsed+3*count(state)+1000<cap;dump('guard-frontier.json',state);raise RuntimeError(reason)
 def ck(x):
  nonlocal maximum
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maximum=max(maximum,b)
  if b>137:stop('source137_before_commit',x)
  return x
 def mul(a,b):
  pending.update(operation='multiply',operand_a=str(a),operand_b=str(b));return ck(a*b)
 def add(a,b):
  pending.update(operation='add',operand_a=str(a),operand_b=str(b));return ck(a+b)
 def dot(a,b):
  nonlocal accumulator
  q=F(0)
  for i,x in b.items():cursor['dot_source_slot918']=i;accumulator=q;q=add(q,mul(a.get(i,F(0)),x))
  return q
 for t in range(2):
  lj=[];ljall.append(lj)
  for b,col in enumerate(J):
   cursor={'phase':'actual_preJ25','channel':t,'basis25':b};q=dot(pre[t],col);D,U=map(int,d['positive25_LJ_envelopes'][t][b]);assert D%q.denominator==0 and abs(q.numerator)*D<=U*q.denominator;lj.append(q)
  a={};alpha.append(a)
  for i in range(102):
   for var in range(9):
    slot=9*i+var;q=pre[t].get(slot,F(0))
    for target,pi in Pi[i].items():
     for b,x in FF[9*target+var].items():
      cursor={'phase':'actual_all918_fixed25_clearing','channel':t,'slot918':slot,'target3':target,'basis25':b};accumulator=q;term=mul(pi,x);term=mul(lj[b],term);q=add(q,-term)
    D,U=env[t].get(slot,(1,0));assert D%q.denominator==0 and abs(q.numerator)*D<=U*q.denominator
    if q:a[slot]=q
 for t in range(2):
  dots=[];Gzero.append(dots)
  for j,col in enumerate(G3):
   cursor={'phase':'all836_actual_G3_annihilation_dots','channel':t,'G3_source418':j};q=dot(alpha[t],col)
   if q:stop('actual_Kfull_G3_annihilation_identity_failed',q)
   dots.append('0')
 ep='actual-two-full918-alpha-covectors-and836-G3-identities.json';ev={'alpha_two_by918_sparse_rows':[rec(v) for v in alpha],'actual_preJ25_two_by25':[[str(x) for x in v] for v in ljall],'all836_G3_literal_dots':Gzero,'all1836_clearing_slots_and_source_envelopes_checked':True,'all_Kfull_G3linear_image_annihilations_exact0':True,'source_slot_convention':'slot=9*F2source102+var9,zA0..6,rA,ell','formula':'alpha=pre-(preJ25)F25PiY3;Krellrow→r,Krrrow→ell;fullPsiincludedinKfull','scope':'ONLYactualalpha918andknownG3dots.NOrho/rightsolve/RHSZ2/b1/upper617orHF/length.'}
 owned=2*count(d)+parsed+count([alpha,ljall,Gzero,ev,g,p,c,pins,cursor])+1000;assert owned<cap;dump(ep,ev)
 report={'id':STEM,'all_passed':True,'actual_two_full918_alpha_and_G3_checks_only':True,'alpha_supports':list(map(len,alpha)),'all836_G3dots_exact0':True,'source_positive_prior137':137,'temporary_prior412':412,'actual_maximum_bits':maximum,'actual_owned_upper':owned,'whole_numeric_cap':cap,'RAW_output_cap':2000000,'no_rho_RHSZ2_b1_or_upper_HF_value':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','alpha_supports','all836_G3dots_exact0','actual_maximum_bits','actual_owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
