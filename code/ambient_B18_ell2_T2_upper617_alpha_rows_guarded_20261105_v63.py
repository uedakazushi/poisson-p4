#!/usr/bin/env python3
"""Prepared ONE STOP-only four point dual rows,384/1153 unchanged."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_INPUT_GATE_20261105_v63/'
GATE_SHA='6354180e36b18d7d2f14ee45fbfacd10edd48a3efbb258eafeddc7db7068ae06'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def parse(v):return {int(i):F(x) for i,x in v if F(x)}
def rec(v):return [[i,str(x)] for i,x in sorted(v.items()) if x]
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True)
 pins={};cursor={};completed=[];rhs={};current={};acc=F(0);maximum=0
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<2_000_000
 g=read(GATE+'report.json');assert g['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'pruned-four-point-dualrow-STOP-only-input.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([g,p,c,pins])<1000 and g['STOP_only_execution_method_not_promised_completion'] and d['STOP_only_no_completion_bound']
 assert d['hard_committed_bits']==384 and d['hard_temporary_bits']==1153 and d['single_uncommitted_offender_maxbits']==769
 K=[[[F(x) for x in v] for v in ch] for ch in d['K99_two_channels_rr_rell']];norm={i:parse(v) for i,v in d['normalized99_basis']};rep={i:parse(v) for i,v in d['normalized99_selected_representations']};Y=list(map(parse,d['selected99_Y0_columns102']))
 assert len(norm)==len(rep)==len(Y)==len(K)==99 and set(norm)==set(rep) and set(range(102))-set(norm)==set(d['free_Y3_point_axes'])=={1,36,38}
 for i,v in norm.items():assert v.get(i)==1 and all(j>=i for j in v)
 allinput=[x for ch in K for v in ch for x in v]+[x for m in [norm,rep] for v in m.values() for x in v.values()]+[x for v in Y for x in v.values()]
 assert max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in allinput)==g['initial_cached_operand_maximum_bits']<=384
 del allinput
 parsed=count([K,norm,rep,Y]);cap=g['future_numeric_cap'];assert cap==20406 and 2*count(d)+parsed+6000==g['explicit_future_live_bound']==15516
 dump('phase.json',{'phase':'literal_source_custody_and_STOP_only_bounds_PASS_before_signedrows','gate_sha256':GATE_SHA,'committed_guard':384,'temporary_guard':1153,'binary_temporary_bound':769,'completion_not_promised':True,'whole_numeric_cap':cap,'prearith_live_upper':15516,'no_alpha918_G3dots_or_sourceZ2':True})
 def stop(reason,result,op,a,b):
  state={'reason':reason,'operation_cursor':dict(cursor),'operation':op,'operand_a':str(a),'operand_b':str(b),'uncommitted_fraction':str(result),'uncommitted_bits':max(abs(result.numerator).bit_length(),result.denominator.bit_length()),'current_accumulator_before_operation':str(acc),'completed_point_rows':completed,'all_current_rhs_pivots':rec(rhs),'current_backsubstitution_row':rec(current),'fixed_pivot_order':sorted(norm),'input_sha256':pins,'code_sha256':sha(Path(__file__)),'no_retry_or_cap_raise':True}
  assert state['uncommitted_bits']<=769
  assert 2*count(d)+parsed+3*count(state)+1000<cap
  dump('guard-frontier.json',state);dump('report.json',{'id':STEM,'all_passed':False,'guard_stop_before_commit':True,'no_retry_or_cap_raise':True,'complete_row_prefix_count':len(completed),'source_gate_sha256':GATE_SHA,'guard_frontier_sha256':sha(out/'guard-frontier.json'),'input_sha256':pins,'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start});raise RuntimeError(reason)
 def binary(a,b,op):
  nonlocal maximum
  assert max(abs(a.numerator).bit_length(),a.denominator.bit_length(),abs(b.numerator).bit_length(),b.denominator.bit_length())<=384
  # Only one binary operation occurs before its commit guard. No product is
  # passed to an addition/subtraction until this function returns successfully.
  if op=='multiply':q=a*b
  elif op=='add':q=a+b
  elif op=='subtract':q=a-b
  else:raise AssertionError(op)
  bits=max(abs(q.numerator).bit_length(),q.denominator.bit_length())
  if bits>384:stop('384bit_STOP_before_commit',q,op,a,b)
  maximum=max(maximum,bits);return q
 for t in range(2):
  for pair in range(2):
   rhs={};current={}
   for pivot in sorted(rep):
    acc=F(0)
    for offset,(j,x) in enumerate(sorted(rep[pivot].items())):
     cursor={'stage':'fixed_rhs_K_times_rep','channel':t,'pair_rr_rell':pair,'pivot':pivot,'factor_offset':offset,'selected99_source':j,'pending':'multiply'}
     term=binary(K[j][t][pair],x,'multiply');cursor['pending']='add_rhs';nxt=binary(acc,term,'add');acc=nxt
    rhs[pivot]=acc
   for pivot in reversed(sorted(norm)):
    acc=rhs[pivot]
    for offset,(j,x) in enumerate(sorted(norm[pivot].items())):
     if j==pivot:continue
     cursor={'stage':'fixed_transpose_backsubstitution','channel':t,'pair_rr_rell':pair,'pivot':pivot,'target_coordinate102':j,'factor_offset':offset,'pending':'multiply'}
     term=binary(x,current.get(j,F(0)),'multiply');cursor['pending']='subtract';nxt=binary(acc,term,'subtract');acc=nxt
    if acc:current[pivot]=acc
   assert not set(current)&{1,36,38}
   for kind,columns in [('normalized_basis',[(i,v) for i,v in sorted(norm.items())]),('actual99_source',list(enumerate(Y)))]:
    for col,v in columns:
     acc=F(0)
     for offset,(j,x) in enumerate(sorted(v.items())):
      cursor={'stage':'fixed_point_membership','membership':kind,'channel':t,'pair_rr_rell':pair,'source_or_pivot':col,'target_coordinate102':j,'factor_offset':offset,'pending':'multiply'}
      term=binary(current.get(j,F(0)),x,'multiply');cursor['pending']='add';nxt=binary(acc,term,'add');acc=nxt
     expected=rhs[col] if kind=='normalized_basis' else K[col][t][pair]
     if acc!=expected:stop('fixed_point_source_identity_failed',acc,'identity_compare',acc,expected)
   completed.append({'channel':t,'pair_rr_rell':pair,'point_row102':rec(current),'all99_source_and_normalized_memberships':True})
   dump('phase.json',{'phase':'completed_fixed_point_dual_rows','completed_rows':completed,'maximum_committed_bits':maximum,'no_alpha918_G3dots_or_sourceZ2':True})
 ep='actual-four-specific-B99-point-dual-rows-STOP-only384.json';ev={'actual_four_rows_channel_pair':completed,'source_equations':'s_t,pairYstar0=K[:,t,pair],s_free1,36,38=0','ordinary_pair_mapping_later':'Krellrow→rslots;Krrrow→ellslots;NOclearingdonehere','all396_source_and396_normalized_memberships_exact':True,'maximum_committed_bits':maximum,'no_alpha918_clearing_G3dots_RHSZ2_or_relation_values':True}
 owned=2*count(d)+parsed+count([completed,rhs,current,ev,g,p,c,pins,cursor])+1000;assert owned<cap;dump(ep,ev)
 report={'id':STEM,'all_passed':True,'STOP_only_four_point_dual_rows_completed':True,'completion_not_assumed_in_gate':True,'hard_committed_bits':384,'hard_temporary_bits':1153,'binary_temporary_bound':769,'actual_committed_bits':maximum,'actual_owned_upper':owned,'whole_numeric_cap':cap,'RAW_output_cap':2000000,'no_alpha918_or_G3dot_or_b1_RHSZ2_response':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_committed_bits','actual_owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
