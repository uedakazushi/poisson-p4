#!/usr/bin/env python3
"""Prepared1120 knownD coefficient expressions; no source callback values."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_INPUT_GATE_20261105_v63/'
GATE_SHA='552d60c970aefd792af4ab4ae2fe622441712e1618be43fbb29d908b8ec9e253'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {int(i):F(x) for i,x in v if F(x)}
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};cursor={};pending={};acc=F(0);D=[];maximum=0
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<1_000_000
 g=read(GATE+'report.json');assert g['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'pruned-D-coefficient-functional-and-ten-pair-source-input.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody);assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(Path(__file__))
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([g,p,c,pins])<1000 and g['source_positive_prior_bits']==92 and g['temporary_prior_bits']==277 and g['future_value_numeric_bound']==14334 and count(d)==2851
 C=[[[F(x) for x in row] for row in ch] for ch in d['Cfull_two_by2_by56']];rho=list(map(parse,d['rho_two_by1144_sparse_rows']));M=list(map(parse,d['Mnew0_columns56']));pairs=d['quadratic_pair_order44'];pindex={tuple(v):i for i,v in enumerate(pairs)};ten=list(map(tuple,d['active_pairs10']));assert len(M)==56 and len(ten)==10
 parsed=count([C,rho,M,pindex,ten]);cap=14334;before=2*count(d)+parsed+2240+1000;assert before<cap
 dump('phase.json',{'phase':'source_peer_custody_and_prior_PASS_before_knownD_coefficients','gate_sha256':GATE_SHA,'source_positive_prior92':92,'temporary_prior277':277,'hard_guards':[384,1153],'numeric_cap':cap,'prearith_live_upper':before,'no_RHS_sourceZ2_or_upper_values':True})
 def stop(reason,x):
  state={'reason':reason,'cursor':dict(cursor),'pending_operation':dict(pending),'accumulator_before_operation':str(acc),'uncommitted_fraction':str(x),'completed_D_prefix':D,'input_sha256':pins,'no_retry_or_cap_raise':True};assert 2*count(d)+parsed+3*count(state)+1000<cap;dump('guard-frontier.json',state);raise RuntimeError(reason)
 def ck(x):
  nonlocal maximum
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());maximum=max(maximum,b)
  if b>92:stop('source92_before_commit',x)
  return x
 def mul(a,b):pending.update(operation='multiply',operand_a=str(a),operand_b=str(b));return ck(a*b)
 def sub(a,b):pending.update(operation='subtract',operand_a=str(a),operand_b=str(b));return ck(a-b)
 for t in range(2):
  ch=[];D.append(ch)
  for k,pair in enumerate(ten):
   row=[];ch.append(row);mi=pindex[pair]
   for j,col in enumerate(M):
    q=C[t][1][j] if pair==(7,7) else C[t][0][j] if pair==(7,8) else F(0)
    for target,x in col.items():
     cursor={'phase':'knownD_coefficient_contraction','channel':t,'pair10':k,'source_row56':j,'F1target26':target};acc=q;term=mul(rho[t].get(44*target+mi,F(0)),x);q=sub(q,term)
    DD,U=map(int,d['positive_D_two_by10_by56'][t][k][j]);assert DD%q.denominator==0 and abs(q.numerator)*DD<=U*q.denominator;row.append(str(q))
 ep='actual-D-coefficients-two-by10-by56.json';ev={'D_two_by10_by56':D,'active_pairs10':d['active_pairs10'],'all1120_source_envelopes_checked':True,'formula':d['formula'],'full_upper_source_formula':d['full_upper_source_formula'],'no_source_callback_value_or_honest_second_gradient_evaluated':True,'OLDfour_pivot_g2_frontier_and_native_L2S0_retained':True,'scope':'OnlycoefficientfunctionalD.NORHSZ2/firstcross/secondor3rdfunctions/b1/upper617/HF/length.'}
 owned=2*count(d)+parsed+count([D,ev,g,p,c,pins,cursor,pending])+1000;assert owned<cap;dump(ep,ev)
 report={'id':STEM,'all_passed':True,'coefficient_functional_D_only':True,'source_positive_prior92':92,'temporary_prior277':277,'hard_guards':[384,1153],'actual_maximum_bits':maximum,'actual_owned_upper':owned,'whole_numeric_cap':cap,'RAW_output_cap':1000000,'no_source_RHS_or_upper_response_values':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_maximum_bits','actual_owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
