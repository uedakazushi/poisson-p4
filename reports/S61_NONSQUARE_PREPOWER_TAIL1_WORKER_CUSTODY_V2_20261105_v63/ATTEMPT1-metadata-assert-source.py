#!/usr/bin/env python3
"""Reuse completed literal cut; repair ONLY native runtime linkage and finish gate."""
import argparse,hashlib,json,math,resource,shutil,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TAIL1_WORKER_CUSTODY_V2_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if type(x)is int:return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except ValueError:return 0
 if isinstance(x,list):return sum(map(numeric,x))
 if isinstance(x,dict):return sum(map(numeric,x.values()))
 return 0
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 old=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_WORKER_CUSTODY_20261105_v63'
 failure=json.loads(pin(old/'failed-preparation-receipt.json').read_text())
 names=['selected-plan.json','literal-atoms58.jsonl','native-input.txt','literal-normalization-witness.json']
 for n in names:
  p=pin(old/n);assert sha(p)==failure['byte_hashes'][n];shutil.copyfile(p,out/n);assert sha(p)==sha(out/n)
 plan=json.loads((out/'selected-plan.json').read_text());NN=plan['nodes'];cols=plan['column_roots'];ks=plan['cached_kappa_nodes']
 needed=set();todo=cols[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);todo+=NN[i]['args']
 assert len(needed)==149 and cols==[66,94,122,150]and ks==[151,152,153,154]
 total=maxatom=bits=0;ownplan=numeric(plan)
 with (out/'literal-atoms58.jsonl').open()as source,(out/'native-input.txt').open()as native:
  assert native.readline().strip()=='J1K001 155 149'
  assert list(map(int,native.readline().split()))==cols
  assert list(map(int,native.readline().split()))==ks
  Dfinal=int(native.readline());assert Dfinal==int(NN[161]['QQ_D'])
  assert list(map(int,native.readline().split()))==sorted(needed)
  for i,n in enumerate(NN[:155]):
   fields=native.readline().split();assert int(fields[0])==i and fields[-1]==n['QQ_D']
  assert int(native.readline())==58
  for line in source:
   a=json.loads(line);i=a['node'];n=NN[i];p,den=a['original_num_den'];D=int(n['QQ_D']);u=Q(a['QQ_divisor'])
   assert a['a_shift']==n['a_shift']and a['QQ_divisor']==n['QQ_divisor']and a['reference']==n['reference']
   node,count=map(int,native.readline().split());assert node==i and count==len(p)
   total+=len(p);maxatom=max(maxatom,numeric(a));L=0
   for m,c in p:
    k,v=map(int,native.readline().split());cc=Q(c)
    num=cc.numerator*u.denominator*D;dd=cc.denominator*u.numerator
    assert v*dd==num and k==(m[0]+n['a_shift'])*1048576+m[1]*1024+m[2]
    L+=abs(v)
    h=max(abs(cc.numerator).bit_length()+u.denominator.bit_length()+D.bit_length(),
      cc.denominator.bit_length()+abs(u.numerator).bit_length(),abs(v).bit_length()+cc.denominator.bit_length()+abs(u.numerator).bit_length())+2
    bits=max(bits,h);assert h<2048
   assert L==int(n['integer_l1_upper'])
  assert not native.read().strip()
 assert total==14188
 for i in sorted(needed):
  n=NN[i]
  if not n['args']:continue
  x,y=(NN[j]for j in n['args']);D=int(n['QQ_D']);di,dj=int(x['QQ_D']),int(y['QQ_D']);L=int(n['integer_l1_upper']);li,lj=int(x['integer_l1_upper']),int(y['integer_l1_upper'])
  if n['op']=='mul_bound_only':assert D==di*dj and L==li*lj
  else:assert D%di==D%dj==0 and L==li*(D//di)+lj*(D//dj)
 maxD=max([Dfinal.bit_length()]+[int(n['QQ_D']).bit_length()for n in NN[:155]])
 srcpairs=sum(NN[n['args'][0]]['bound']['support']*NN[n['args'][1]]['bound']['support']for n in NN[:151]if n['op']=='mul_bound_only')
 prep=ownplan+maxatom+17000;assert prep<91676
 cpp=pin(ROOT/'code/s61_nonsquare_tail1_cached_kappa_20261105_v63.cpp')
 driver=pin(ROOT/'code/s61_nonsquare_prePOWER_tail1_cached_kappa_entry_20261105_v63.py')
 me=pin(Path(__file__))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','tail1-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(driver),'--help'])]:
  h=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(h.stdout+h.stderr)
 own=dict(limit91676=91676,coefficient_and_packed_key_cells_per_term=2,ALL_native_metadata_and_scalar_scratch_reserve=15100,
  native_live_terms_limit=38288,source_column_nodes149=149,input_literal_terms14188=14188,
  column_terms9268=9268,cached_kappa_terms4920=4920,source_preparation_owned_upper=prep,
  literal_integer_normalization_temp_bits_upper=bits,native_coefficient_bits_upper=359,
  native_denominator_temp_bits_upper=2*maxD+1,preparation_hard_integer_bound2048=2048,
  source_mul_pairs_upper=srcpairs,dot_pairs_one_pass_upper=14283624,dot_pairs_three_pass_upper=42850872,
  verification='forward; reverse subtraction to literal ZERO; forward reproduction, one accumulator, no product map',
  prefix_version='J1KP001',prefix_parser='entry driver --inspect-prefix FILE; streaming, no auto-resume',
  guard_cursor='source149-node operation or fixed pass/pair and inner product indices; before-commit stop',
  soft_seconds45=45,hard_seconds60=60,RAM_bytes=2000000000,
  NO_column_or_contraction_value_run=True)
 dump('ownership-and-worker-proof.json',own)
 raw=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in [cpp,driver,me])
 future=raw+4700000+100000;assert future<9517864
 dump('RAW-closure.json',dict(actual_prepared_bytes_plus_three_codes=raw,prefix_OR_success_upper=4700000,
  all_future_reports_logs_receipts_reserve=100000,total_upper=future,hard_limit=9517864,
  no_atomic_doublebuffer_no_second_accumulator=True,no_old_checkpoint_Q2_or572_read_by_value_worker=True))
 dump('failure-repair.json',dict(old_failure_receipt=str((old/'failed-preparation-receipt.json').relative_to(ROOT)),
  old_failure_sha256=sha(old/'failed-preparation-receipt.json'),completed_source4_files_reused_byte_for_byte=names,
  repair_only='native executable rpath to the already installed libflint.so.19 directory; source arithmetic code unchanged',
  no_native_coefficient_arithmetic_before_or_after_repair=True,no_failed_bytes_changed=True))
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,strict_owned_limit91676=91676,
  literal_column_terms=9268,cached_kappa_terms=4920,column_DAG_nodes=149,
  worker_and_driver_help_passed=True,worker_compiled=True,value_worker_NOT_run=True,
  native_coefficient_bits359=359,preparation_temp_bits=bits,RAW_total_upper=future,RAW_hard_limit=9517864,
  source_integer_scalar_pair_upper=srcpairs+42850872,independent_source_code_review_pending=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
