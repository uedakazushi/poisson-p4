#!/usr/bin/env python3
"""ONE fixed primitive-B4 quotient. Run ONLY after fresh ROOT audit/ONE GO."""
from pathlib import Path
import argparse,hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_WORKER_CUSTODY_20261105_v63'
ID='S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def inspect(p):
 def poly(f):
  n=int(f.readline());prev=None
  for _ in range(n):
   key,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518 and(prev is None or prev<key);prev=key
  return n
 with p.open()as f:
  assert f.readline().strip()=='Q1B4QP001';cursor=f.readline().strip();assert len(cursor.split())==16
  why=f.readline().strip();terms=list(map(int,f.readline().split()));digests=list(map(int,f.readline().split()));assert len(terms)==len(digests)==4
  D,L,Lr,active=map(int,f.readline().split());sizes=[poly(f)for _ in range(4)];assert not f.read().strip()
 return dict(version='Q1B4QP001',cursor=cursor,reason=why,stage_terms=terms,stage_telemetry=digests,
  R=sizes[0],Q=sizes[1],fixed_power=sizes[2],next_power=sizes[3],owned_current=2*sum(sizes)+16000,
  cursor_fields='phase,pass,Bterm,factorh,sourcei,active_ready,activeQkey,factorkey,verifyQkey,verifyFkey,products,scalars,updates,peak,inputFNV,rereadFNV',
  active_ready='0 noactiveQ;1 quotientscalarchecked/notinserted;2 Qinserted and Bterm is nextpending update',
  restore='R/Q/F/nextF,D/L/Lr/activequotient andallcursors saved; immutable4811input reusedonlyforpendingread; fixed6Bterms fromliteralcode',
  FNV='metadata telemetry only, not a proof hash; inputSHA256 rechecked before/after actualjob',no_auto_resume=True)
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);r=json.loads((BASE/'report.json').read_text());assert r['all_passed']and r['value_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for n,h in r['input_sha256'].items():
  p=ROOT/n;assert sha(p)==h;pins[n]=h
 for n,h in r['evidence_sha256'].items():
  p=BASE/n;assert sha(p)==h;pins[str(p.relative_to(ROOT))]=h
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+20000;assert raw+6000000+200000<=8000000
 (out/'phase.json').write_text(json.dumps(dict(phase='BEFORE_B4quotient',input_sha256=pins,only_known_B4=True),indent=2)+'\n')
 with(out/'native.log').open('w')as f:
  run=subprocess.run([str(BASE/'native-worker'),'--run',str(BASE/'native-input.txt'),str(out),str(raw)],stdout=f,stderr=subprocess.STDOUT,timeout=58)
 # The native reread is bound to the same immutable SHA input, not an oracle.
 for n,h in pins.items():assert sha(ROOT/n)==h,n
 n=json.loads((out/'native-report.json').read_text())
 if run.returncode:
  assert n.get('prefix_only');n['prefix_inspection']=inspect(out/'Q1B4QP001-prefix.txt');assert n['prefix_inspection']['owned_current']<=52826
 else:
  assert n['all_passed']and n['all_four_stage_remainders_zero']and n['full_B4_multiplyback_zero']and n['original_input_reread_matched']
  assert len(n['stage_terms'])==4 and n['stage_terms'][-1]==n['quotient_terms']
  with(BASE/'native-input.txt').open()as f:f.readline();D0=int(f.readline())
  with(out/'quotient-integer-common-D.txt').open()as f:
   assert f.readline().strip()=='Q1B4QV001';D=int(f.readline());assert D==D0;count=int(f.readline());assert count==n['quotient_terms'];old=None
   for _ in range(count):
    key,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518 and(old is None or old<key);old=key
    a,x=divmod(key,1048576);e,d=divmod(x,1024);assert 0<=a<=12 and 0<=e<=25 and 0<=d<=14 and a+e+d<=38
   assert not f.read().strip()
  n['all_final_quotient_coefficients_readback']=True
 total=raw+sum(p.stat().st_size for p in out.iterdir())+20000;assert total<=8000000
 result=dict(id=ID,all_passed=n['all_passed'],frozen=True,only_known_primitive_B4_quotient=True,result=n,
  Q2_definition='Q2=Q1/B4=T1/(a3*C4*B4), numeratorinteger/originalD197; notoriginalS61quadricq2',
  no_GCD_other_factor_tail_or_bad_locus_claim=True,RAW_with_inputs_upper=total,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--inspect-prefix',type=Path);a=p.parse_args()
 if a.inspect_prefix:print(json.dumps(inspect(a.inspect_prefix),indent=2))
 elif a.output:main(a.output.resolve())
 else:p.error('--output or --inspect-prefix required')
