#!/usr/bin/env python3
"""Actual ONE T0 C-Hasse job; NOT to execute before independent audit and GO."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63'
ID='S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def inspect(p):
 with p.open()as f:
  assert f.readline().strip()=='T0CHP001';cursor=f.readline().strip();why=f.readline().strip();D,input_l1=map(int,f.readline().split());n=int(f.readline())
  for _ in range(n):key,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518
  groups=int(f.readline());power_cells=0
  for _ in range(groups):
   m=int(f.readline());power_cells+=m
   for __ in range(m):assert abs(int(f.readline())).bit_length()<=518
  m=int(f.readline())
  for _ in range(m):key,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518
  assert not f.read().strip()
 return dict(version='T0CHP001',cursor=cursor,reason=why,input_records=n,power_slots=power_cells,accumulator_terms=m,
  immutable_instructions=str(BASE/'native-input.txt'),no_auto_resume=True)
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);r=json.loads((BASE/'report.json').read_text());assert r['all_passed']and r['value_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for name,h in r['evidence_sha256'].items():
  p=BASE/name;assert sha(p)==h;pins[str(p.relative_to(ROOT))]=h
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+20000;assert raw+4600000<8000000
 (out/'phase.json').write_text(json.dumps(dict(phase='before_C_valuation',source_pins=pins,only_T0=True),indent=2)+'\n')
 with (out/'native.log').open('w')as f:run=subprocess.run([str(BASE/'native-worker'),'--run',str(BASE/'native-input.txt'),str(out),str(raw)],stdout=f,stderr=subprocess.STDOUT,timeout=58)
 n=json.loads((out/'native-report.json').read_text())
 if run.returncode:
  assert n.get('prefix_only');n['prefix_inspection']=inspect(out/'T0CHP001-prefix.txt')
 else:
  assert n['all_passed']and n['reverse_first_layer_zero']and n['reproduced_nonzero']
  k=n['C_valuation'];assert 0<=k<=20
  assert (out/'completed-zero-orders.txt').read_text().splitlines()==[f'{i} ZERO'for i in range(k)]
  with (BASE/'native-input.txt').open()as f:f.readline();D0=int(f.readline())
  with (out/'first-nonzero-C-Hasse.txt').open()as f:
   assert f.readline().strip()=='T0CHV001';kk,D=map(int,f.readline().split());assert kk==k and D==D0*100**20
   m=int(f.readline());assert m==n['first_nonzero_terms'];prev=None
   for _ in range(m):
    e,d,c=map(int,f.readline().split());assert 0<=e<=80 and 0<=d<=18 and c and abs(c).bit_length()<=518
    assert prev is None or prev<(e,d);prev=(e,d)
   assert not f.read().strip()
  n['full_first_layer_readback']=True
 total=raw+sum(p.stat().st_size for p in out.iterdir())+20000;assert total<8000000
 result=dict(id=ID,all_passed=n['all_passed'],frozen=True,only_T0_known_C_valuation=True,result=n,
  no_quotient_GCD_other_tail_or_bad_locus_claim=True,RAW_with_inputs_upper=total,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--inspect-prefix',type=Path);a=p.parse_args()
 if a.inspect_prefix:print(json.dumps(inspect(a.inspect_prefix),indent=2))
 elif a.output:main(a.output.resolve())
 else:p.error('--output or --inspect-prefix required')
