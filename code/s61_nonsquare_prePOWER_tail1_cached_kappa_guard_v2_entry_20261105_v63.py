#!/usr/bin/env python3
"""ONE fixed tail1 value, usable only after source/code acceptance and parent GO."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63'
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_GUARD_CONTROL_V2_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def inspect(p):
 with p.open()as f:
  assert f.readline().strip()=='J1KP001';cursor=f.readline().strip();reason=f.readline().strip();n=int(f.readline());total=0
  for i in range(n):
   node,use=map(int,f.readline().split());assert node==i;m=int(f.readline());total+=m;old=None
   for _ in range(m):
    k,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=359
    assert old is None or old<k;old=k
  assert not f.read().strip()
 return dict(version='J1KP001',cursor=cursor,reason=reason,all_maps_parsed=True,coefficient_key_pairs=total,
  exact_pending_operation='phase0 literal atom cursor; phase1 pair=-2 means node complete and uses released, so advance to next sorted source node; phase1 pair=-1 is current uncommitted operation; phases2/3/4 fixed dot pass and pair indices before commit.',
  immutable_instruction_file=str(BASE/'native-input.txt'),no_resume_authorized=True)
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 report=json.loads((BASE/'report.json').read_text());assert report['all_passed']and report['value_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for name,h in report['evidence_sha256'].items():
  p=BASE/name;assert sha(p)==h;pins[str(p.relative_to(ROOT))]=h
 me=Path(__file__);pins[str(me.relative_to(ROOT))]=sha(me)
 base=sum((ROOT/p).stat().st_size for p in pins)+len(json.dumps(pins))+10000
 assert base+4700000+100000<9517864
 (out/'phase.json').write_text(json.dumps(dict(phase='prepared_inputs_verified_before_value',input_sha256=pins,
  expression='T1=B*kappa*M15[rows2,4,5,8;original32column1]',single_target=True,base_RAW=base),indent=2)+'\n')
 cmd=[str(BASE/'native-worker'),'--run',str(BASE/'native-input.txt'),str(out),str(base)]
 with (out/'native.log').open('w')as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=58)
 result=json.loads((out/'native-report.json').read_text())
 if r.returncode!=0:
  assert result.get('prefix_only');result['prefix_inspection']=inspect(out/'J1KP001-prefix.txt')
 else:
  assert result['all_passed']and result['reverse_all_four_products_zero']and result['three_passes_complete']
  counts=[]
  with (out/'integer-column4-and-result.txt').open()as f:
   assert f.readline().strip()=='J1KV001'
   for expect in [66,94,122,150,155]:
    node,D=map(int,f.readline().split());assert node==expect and D>0;m=int(f.readline());old=None;first=None
    for _ in range(m):
     k,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=359
     assert old is None or old<k;old=k
     if first is None:first=[k,str(c)]
    counts.append(dict(node=node,D=str(D),terms=m,first=first))
   assert not f.read().strip()
  assert counts[-1]['terms']==result['result_terms'];result['all5_outputs_readback']=counts
 total=base+sum(p.stat().st_size for p in out.iterdir())+10000;assert total<9517864
 final=dict(id=ID,all_passed=result['all_passed'],frozen=True,only_one_tail1_value=True,
  no_full15_rank_or_actual_bad_ideal_equality=True,result=result,input_sha256=pins,
  total_RAW_with_inputs_upper=total,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(final,indent=2)+'\n');print(json.dumps({k:v for k,v in final.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--inspect-prefix',type=Path);a=p.parse_args()
 if a.inspect_prefix:print(json.dumps(inspect(a.inspect_prefix),indent=2))
 elif a.output:main(a.output.resolve())
 else:p.error('--output or --inspect-prefix required')
