#!/usr/bin/env python3
"""ONE known C4 quotient, to run ONLY after fresh source/code audit and parent GO."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63'
ID='S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def inspect(p):
 def readpoly(f):
  n=int(f.readline());old=None
  for _ in range(n):
   k,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518 and(old is None or old<k);old=k
  return n
 with p.open()as f:
  assert f.readline().strip()=='T1C4QP001';cursor=f.readline().strip();why=f.readline().strip();D,L=map(int,f.readline().split());n=int(f.readline())
  for _ in range(n):k,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518
  nr=int(f.readline());root=0
  for _ in range(nr):
   m=int(f.readline());root+=m
   for __ in range(m):assert abs(int(f.readline())).bit_length()<=518
  nc=int(f.readline());const=sum(readpoly(f)for _ in range(nc));q,s,m=[readpoly(f)for _ in range(3)]
  assert not f.read().strip()
 return dict(version='T1C4QP001',cursor=cursor,reason=why,input=n,root_slots=root,constant_terms=const,Q=q,slab=s,remainder=m,
  pending='phase2 formulai/inputj/rootr; phase3 whole slab complete and next transferkey uncommitted; phase4 constantC power cursor; phase5 next input scale; phase6 pendingQkey/Ckey; phase7 all multiplyback complete beforeoutput',
  fixed_metadata='100powers0..20 andchoose(h+3,3)0..16 reconstructed from immutable literal operator, not unknown cached values',no_auto_resume=True)
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);r=json.loads((BASE/'report.json').read_text());assert r['all_passed']and r['value_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for name,h in r['evidence_sha256'].items():
  p=BASE/name;assert sha(p)==h;pins[str(p.relative_to(ROOT))]=h
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 base=sum((ROOT/p).stat().st_size for p in pins)+20000;assert base+8410960+200000<=10610960
 (out/'phase.json').write_text(json.dumps(dict(phase='BEFOREquotient',only_T1_C4=True,input_sha256=pins),indent=2)+'\n')
 with (out/'native.log').open('w')as f:run=subprocess.run([str(BASE/'native-worker'),'--run',str(BASE/'native-input.txt'),str(out),str(base)],stdout=f,stderr=subprocess.STDOUT,timeout=58)
 n=json.loads((out/'native-report.json').read_text())
 if run.returncode:
  assert n.get('prefix_only');n['prefix_inspection']=inspect(out/'T1C4QP001-prefix.txt')
 else:
  assert n['all_passed']and n['full_C4_multiplyback_zero']and n['all_completed_slabs_checked']
  with (BASE/'native-input.txt').open()as f:f.readline();D1=int(f.readline())
  with (out/'quotient-integer-common-D.txt').open()as f:
   assert f.readline().strip()=='T1C4QV001';D=int(f.readline());assert D==D1*100**20;m=int(f.readline());assert m==n['quotient_terms']
   old=None
   for _ in range(m):
    k,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518 and(old is None or old<k);old=k
    a,x=divmod(k,1048576);e,d=divmod(x,1024);assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50
   assert not f.read().strip()
  n['all_output_coefficients_readback']=True
 total=base+sum(p.stat().st_size for p in out.iterdir())+20000;assert total<=10610960
 result=dict(id=ID,all_passed=n['all_passed'],frozen=True,only_T1_known_C4_quotient=True,result=n,
  no_GCD_other_factor_tail_or_bad_locus_claim=True,RAW_with_inputs_upper=total,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--inspect-prefix',type=Path);a=p.parse_args()
 if a.inspect_prefix:print(json.dumps(inspect(a.inspect_prefix),indent=2))
 elif a.output:main(a.output.resolve())
 else:p.error('--output or --inspect-prefix required')
