#!/usr/bin/env python3
"""ONE canonical Q1/B first-nonzero job; execute ONLY after fresh ROOT audit/ONE GO."""
from pathlib import Path
import argparse,hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_WORKER_CUSTODY_20261105_v63'
ID='S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def inspect(p):
 def poly(f):
  n=int(f.readline());prev=None
  for _ in range(n):
   key,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518 and(prev is None or prev<key);prev=key
  return n
 with p.open()as f:
  assert f.readline().strip()=='Q1BHP001';cursor=f.readline().strip();assert len(cursor.split())==15
  why=f.readline().strip();D,L,activebase=map(int,f.readline().split());n=int(f.readline())
  for _ in range(n):key,c=map(int,f.readline().split());assert c and abs(c).bit_length()<=518
  groups=int(f.readline());nr=0
  for _ in range(groups):
   m=int(f.readline());nr+=m
   for __ in range(m):assert abs(int(f.readline())).bit_length()<=518
  np=poly(f);nn=poly(f);groupsS=int(f.readline());ns=sum(poly(f)for _ in range(groupsS));na=poly(f)
  assert not f.read().strip()
 return dict(version='Q1BHP001',cursor=cursor,reason=why,input=n,Rslots=nr,Pcurrent=np,Pnext=nn,Scache=ns,accumulator=na,
  owned_current=2*(n+nr+np+nn+ns+na)+16000,
  cursor_fields='phase,t,h,Rindex,Rterm,b1term,inputindex,k,Pkey,Skey,base_ready,products,scalings,updates,peak',
  restoration='allinput/R/Pcurrent/Pnext/S/A/D/L andactivebase stored; choose19x19 regeneratedfromliteral constants only; mathematical pendingproduct has notcommitted',
  denominator='D197*b1^18 SYMBOLIC; b1 NOT assumed FULLUunit',no_auto_resume=True)
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);r=json.loads((BASE/'report.json').read_text());assert r['all_passed']and r['value_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for n,h in r['input_sha256'].items():
  p=ROOT/n;assert sha(p)==h;pins[n]=h
 for n,h in r['evidence_sha256'].items():
  p=BASE/n;assert sha(p)==h;pins[str(p.relative_to(ROOT))]=h
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+20000;assert raw+5000000+200000<=8000000
 (out/'phase.json').write_text(json.dumps(dict(phase='BEFORE_Bvaluation',source_pins=pins,canonical_Q1_ONLY=True),indent=2)+'\n')
 with(out/'native.log').open('w')as f:
  run=subprocess.run([str(BASE/'native-worker'),'--run',str(BASE/'native-input.txt'),str(out),str(raw)],stdout=f,stderr=subprocess.STDOUT,timeout=58)
 n=json.loads((out/'native-report.json').read_text())
 if run.returncode:
  assert n.get('prefix_only');n['prefix_inspection']=inspect(out/'Q1BHP001-prefix.txt');assert n['prefix_inspection']['owned_current']<=37436
 else:
  assert n['all_passed']and n['reverse_first_layer_zero']and n['reproduced_nonzero']and n['b1_NOT_assumed_FULLU_unit']
  k=n['B_valuation'];assert 0<=k<=18
  assert(out/'completed-zero-orders.txt').read_text().splitlines()==[f'{j} ZERO'for j in range(k)]
  with(BASE/'native-input.txt').open()as f:f.readline();D0=int(f.readline())
  with(out/'first-nonzero-B-Hasse.txt').open()as f:
   assert f.readline().strip()=='Q1BHV001';kk,D,be=map(int,f.readline().split());assert kk==k and D==D0 and be==18
   m=int(f.readline());assert m==n['first_nonzero_terms'];old=None
   for _ in range(m):
    a,e,c=map(int,f.readline().split());assert 0<=a<=34 and 0<=e<=69 and a+e<=86-k and c and abs(c).bit_length()<=518
    assert old is None or old<(a,e);old=(a,e)
   assert not f.read().strip()
  n['full_first_layer_and_symbolic_denominator_readback']=True
 total=raw+sum(p.stat().st_size for p in out.iterdir())+20000;assert total<=8000000
 result=dict(id=ID,all_passed=n['all_passed'],frozen=True,only_canonical_Q1_known_B_valuation=True,result=n,
  denominator='D197*b1^18 symbolic; ordinaryBvaluation byprimitiveGauss, no pointwiseb1chart',
  no_Bquotient_GCD_other_tail_or_bad_locus_claim=True,RAW_with_inputs_upper=total,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--inspect-prefix',type=Path);a=p.parse_args()
 if a.inspect_prefix:print(json.dumps(inspect(a.inspect_prefix),indent=2))
 elif a.output:main(a.output.resolve())
 else:p.error('--output or --inspect-prefix required')
