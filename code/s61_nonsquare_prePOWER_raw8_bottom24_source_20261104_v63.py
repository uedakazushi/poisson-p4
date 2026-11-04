#!/usr/bin/env python3
"""Expose only the already-computed24 source entries, from saved native input."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports/S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_ENTRY_20261104_v63';old=read(base/'report.json');assert old['all_passed']and old['normalized_minor_is_zero']
 opath=ROOT/'reports/S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_FUSED_PREFLIGHT_20261104_v63/operator-and-cost.json';op=read(opath);N=op['scalar_nodes'];roots=[i for row in op['matrix8x8_after_row_a_normalization'][4:]for i in row[:6]];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);todo+=N[i]['args']
 assert len(need)==190;order=sorted(need);index={x:i for i,x in enumerate(order)};ss=lambda i:N[i]['bound']['support']if N[i]['bound']else 0
 bound=sum(map(ss,need))+5*max(map(ss,need))+15000;assert bound==99960 and bound<150000
 tokens=iter((base/'native-input.txt').read_text().split());assert next(tokens)=='8048'and next(tokens)=='2094';records={}
 for i in range(2094):
  typ=next(tokens)
  if typ=='x':continue
  upper,den=next(tokens),next(tokens)
  if typ=='p':
   sh,unit=next(tokens),next(tokens);n=int(next(tokens));p=[[next(tokens)for _ in range(4)]for j in range(n)];record=(typ,sh,unit,p)
  elif typ=='k':record=(typ,int(next(tokens)),next(tokens))
  else:record=(typ,int(next(tokens)),int(next(tokens)))
  if i in need:records[i]=record
 assert len(records)==190
 with(out/'input.txt').open('w')as f:
  f.write('190\n')
  for i in order:
   r=records[i];f.write(f'{i} {r[0]} ')
   if r[0]=='p':
    f.write(f'{r[1]} {r[2]}\n{len(r[3])}\n')
    for row in r[3]:f.write(' '.join(row)+'\n')
   elif r[0]=='k':f.write(f'{index[r[1]]} {r[2]}\n')
   else:f.write(f'{index[r[1]]} {index[r[2]]}\n')
  f.write(' '.join(str(index[i])for i in roots)+'\n')
 dump(out/'source-binding.json',dict(saved_native_input_sha256=sha(base/'native-input.txt'),source_roots24=roots,minimal_closure190=order,source_current15_rows=[2,4,5,8],original32_columns=[22,11,13,25,15,19],all_row_a_shifts_zero=True,kept_all190_plus_scratch_guard=bound,no_minor_replay=True))
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=55);(out/'compile.log').write_text(run.stdout);assert run.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(out/'input.txt')+' '+str(out/'native.json');run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,57-(time.monotonic()-start)));(out/'native.log').write_text(run.stdout);assert run.returncode==0
 d=read(out/'native.json');assert d['all_passed'];paths=[base/'report.json',base/'native-input.txt',opath,Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md'];size=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p!=exe);assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_entry_values_only=True,source_nodes=190,source_atoms=46,explicit_output_terms=sum(len(p)for row in d['bottom4x6']for p in row),forecast_terms=bound,actual_peak_terms=d['peak_terms_with15000reserve'],no_other_column_minor_rank_or_point=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-start);dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
