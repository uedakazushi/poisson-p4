#!/usr/bin/env python3
"""Regenerate all scalar142 normal forms and actions from the original135 equations."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
from fractions import Fraction
from itertools import combinations_with_replacement as cwr
from s61_flag_kernel_residue_raw_engine_20261102 import raw_residual
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_RESIDUE_ACTIONS_20261102.md'

def export_scalar(destination):
    start=time.monotonic();root=ROOT
    initial=json.loads((root/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json').read_text());lifts=json.loads((root/'reports/S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102.json').read_text())['actual_generator_lifts']
    mons=[m for d in range(7)for m in cwr(range(7),d)];ix={m:i for i,m in enumerate(mons)};lowN=792;assert len(mons)==1716 and len(mons[lowN:])==924
    primary=json.loads((ROOT/'reports/S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102.json').read_text())
    assert primary['all_passed']
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in primary['input_sha256'].items())
    with tempfile.TemporaryDirectory(prefix='zero12-top-projection-')as tmp:
     td=Path(tmp);raw,raw_log=raw_residual(initial,td);polys=[]
     for f in initial['initial_relations']:
      terms=[]
      for mon,c in f:
       if any(j in(0,1)for j in mon):continue
       terms.append('('+c+')*'+'*'.join('x_'+str(j-2)for j in mon))
      if terms:polys.append('+'.join(terms))
     out=td/'projection.txt';program='R=QQ[x_0..x_6];J=ideal('+','.join(polys)+');G=gb J;assert(numColumns basis(6,R/J)==15);f=openOut "'+str(out)+'";\n'
     for m in mons[lowN:]:
      e='*'.join('x_'+str(j)for j in m)
      program+='a=listForm(('+e+')%G);f<<toString(#a)<<" ";scan(a,z->{f<<toString(last z)<<" ";scan(first z,v->f<<toString(v)<<" ");});f<<endl;\n'
     program+='close f;';src=td/'projection.m2';src.write_text(program)
     env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(root/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
     r=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(src)],env=env,capture_output=True,text=True,check=True,timeout=30)
     top=[];keys=set()
     for line in out.read_text().splitlines():
      z=line.split();assert len(z)==1+8*int(z[0]);p={}
      for i in range(int(z[0])):
       m=tuple(map(int,z[2+8*i:9+8*i]));p[m]=Fraction(z[1+8*i]);keys.add(m)
      top.append(p)
     keys=sorted(keys);assert len(keys)==15;ki={k:lowN+i for i,k in enumerate(keys)}
     projected=[{i:Fraction(1)}for i in range(lowN)]+[{ki[k]:c for k,c in f.items()}for f in top]
     N=lowN+15
     def project(p):
      ans={}
      for j,c in p.items():
       for k,d in projected[j].items():ans[k]=ans.get(k,Fraction(0))+c*d
      return{k:c for k,c in ans.items()if c}
     def line(p):return str(len(p))+' '+' '.join(f'{i} {c}'for i,c in sorted(p.items()))+'\n'
     ma=td/'maps.txt'
     with ma.open('w')as f:
      f.write(f'7 {N} {lowN} 135\n'+' '.join(map(str,[len(m)for m in mons[:lowN]]+[6]*15))+'\n')
      for j in range(7):
       for k in range(N):f.write(line(projected[ix[tuple(sorted(mons[k]+(j,)))]]if k<lowN else {}))
     projectedraw=[]
     for row in raw:
      rr=project({j:c for j,c in enumerate(row,8)if c});projectedraw.append(rr)
     qq=td/'qq.txt';qq.write_text('\n'.join(' '.join(str(rr.get(j,0))for j in range(8,N))for rr in projectedraw[:135])+'\n')
     mod=td/'mod.txt';prime=32003
     with mod.open('w')as f:
      f.write(f'{prime} 7 {lowN} 135\n'+' '.join(str(len(m))for m in mons[:lowN])+'\n')
      for j in range(7):f.write(' '.join(str(ix[tuple(sorted(m+(j,)))])if len(m)<5 else '-1'for m in mons[:lowN])+'\n')
      for rr in projectedraw[:135]:
       row={j:c.numerator*pow(c.denominator,-1,prime)%prime for j,c in rr.items()if j<lowN};row={j:c for j,c in row.items()if c};f.write(line(row))
      f.write(' '.join('0'for _ in range(lowN))+'\n')
     binaries=[]
     for src in [root/'code/truncated_local_ideal_filtered_trace_20261101.cpp',root/'code/s61_flag_kernel_residue_export_engine_20261102.cpp']:
      binary=td/src.stem;subprocess.run(['c++','-O3','-w','-std=c++17','-I/opt/homebrew/include',str(src),'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(binary)],check=True);binaries.append(binary)
     trace,order,first=[td/n for n in ['trace.txt','order.txt','first.txt']]
     r=subprocess.run([str(binaries[0]),str(mod),str(trace),str(order),str(first)],capture_output=True,text=True,check=True,timeout=30);assert r.stdout.strip()=='665';print('lower trace665, cols807',time.monotonic()-start,flush=True)
     ev=td/'evaluation.txt'
     with ev.open('w')as f:
      f.write('1716\n')
      for i,m in enumerate(mons):f.write((f'{ix[m[:-1]]} {m[-1]} 'if m else '-1 -1 ')+line(projected[i]))
     qfile=td/'queries.txt';qfile.write_text('0\n');answer=td/'answers.txt'
     try:r=subprocess.run([str(binaries[1]),str(ma),str(qq),str(trace),str(ev),str(qfile),str(answer)],capture_output=True,text=True,check=True,timeout=90)
     except subprocess.TimeoutExpired as error:
      raise RuntimeError('Projected exact replay timed out; no report was written') from error
     print(r.stderr,flush=True)

     import shutil
     shutil.copyfile(str(answer)+'.full',destination)
     print('EXPORTED',destination,time.monotonic()-start,flush=True)
     return dict(all_passed=True,seconds=time.monotonic()-start)
