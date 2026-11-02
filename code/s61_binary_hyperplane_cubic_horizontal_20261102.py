#!/usr/bin/env python3
"""Exact d-degree79 and B^12 divisibility of the saved full cubic minor."""
import argparse,functools,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import numpy as np
import sympy as sp
from sympy.polys.rings import ring
from s61_binary_hyperplane_cubic_matrix_20261102 import ROOT,build_matrix,compile_native,digest
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
PRIME=32003
@functools.lru_cache(None)
def mod(c):
 q=sp.Rational(c);assert int(q.q)%PRIME;return int(q.p)%PRIME*pow(int(q.q)%PRIME,-1,PRIME)%PRIME

def read_polys(path,nv):
 lines=Path(path).read_text().splitlines();shape=list(map(int,lines[0].split()));n,m=shape[:2];data=[]
 offset=3 if shape==[42,70,21] else 1
 for z in lines[offset:offset+n*m]:
  z=z.split();nt=int(z[0]);assert len(z)==1+(nv+1)*nt
  data.append([(tuple(map(int,z[2+(nv+1)*k:2+(nv+1)*k+nv])),z[1+(nv+1)*k]) for k in range(nt)])
 return shape,data

def write_matrix(path,shape,data):
 with Path(path).open('w') as f:
  f.write(' '.join(map(str,shape))+'\n')
  for q in data:f.write(str(len(q))+''.join(' '+str(c)+' '+' '.join(map(str,m)) for m,c in q)+'\n')

def leading_input(path,out):
 (n,m),data=read_polys(path,2)
 A=np.array([sum(mod(c)*pow(2,ex[0],PRIME)*pow(3,ex[1],PRIME) for ex,c in q)%PRIME for q in data],dtype=np.int64).reshape(n,m)
 rank,cols=elimination_rank(A,PRIME);_,rows=elimination_rank(A[:,cols].T,PRIME);free=[i for i in range(m) if i not in cols]
 with Path(out).open('w') as f:
  f.write(f'{n} {m} {rank} {len(free)}\n'+' '.join(map(str,rows))+'\n'+' '.join(map(str,cols))+'\n'+' '.join(map(str,free))+'\n')
  for q in data:f.write(str(len(q))+''.join(' '+c+' '+' '.join(map(str,ex)) for ex,c in q)+'\n')
 return dict(rank=rank,rows=rows,columns=cols,free=free,minor_modp=int(determinant_mod(A[np.ix_(rows,cols)],PRIME)))

def transpose(path,out):
 (n,m),data=read_polys(path,2);write_matrix(out,[m,n],[data[i*m+j] for j in range(m) for i in range(n)])

def hungarian(degree):
 n=len(degree);u=[0]*(n+1);v=[0]*(n+1);match=[0]*(n+1);way=[0]*(n+1)
 for i in range(1,n+1):
  match[0]=i;j0=0;minv=[10**9]*(n+1);used=[False]*(n+1)
  while True:
   used[j0]=True;i0=match[j0];delta=10**9;j1=0
   for j in range(1,n+1):
    if not used[j]:
     cur=-degree[i0-1][j-1]-u[i0]-v[j]
     if cur<minv[j]:minv[j]=cur;way[j]=j0
     if minv[j]<delta:delta=minv[j];j1=j
   for j in range(n+1):
    if used[j]:u[match[j]]+=delta;v[j]-=delta
    else:minv[j]-=delta
   j0=j1
   if match[j0]==0:break
  while True:
   j1=way[j0];match[j0]=match[j1];j0=j1
   if j0==0:break
 rw=[-x for x in u[1:]];cw=[-x for x in v[1:]]
 assert all(rw[i]+cw[j]>=degree[i][j] for i in range(n) for j in range(n))
 assert sum(rw)+sum(cw)==sum(degree[match[j]-1][j-1] for j in range(1,n+1))
 return rw,cw

def audit_B_factors(path):
 R,a,e,d=ring('a,e,d',sp.QQ);B=50*a*d+18*d*e**2-9*d*e-2*d+60*e-40
 lines=iter(Path(path).read_text().splitlines());head=next(lines).split();assert head[:2]==['MINOR','0'];net=0;rows=[]
 for i in range(int(head[3])):
  sign=next(lines).split();assert sign[0]=='EXP';mult=int(sign[1]);z=next(lines).split();q=R.from_dict({tuple(map(int,z[2+4*k:5+4*k])):sp.QQ(z[1+4*k]) for k in range(int(z[0]))});order=0
  while not q.rem(B):q=q.exquo(B);order+=1
  net+=mult*order;rows.append(dict(exponent=mult,B_order=order,remaining_terms=len(q)))
 assert net==2
 return dict(net_B_order=net,entries=rows)

def run_horizontal(wd,matrix=None):
 start=time.time();wd=Path(wd);wd.mkdir(parents=True,exist_ok=True)
 source=build_matrix(wd/'source') if matrix is None else None;matrix=wd/'source/schur' if matrix is None else Path(matrix)
 certpath=ROOT/'code/certificates/s61_binary_hyperplane_cubic_first_minor_20261102.json';cert=json.loads(certpath.read_text());assert digest(matrix)==cert['source_schur_sha256']
 (nr,nc,k),allq=read_polys(matrix,3);assert (nr,nc)==(42,70);rr=cert['rows'];cc=cert['columns'];A=[allq[i*nc+j] for i in rr for j in cc];n=36
 write_matrix(wd/'original',[n,n],A);write_matrix(wd/'partial_input',[n,n,n],A)
 with (wd/'partial_input').open('a') as f:f.write('1\n'+' '.join(map(str,list(range(n))*2))+'\n')
 names=['leading_kernel','kernel_peel','partial_determinant','B_jet','B_cokernel_project','leading_primitive'];native={name:ROOT/f'code/s61_binary_hyperplane_{name}_20261102.cpp' for name in names}
 for name,p in native.items():compile_native(p,wd/name)
 env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'};logs={};timings={}
 def run(name,args,label,timeout=180):
  t=time.time();cp=subprocess.run([str(wd/name),*map(str,args)],env=env,capture_output=True,text=True,check=True,timeout=timeout);timings[label]=time.time()-t;logs[label]=cp.stderr;print(label,round(timings[label],3),cp.stderr[-180:],flush=True)
 degree=[[max([m[2] for m,c in A[i*n+j]],default=-100000) for j in range(n)] for i in range(n)];rw,cw=hungarian(degree);assert sum(rw)+sum(cw)==107
 transformed=[[(m[:2]+(rw[i]+cw[j]-m[2],),c) for m,c in A[i*n+j]] for i in range(n) for j in range(n)]
 write_matrix(wd/'d_full0',[n,n],transformed);write_matrix(wd/'d_leading0',[n,n],[[(m[:2],c) for m,c in q if m[2]==0] for q in transformed])
 profiles=[];current=0
 for nullity in [14,9,5]:
  info=leading_input(wd/f'd_leading{current}',wd/f'd_input{current}');assert n-info['rank']==nullity;profiles.append(info)
  run('leading_kernel',[wd/f'd_input{current}',wd/f'd_kernel{current}'],f'd_kernel_{current}')
  nxt=current+nullity;run('kernel_peel',[wd/f'd_full{current}',wd/f'd_input{current}',wd/f'd_kernel{current}',wd/f'd_full{nxt}',wd/f'd_leading{nxt}'],f'd_peel_{current}');current=nxt
 final=leading_input(wd/'d_leading28',wd/'d_input28');assert final['rank']==36 and final['minor_modp']==21773;profiles.append(final)
 run('partial_determinant',[wd/'partial_input',wd/'partial_factors','1',wd/'partial18','18'],'partial18')
 bag=audit_B_factors(wd/'partial_factors')
 # One B-jet is sufficient: only the next constant matrix is used.
 run('B_jet',[wd/'partial18',wd/'B_jet1',wd/'B_leading0','1'],'B_jet1')
 B0=leading_input(wd/'B_leading0',wd/'B_input0');assert B0['rank']==11
 run('leading_kernel',[wd/'B_input0',wd/'B_kernel0'],'B_right_kernel')
 transpose(wd/'B_leading0',wd/'B_left_leading0');left=leading_input(wd/'B_left_leading0',wd/'B_left_input0');assert left['rank']==11
 run('leading_kernel',[wd/'B_left_input0',wd/'B_left_kernel0'],'B_left_kernel')
 run('kernel_peel',[wd/'B_jet1',wd/'B_input0',wd/'B_kernel0',wd/'B_stage7',wd/'B_stage7_leading'],'B_first_peel')
 run('B_cokernel_project',[wd/'B_stage7_leading',wd/'B_left_kernel0',wd/'B_projected'],'B_cokernel_projection')
 run('leading_primitive',[wd/'B_projected',wd/'B_projected_primitive'],'B_projected_primitive')
 small=leading_input(wd/'B_projected_primitive',wd/'B_small_input');assert small['rank']==4
 run('leading_kernel',[wd/'B_small_input',wd/'B_small_kernel'],'B_small_kernel')
 paths=[Path(__file__),certpath,*native.values(),ROOT/'code/s61_binary_deeper_cubic_field_20261101.hpp',ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py',ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py']
 hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
 if source:hashes.update(source['input_sha256'])
 return dict(all_passed=True,scope='One specified maximal minor: exact QQ d-degree79 with p-good leading coefficient; QQ polynomial divisibility by B^12. No all-core rank or S61 exhaustion.',prime=PRIME,source_matrix=source,source_matrix_sha256=digest(matrix),first_minor=cert,hungarian_row_weights=rw,hungarian_column_weights=cw,hungarian_upper=107,d_kernel_coranks=[14,9,5],d_profiles=profiles,exact_d_degree=79,B_partial_factor_audit=bag,B_initial= B0,B_left_initial=left,B_projected=small,B_partial_order_lower=10,B_total_order_lower=12,B_first_jet_only=True,generated_sha256={p.name:digest(p) for p in wd.iterdir() if p.is_file() and p.name not in names},native_logs=logs,stage_seconds=timings,input_sha256=hashes,seconds=time.time()-start)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists()
 with tempfile.TemporaryDirectory(prefix='poisson-hyperplane-horizontal-') as tmp:result=run_horizontal(tmp)
 with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k] for k in ['all_passed','exact_d_degree','B_total_order_lower','seconds']}))
