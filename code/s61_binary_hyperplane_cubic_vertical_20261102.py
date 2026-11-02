#!/usr/bin/env python3
"""Exact vertical degree, leading-coefficient and known-factor certificates."""
import argparse,hashlib,json,os,subprocess,sys,tempfile,time
from fractions import Fraction
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from s61_binary_hyperplane_cubic_matrix_20261102 import build_matrix,compile_native,digest
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
P=32003

def assignment(cost):
 n=len(cost);u=[0]*(n+1);v=[0]*(n+1);match=[0]*(n+1);way=[0]*(n+1)
 for i in range(1,n+1):
  match[0]=i;j0=0;mv=[10**20]*(n+1);used=[False]*(n+1)
  while True:
   used[j0]=True;i0=match[j0];delta=10**20;j1=0
   for j in range(1,n+1):
    if not used[j]:
     cur=cost[i0-1][j-1]-u[i0]-v[j]
     if cur<mv[j]:mv[j]=cur;way[j]=j0
     if mv[j]<delta:delta=mv[j];j1=j
   for j in range(n+1):
    if used[j]:u[match[j]]+=delta;v[j]-=delta
    else:mv[j]-=delta
   j0=j1
   if match[j0]==0:break
  while True:
   j1=way[j0];match[j0]=match[j1];j0=j1
   if j0==0:break
 ans=[0]*n
 for j in range(1,n+1):ans[match[j]-1]=j-1
 assert all(cost[i][j]>=u[i+1]+v[j+1] for i in range(n) for j in range(n))
 assert sum(cost[i][ans[i]] for i in range(n))==sum(u[1:])+sum(v[1:])
 return -sum(u[1:])-sum(v[1:]),[-x for x in u[1:]],[-x for x in v[1:]]

def mod(q):
 q=Fraction(q);assert q.denominator%P;return q.numerator%P*pow(q.denominator%P,-1,P)%P

def materialize_minor_input(schur,destination):
 candidate=json.loads((ROOT/'code/certificates/s61_binary_hyperplane_cubic_first_minor_20261102.json').read_text())
 assert digest(schur)==candidate['source_schur_sha256'];lines=Path(schur).read_text().splitlines();assert lines[0]=='42 70 21';assert len(lines)==3+42*70
 for line in lines[3:]:
  z=line.split();assert len(z)==1+4*int(z[0])
  for k in range(int(z[0])):mod(z[1+4*k])
 Path(destination).write_text('42 70 36\n'+'\n'.join(lines[3:])+'\n1\n'+' '.join(map(str,candidate['rows']+candidate['columns']))+'\n')
 return candidate

def profile(source,destination):
 lines=Path(source).read_text().splitlines();n,m=map(int,lines[0].split());assert len(lines)==1+n*m;M=[]
 for line in lines[1:]:
  z=line.split();nt=int(z[0]);assert len(z)==1+3*nt
  M.append(sum(mod(z[1+3*i])*pow(2,int(z[2+3*i]),P)*pow(3,int(z[3+3*i]),P) for i in range(nt))%P)
 M=np.array(M,dtype=np.int64).reshape(n,m);rank,cc=elimination_rank(M,P);_,rr=elimination_rank(M[:,cc].T,P);free=[j for j in range(m) if j not in cc];det=determinant_mod(M[np.ix_(rr,cc)],P)
 Path(destination).write_text(f'{n} {m} {rank} {len(free)}\n'+' '.join(map(str,rr))+'\n'+' '.join(map(str,cc))+'\n'+' '.join(map(str,free))+'\n'+'\n'.join(lines[1:])+'\n')
 return dict(rank=rank,rows=rr,columns=cc,free=free,minor_at_a2_e3=det)

def e_infinity_input(source,wd):
 with Path(source).open() as f:
  nr,nc,n=map(int,f.readline().split());lines=[f.readline().split() for _ in range(nr*nc)];assert int(f.readline())==1;choice=list(map(int,f.readline().split()));rr=choice[:n];cc=choice[n:]
 A=[[{} for _ in range(n)] for _ in range(n)]
 for i,r in enumerate(rr):
  for j,c in enumerate(cc):
   z=lines[r*nc+c];q=A[i][j]
   for h in range(int(z[0])):
    m=tuple(map(int,z[2+4*h:4+4*h]));q[m]=q.get(m,0)+Fraction(z[1+4*h])*3**int(z[4+4*h])
   A[i][j]={m:c for m,c in q.items() if c}
 degree=[[max([m[1] for m in q],default=-100000) for q in row] for row in A];upper,rw,cw=assignment([[-v for v in row] for row in degree])
 full=wd/'e-stage0';lead=wd/'e-stage0-leading'
 with full.open('w') as f,lead.open('w') as g:
  f.write(f'{n} {n}\n');g.write(f'{n} {n}\n')
  for i,row in enumerate(A):
   for j,q in enumerate(row):
    f.write(str(len(q)))
    l=[]
    for m,c in sorted(q.items(),reverse=True):
     exponent=rw[i]+cw[j]-m[1];assert exponent>=0;f.write(f' {c} {m[0]} 0 {exponent}')
     if not exponent:l.append((m[0],c))
    f.write('\n');g.write(str(len(l))+''.join(f' {c} {k} 0' for k,c in l)+'\n')
 return dict(upper=upper,row_weights=rw,column_weights=cw)

def verify(matrix,wd):
 start=time.time();wd=Path(wd);wd.mkdir(exist_ok=True,parents=True);choice=materialize_minor_input(matrix,wd/'minor');env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
 sources=[ROOT/'code/s61_binary_hyperplane_leading_kernel_20261102.cpp',ROOT/'code/s61_binary_hyperplane_kernel_peel_20261102.cpp',ROOT/'code/s61_binary_hyperplane_vertical_factor_export_20261102.cpp',ROOT/'code/s61_binary_hyperplane_constant_leading_peel_20261102.cpp']
 for i,src in enumerate(sources):compile_native(src,wd/f'native{i}')
 logs={}
 def run(name,cmd):
  cp=subprocess.run(list(map(str,cmd)),env=env,capture_output=True,text=True,check=True,timeout=180);logs[name]=cp.stderr;return cp
 def chain(stem):
  value=0;steps=[]
  for stage in range(20):
   full=wd/f'{stem}-stage{stage}';lead=wd/f'{stem}-stage{stage}-leading';inp=wd/f'{stem}-stage{stage}-input';pr=profile(lead,inp);k=36-pr['rank'];steps.append(pr)
   print('CHAIN',stem,stage,pr['rank'],value,flush=True)
   if not k:return dict(valuation=value,stages=steps,final_leading=str(lead),kernel_files=[str(wd/f'{stem}-kernel{i}') for i in range(stage)])
   ko=wd/f'{stem}-kernel{stage}';run(f'{stem}-kernel{stage}',[wd/'native0',inp,ko]);run(f'{stem}-peel{stage}',[wd/'native1',full,inp,ko,wd/f'{stem}-stage{stage+1}',wd/f'{stem}-stage{stage+1}-leading']);value+=k
  raise AssertionError('finite stage bound exceeded')
 ein=e_infinity_input(wd/'minor',wd);eproof=chain('e');assert ein['upper']==308 and eproof['valuation']==95
 # Certify the a-leading coefficient of the final e-leading determinant, and
 # all diagonal kernel denominators, are p-units. This yields LC_lex(e,a)(P3).
 lines=Path(eproof['final_leading']).read_text().splitlines();n,m=map(int,lines[0].split());polys=[z.split() for z in lines[1:]];degree=[]
 for z in polys:
  assert all(int(z[3+3*k])==0 for k in range(int(z[0])));degree.append(max([int(z[2+3*k]) for k in range(int(z[0]))],default=-100000))
 aupper,rw,cw=assignment([[-degree[i*n+j] for j in range(n)] for i in range(n)])
 with (wd/'lc-input').open('w') as f:
  f.write(str(n)+'\n')
  for i in range(n):
   for j in range(n):
    z=polys[i*n+j];f.write(z[0])
    for k in range(int(z[0])):
     ex=rw[i]+cw[j]-int(z[2+3*k]);assert ex>=0;f.write(' '+z[1+3*k]+' '+str(ex))
    f.write('\n')
 den=[]
 for path in eproof['kernel_files']:
  ls=Path(path).read_text().splitlines();idx=0
  while idx<len(ls):
   assert ls[idx].startswith('KERNEL ');z=ls[idx+1].split();terms=[]
   for k in range(int(z[0])):
    q=Fraction(z[1+3*k]);ex=int(z[2+3*k]);assert int(z[3+3*k])==0 and q.denominator%P;terms.append((ex,q))
   ex,q=max(terms);assert mod(q);den.append(dict(degree=ex,leading_coefficient=str(q)));idx+=n+2
 run('a-leading',[wd/'native3',wd/'lc-input',wd/'lc-proof',90]);lcproof=(wd/'lc-proof').read_text();last=[z for z in lcproof.splitlines() if z.startswith('PASS LEADING_P_UNIT ')];assert len(last)==1;z=last[0].split();residue=int(z[2]);valuation=int(z[4]);assert aupper==156 and valuation==85 and residue==11837
 factors={}
 for factor,target in [('a',37),('F',10),('C',26)]:
  run(f'export-{factor}',[wd/'native2',wd/'minor',factor,wd/f'{factor}-stage0',wd/f'{factor}-stage0-leading']);v=chain(factor);assert v['valuation']==target;factors[factor]=dict(valuation=target,stages=v['stages'])
 note=ROOT/'research/notes/S61_BINARY_HYPERPLANE_CUBIC_VERTICAL_20261102.md'
 paths=[Path(__file__),note,ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',ROOT/'code/certificates/s61_binary_hyperplane_cubic_first_minor_20261102.json',ROOT/'code/s61_binary_deeper_cubic_field_20261101.hpp',*sources]
 return dict(id='S61_BINARY_HYPERPLANE_CUBIC_VERTICAL_20261102',all_passed=True,prime=P,schur_matrix_sha256=digest(matrix),specialized_d=3,e_infinity=ein,e_valuation=eproof['valuation'],e_degree=ein['upper']-eproof['valuation'],e_stages=eproof['stages'],final_a_leading_assignment_bound=aupper,final_a_leading_valuation=valuation,final_a_leading_residue=residue,kernel_denominator_leading_coefficients=den,lex_e_a_leading_coefficient_is_p_unit=True,known_factor_valuations=factors,raw_matrix_coefficients_p_integral=True,leading_monomial_a_e=[aupper-valuation-sum(v['degree'] for v in den),ein['upper']-eproof['valuation']],proof_note=str(note.relative_to(ROOT)),native_logs=logs,input_sha256={str(p.relative_to(ROOT)):digest(p) for p in paths},limitations=['This is an exact certificate for one specified maximal minor specialized at d=3.','No cubic codimension conclusion without the separate modular support and horizontal certificates.'],seconds=time.time()-start)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--matrix',type=Path);args=ap.parse_args();assert not args.output.exists()
 with tempfile.TemporaryDirectory(prefix='poisson-hyperplane-vertical-') as td:
  wd=Path(td);matrix=args.matrix
  if matrix is None:
   build=build_matrix(wd/'matrix');matrix=wd/'matrix/schur'
  result=verify(matrix,wd/'vertical')
  if args.matrix is None:result['matrix_reconstruction']=build
 with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k] for k in ['all_passed','e_degree','lex_e_a_leading_coefficient_is_p_unit','seconds']}))
