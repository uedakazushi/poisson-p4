#!/usr/bin/env python3
"""Exact ordinary-flag fourth-order core cover: active3."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1])
import json,os,subprocess,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
import numpy as np
st=time.time();data=json.load(open(WD/'flag-global-low-quad.json'));r=s.Symbol('b2');A=s.Matrix(data['basis']).applyfunc(s.sympify);m2=list(cwr(range(4),2));m3=list(cwr(range(4),3));B=s.zeros(36,20)
for i in range(9):
 for j in range(4):
  for k,m in enumerate(m2):B[4*i+j,m3.index(tuple(sorted(m+(j,))))]+=A[i,k]
g=s.expand((r+s.Rational(119,20))**2*(r+s.Rational(91,400)));degree=max(s.Poly(p,r).degree()for p in B);terms={(i,j):s.Poly(B[i,j],r).terms()for i in range(36)for j in range(20)if B[i,j]};mod=32003;env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
for bound in range(11):
 n=36*(bound+1);nr=20*(degree+bound+1);dd={}
 for(i,j),pol in terms.items():
  for(e,),c in pol:
   for k in range(bound+1):dd[j*(degree+bound+1)+e+k,i*(bound+1)+k]=c
 for j in range(20):
  for(e,),c in s.Poly(g,r).terms():dd[j*(degree+bound+1)+e,n+j]=c
 M=np.zeros((nr,n+20),dtype=np.int64)
 for(i,j),c in dd.items():M[i,j]=int(c.p)*pow(int(c.q),-1,mod)%mod
 rank=0
 for j in range(n):
  ix=np.flatnonzero(M[rank:,j])
  if not len(ix):continue
  k=rank+int(ix[0]);M[[rank,k]]=M[[k,rank]];M[rank]=M[rank]*pow(int(M[rank,j]),-1,mod)%mod
  ix=np.flatnonzero(M[rank+1:,j])+rank+1
  if len(ix):M[ix]=(M[ix]-M[ix,j,None]*M[rank,None,:])%mod
  rank+=1
 consistent=not np.any(M[rank:,n:]);print('bound',bound,'matrix',nr,n,'rank',rank,'consistent',consistent,time.time()-st,flush=True)
 if not consistent:continue
 inp=WD/'active-cubic-linear-input.txt';out=WD/'active-cubic-linear-output.txt'
 with inp.open('w')as f:
  f.write(f'{nr} {n} 20\n')
  for i in range(nr):f.write(' '.join(str(dd.get((i,j),0))for j in range(n+20))+'\n')
 subprocess.run([sys.argv[2],str(inp),str(out)],env=env,check=True,timeout=120)
 with out.open()as f:
  nn,kk,rank,good=map(int,next(f).split());assert(nn,kk)==(n,20)
  if not good:continue
  X=s.zeros(20,36)
  for i in range(n):
   q,e=divmod(i,bound+1)
   for j in range(20):X[j,q]+=s.Rational(next(f).strip())*r**e
 assert (X*B-g*s.eye(20)).applyfunc(s.expand)==s.zeros(20)
 enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
 Path(WD/'flag-global-low-active-cubics.json').write_text(json.dumps(dict(lift_numerator=enc(X),denominator=str(g),polynomial_degree=bound,source_basis='primitive nine quadratic rows times four active variables',seconds=time.time()-st),indent=2)+'\n');print('LOW DEGREE ALL ACTIVE CUBIC LIFT PASS',time.time()-st,flush=True);break
else:raise RuntimeError('no inverse through bound10')
