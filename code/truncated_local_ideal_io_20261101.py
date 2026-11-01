"""Sparse input/output for work/truncated_polynomial_closure.cpp (exploration)."""
import json
from itertools import combinations_with_replacement
from pathlib import Path

def monomials(n,maximum_degree):
 return [m for d in range(maximum_degree+1) for m in combinations_with_replacement(range(n),d)]

def write_input(path,prime,n,jets,t_jets=None,maximum_degree=None):
 """jets[d]: generators x degree-d monomials in combinations_with_replacement order."""
 jets={int(k):v for k,v in jets.items()};t_jets={} if t_jets is None else {int(k):v for k,v in t_jets.items()};D=max(jets) if maximum_degree is None else maximum_degree
 mons=monomials(n,D);idx={m:i for i,m in enumerate(mons)};count=len(next(iter(jets.values())));F=[{} for _ in range(count)]
 for d,A in jets.items():
  if d>D:continue
  mm=list(combinations_with_replacement(range(n),d));assert len(A)==count
  for row,coefs in enumerate(A):
   assert len(coefs)==len(mm)
   for mon,c in zip(mm,coefs):
    if int(c)%prime:F[row][idx[mon]]=int(c)%prime
 with Path(path).open('w') as f:
  f.write(f'{prime} {n} {len(mons)} {count}\n');f.write(' '.join(str(len(m)) for m in mons)+'\n')
  for j in range(n):f.write(' '.join(str(idx[tuple(sorted(m+(j,)))]) if len(m)<D else '-1' for m in mons)+'\n')
  for row in F:f.write(str(len(row))+' '+' '.join(f'{i} {v}' for i,v in sorted(row.items()))+'\n')
  tt={}
  for d,coefs in t_jets.items():
   if d>D:continue
   for mon,v in zip(combinations_with_replacement(range(n),d),coefs):tt[idx[mon]]=int(v)%prime
  f.write(' '.join(str(tt.get(i,0)) for i in range(len(mons)))+'\n')
 return dict(prime=prime,variables=n,maximum_degree=D,columns=len(mons),generators=count)

def read_initial(path,n,maximum_degree):
 mons=monomials(n,maximum_degree);lines=Path(path).read_text().splitlines();assert int(lines[0])==len(lines)-1;polys=[]
 for line in lines[1:]:
  row=list(map(int,line.split()));assert len(row)==1+2*row[0];polys.append({mons[i]:c for i,c in zip(row[1::2],row[2::2])})
 return polys

def m2_initial(path,n,maximum_degree,prime):
 polys=read_initial(path,n,maximum_degree)
 expr=lambda f:'+'.join(str(c)+'*'+'*'.join('v_'+str(j) for j in mon) if mon else str(c) for mon,c in f.items()) or '0'
 return f'R=ZZ/{prime}[v_0..v_{n-1}];\nJ=ideal('+','.join(expr(f) for f in polys)+');\nprint("dim",dim(R/J));print("degree",degree(R/J));print("hilbert",hilbertSeries(R/J));\n'
