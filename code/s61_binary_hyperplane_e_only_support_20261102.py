#!/usr/bin/env python3
"""Exclude unknown e-only divisors of the B-free quartic prolongation module."""
import argparse,functools,hashlib,json,math,os,subprocess,tempfile,time
from fractions import Fraction
from pathlib import Path
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_HYPERPLANE_E_ONLY_SUPPORT_20261102'
BASE='reports/S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102.json'
CERT='code/certificates/s61_binary_hyperplane_e_only_selection_20261102.json'
SCHUR='code/s61_binary_hyperplane_e_only_schur_20261102.cpp'
MODULE='code/s61_binary_hyperplane_vertical_rowmodule_20261102.cpp'
NOTE=f'research/notes/{STEM}.md'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
@functools.lru_cache(None)
def rational(s):return Fraction(s)
def run_certificate(core,wd):
 from s61_binary_hyperplane_Bfree_horizontal_20261102 import compile_native
 st=time.monotonic();core=Path(core);wd=Path(wd);wd.mkdir(exist_ok=True,parents=True)
 base=json.loads((ROOT/BASE).read_text());key='hyperplane-Bfree-t-targets-exact-core-20261102.in';assert sha(core)==base['generated_sha256'][key]
 s=json.loads((ROOT/CERT).read_text());rows=s['rows'];cols=s['columns'];rb=s['row_bounds'];cb=s['column_bounds'];ranks=s['ranks'];prec=len(ranks)-1;n=238;p=32003
 assert len(rows)==len(set(rows))==n and min(rows)>=0 and max(rows)<321 and sorted(cols)==list(range(n))and sum(ranks)==n
 assert sum(rb)+sum(cb)==s['degree_upper']==1132
 ri={x:i for i,x in enumerate(rows)};ci={x:i for i,x in enumerate(cols)};selected=[[None]*n for i in rows];Ms={v:[{}for _ in range(321)] for v in (1,2)}
 @functools.lru_cache(None)
 def mod(s):
  q=rational(s);assert q.denominator%p;return q.numerator*pow(q.denominator,-1,p)%p
 count=0
 with core.open()as f:
  assert list(map(int,f.readline().split()))==[328,238,321]
  for i in range(328):
   for j in range(n):
    l=f.readline();z=l.split();nt=int(z[0]);assert len(z)==1+4*nt
    if i>=321:continue
    w1,w2={},{};val={}
    for c,aa,ee,dd in zip(z[1::4],z[2::4],z[3::4],z[4::4]):
     a,e,d=int(aa),int(ee),int(dd);v=mod(c)*pow(3,d,p)%p;count+=1
     w1[e]=(w1.get(e,0)+v)%p;w2[e]=(w2.get(e,0)+v*pow(2,a,p))%p
     if i in ri:val[e]=val.get(e,0)+rational(c)*3**d
    for t,w in [(1,w1),(2,w2)]:
     w={k:x for k,x in w.items()if x}
     if w:Ms[t][i][j]=[w.get(k,0)for k in range(max(w)+1)]
    if i in ri:
     u,v=ri[i],ci[j];val={k:x for k,x in val.items()if x};assert all(k<=rb[u]+cb[v]for k in val)
     selected[u][v]={rb[u]+cb[v]-k:x for k,x in val.items()if rb[u]+cb[v]-k<=prec}
  assert not f.read().strip()
 env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
 inp=wd/'schur.in';out=wd/'schur.out';scalar=1
 with inp.open('w')as f:
  f.write(f'{n} {prec} {len(ranks)-1}\n'+' '.join(map(str,ranks[:-1]))+'\n')
  for row in selected:
   den=math.lcm(*(v.denominator for q in row for v in q.values()));g=math.gcd(*(int(v*den)for q in row for v in q.values()));assert g%p and den%p
   scalar=scalar*den*pow(g,-1,p)%p
   for q in row:f.write(str(len(q))+' '+' '.join(f'{k} {int(v*den)//g}'for k,v in sorted(q.items()))+'\n')
 compile_native(ROOT/SCHUR,wd/'schur');compile_native(ROOT/MODULE,wd/'module')
 def native(args):
  r=subprocess.run(list(map(str,args)),env=env,check=True,capture_output=True,text=True);return r.stdout+r.stderr
 slog=native([wd/'schur',inp,out]);lines=out.read_text().splitlines();assert len(lines)==len(ranks)
 stages=[]
 for line,r in zip(lines[:-1],ranks):
  q=line.split();assert q[0]=='STAGE' and int(q[3])==r and int(q[7])%p;stages.append(dict(size=int(q[1]),rank=r,valuation=int(q[5]),pivot_unit=int(q[7])))
 z=lines[-1].split();assert z[0]=='TOTAL';valuation,last,lc=map(int,z[1::2]);assert valuation==sum(b['size']-b['rank']for b in stages)==190 and last%p
 sign=(-1)**sum(cols[i]>cols[j]for i in range(n)for j in range(i+1,n));origlc=lc*pow(scalar,-1,p)*sign%p;assert origlc==s['predicted_LC_mod']==2514
 e=sp.Symbol('e');pols=[];modules={};gen={inp.name:sha(inp),out.name:sha(out)}
 for a,M in Ms.items():
  v=wd/f'rows_a{a}.in';w=wd/f'rows_a{a}.out'
  with v.open('w')as f:
   f.write(f'{p} 321 238 0\n0 1\n')
   for row in M:
    f.write(str(len(row))+'\n')
    for j,pol in row.items():f.write(f'{j} {len(pol)} '+' '.join(map(str,pol))+'\n')
  log=native([wd/'module',v,w]);ll=w.read_text().splitlines();z=next(l.split()for l in ll if l.startswith('ROWMODULE'));assert int(z[1])==n;deg=int(z[3])
  z=next(l.split()for l in ll if l.startswith('DETERMINANT_COEFFICIENTS'));vco=list(map(int,z[2:]));assert int(z[1])==len(vco)
  pol=sp.Poly.from_list(list(reversed(vco)),e,modulus=p).monic();assert pol.degree()==deg==522
  pols.append(pol);modules[str(a)]={'rank':n,'degree':deg,'monic_coefficients_ascending':[int(pol.nth(i))%p for i in range(deg+1)]};gen.update({v.name:sha(v),w.name:sha(w)})
 g=sp.gcd(*pols);assert g.degree()==0
 paths=[Path(__file__).relative_to(ROOT),BASE,CERT,SCHUR,MODULE,NOTE,'code/s61_binary_hyperplane_Bfree_horizontal_20261102.py'];hs={str(q):sha(ROOT/q)for q in paths}
 return dict(all_passed=True,scope='No irreducible polynomial of e alone divides all maximal minors of the three-core 321x238 B-free quartic prolongation. This is independent of the two-cut vertical patch; no assertion of rank at every core point or ambient removal.',static_input_sha256=hs,input_sha256=hs,generated_core_sha256=sha(core),generated_outputs_sha256=gen,coefficient_terms_checked=count,rows=n,prime=p,selected_degree_bound=1132,reciprocal_precision=prec,QQ_Schur_stages=stages,reciprocal_valuation=valuation,last_unit=last,row_scale_mod=scalar,column_sign=sign,original_degree=1132-valuation,leading_coefficient_mod=origlc,modular_modules=modules,modular_gcd_degree=g.degree(),seconds=time.monotonic()-st)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--core',type=Path,required=True,help='Generated matrix with SHA fixed by BFREE_HORIZONTAL');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workdir',type=Path);a=ap.parse_args()
 assert not a.output.exists()
 if a.workdir:r=run_certificate(a.core,a.workdir)
 else:
  with tempfile.TemporaryDirectory(prefix='poisson-Bfree-e-only-')as td:r=run_certificate(a.core,td)
 with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print(f"PASS degree {r['original_degree']}, leading mod {r['leading_coefficient_mod']}, gcd {r['modular_gcd_degree']}; {r['seconds']:.3f} seconds")
