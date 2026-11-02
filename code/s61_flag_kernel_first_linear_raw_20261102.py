"""Exact QQ residue module M/(v0,v1)M and its h,t actions, via raw135 six-jets."""
import sys,json,subprocess,tempfile,time
from pathlib import Path
from fractions import Fraction
from itertools import combinations_with_replacement as cwr
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize

def raw_first_linear(destination):
 started=time.monotonic();root=Path(__file__).resolve().parents[1]
 initial=json.loads((root/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json').read_text())
 point=s.Matrix(initial['point61']).applyfunc(s.Rational);model=S61()
 one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
 local=calculate(model,True,point);J=s.Matrix(local['jacobian']).applyfunc(s.Rational);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
 _,eq=specialize(model,dict(zip(one,[point[i]for i in one])),True)
 chart=local['chart_columns'];fixed=next(i for i in range(45)if i not in chart)
 terms=[]
 for f in eq:
  row=[]
  for mon,c in f.items():
   if fixed not in mon:
    a,b=[chart.index(i)for i in mon]+[44]*(2-len(mon));row.append((a,b,c))
  terms.append(row)
 piv=list(J.to_DM().rref()[1]);selected=list(J[:,piv].T.to_DM().rref()[1]);U=s.zeros(45,36);inverse=J[selected,piv].inv(method='DM')
 for i,col in enumerate(piv):U[col,:]=inverse[i,:]
 assert J*T==s.zeros(135,9)and J[selected,:]*U==s.eye(36)
 scales=[s.ilcm(*[s.denom(c)for c in list(J[i,:])+[c for a,b,c in terms[i]]])for i in range(135)]
 Jn=s.diag(*scales)*J;Un=U*s.diag(*[1/s.Integer(scales[i])for i in selected]);den=s.ilcm(*[s.denom(c)for c in Un]);Un=(den*Un).applyfunc(s.cancel)
 td=s.ilcm(*[s.denom(c)for c in T]);Tn=td*T
 ms={d:[m for m in cwr(range(9),d)if sum(i in(0,1)for i in m)<=1]for d in range(8)}
 ix={d:{m:i for i,m in enumerate(mm)}for d,mm in ms.items()};mons=sum((ms[d]for d in range(8)),[]);index={m:i for i,m in enumerate(mons)}
 assert len(mons)==6864;prime=32003
 with tempfile.TemporaryDirectory(prefix='zero12-first-kernel-action-')as tmp:
  temp=Path(tmp)
  binaries=[]
  for src in [root/'code/monomial_quotient_implicit_jets_20261102.cpp']:
   exe=temp/src.stem;subprocess.run(['c++','-O3','-w','-std=c++17','-I/opt/homebrew/include',str(src),'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(exe)],check=True);binaries.append(exe)
  inp=temp/'raw.txt'
  with inp.open('w')as out:
   out.write(f'135 9 7 36\n{den} {td}\n'+' '.join(str(len(ms[d]))for d in range(8))+'\n'+' '.join(map(str,selected))+'\n')
   for M in [Jn,Tn,Un]:
    for row in M.tolist():out.write(' '.join(map(str,row))+'\n')
   for i,row in enumerate(terms):out.write(str(len(row))+' '+' '.join(f'{a} {b} {c*scales[i]}'for a,b,c in row)+'\n')
   for a in range(1,7):
    for b in range(1,8-a):out.write(' '.join(str(ix[a+b].get(tuple(sorted(x+y)),-1))for x in ms[a]for y in ms[b])+'\n')
  qq=temp/'qq.txt';run=subprocess.run([str(binaries[0]),str(inp),str(qq)],capture_output=True,text=True,check=True,timeout=60);print(run.stderr,flush=True)
  rows=[line.split()for line in qq.read_text().splitlines()];assert len(rows)==136 and all(len(r)==6854 for r in rows)
  import shutil
  shutil.copyfile(qq,destination)
  print('RAW7 BASE1 SAVED',time.monotonic()-started,flush=True)
  return dict(all_passed=True,seconds=time.monotonic()-started,log=run.stderr)
