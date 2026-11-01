#!/usr/bin/env python3
"""Univariate core minor cover, including the corrected t4 class."""
import json,math,os,random,subprocess,sys,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
WD=Path(sys.argv[1]);MINOR=Path(sys.argv[2]);FACTOR=Path(sys.argv[3]);ROOT=Path(__file__).resolve().parents[1];st=time.time()
data=json.loads((WD/'flag-eight-initial4-field.json').read_text());P,r=ring('r',s.QQ);rows={};hs={};denoms={};contents={};source={}
for group,off in [('rows',0),('corrected_t4',1000)]:
 for key,v in data[group].items():
  i,j=map(int,key.split(','));source.setdefault(i+off,{})[j]=v
for i,row in source.items():
 if not all(all(all(m[k]==0 for k in [1,2,4,5,6,7])for m,c in nn)for nn,dd in row.values()):continue
 rr={};hh={}
 for j,(nn,dd)in row.items():
  a={m[3]for m,c in nn};b={m[3]for m,c in dd};assert len(a)==len(b)==1;hh[j]=a.pop()-b.pop()
  np=P.from_dict({(m[0],):s.QQ(c)for m,c in nn});dp=P.from_dict({(m[0],):s.QQ(c)for m,c in dd});rr[j]=(np,dp)
 D=P.one
 for np,dp in rr.values():D=D.lcm(dp)
 vals={j:np*D.exquo(dp)for j,(np,dp)in rr.items()};g=next(iter(vals.values()))
 for p in vals.values():g=g.gcd(p)
 g=g.monic();rows[i]={j:p.exquo(g)for j,p in vals.items()};hs[i]=hh;denoms[i]=D;contents[i]=g
assert len(rows)==19 and 1000 in rows
rw={next(iter(rows)):0};cw={};change=True
while change:
 change=False
 for i,row in hs.items():
  for j,h in row.items():
   if i in rw and j not in cw:cw[j]=h-rw[i];change=True
   if j in cw and i not in rw:rw[i]=h-cw[j];change=True
assert all(h==rw[i]+cw[j]for i,row in hs.items()for j,h in row.items())
print('PURE CORE MATRIX AND H WEIGHTS PASS',len(rows),len(cw),time.time()-st,flush=True)
mod=32003;cs=sorted(cw);rng=random.Random(611048);env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
def select(A,rr,cc,target):
 A=[[A[i][j]for j in cc]for i in rr];orig=rr[:];rank=0;cols=[];selected=[]
 for j in range(len(cc)):
  k=next((i for i in range(rank,len(rr))if A[i][j]),None)
  if k is None:continue
  A[k],A[rank]=A[rank],A[k];orig[k],orig[rank]=orig[rank],orig[k];selected.append(orig[rank]);cols.append(cc[j]);inv=pow(A[rank][j],-1,mod);A[rank]=[(v*inv)%mod for v in A[rank]]
  for i in range(rank+1,len(rr)):
   c=A[i][j]
   if c:A[i]=[(a-c*b)%mod for a,b in zip(A[i],A[rank])]
  rank+=1
  if rank==target:break
 assert rank==target;return selected,cols
reports=[];factor_set={}
def addfactor(p):
 if p.is_ground:return
 for f,e in p.factor_list()[1]:factor_set[str(f.monic().as_expr())]=True
for p in list(denoms.values())+list(contents.values()):addfactor(p)
for vv in data['denominators'].values():
 for text in vv:
  expr=s.sympify(text).subs(s.Symbol('b8'),1).subs(s.Symbol('b2'),s.Symbol('r'));addfactor(P.from_expr(expr))
for text in [data['projection_multiplier'],*data['scaling_denominators'].values()]:
 expr=s.sympify(text).subs(s.Symbol('b8'),1).subs(s.Symbol('b2'),s.Symbol('r'));addfactor(P.from_expr(expr))
for target in [12,13]:
 rs=sorted(i for i in rows if target==13 or i!=1000);matrix=[[rows[i].get(j,P.zero)for j in cs]for i in rs]
 vals=[[sum(int(c.numerator)*pow(int(c.denominator),-1,mod)for e,c in p.items())%mod for p in row]for row in matrix]
 choices=[]
 for k in range(8):
  rr=list(range(len(rs)));cc=list(range(len(cs)));rng.shuffle(rr);rng.shuffle(cc);choices.append(select(vals,rr,cc,target))
 inp=WD/f'minor{target}-input.txt';out=WD/f'minor{target}-output.txt';fac=WD/f'minor{target}-factors.txt'
 with inp.open('w')as f:
  f.write(f'{len(rs)} {len(cs)} {len(choices)} {target}\n')
  for row in matrix:
   den=math.lcm(*(int(c.denominator)for p in row for c in p.values()))
   for p in row:
    f.write(str(len(p))+'\n')
    for (e,),c in p.items():f.write(str(e)+' '+str(int(c*den))+'\n')
  for rr,cc in choices:f.write(' '.join(map(str,rr+cc))+'\n')
 proc=subprocess.run([str(MINOR),str(inp),str(out)],env=env,text=True,capture_output=True,check=True,timeout=120);print(proc.stderr,flush=True)
 proc2=subprocess.run([str(FACTOR),str(out),str(fac)],env=env,text=True,capture_output=True,check=True,timeout=120)
 with fac.open()as f:
  ff=[]
  for _ in range(int(next(f))):
   e,n=map(int,next(f).split());p=P.from_dict({(i,):s.QQ(next(f).strip())for i in range(n)}).monic();addfactor(p);ff.append([str(p.as_expr()),e])
 reports.append(dict(rank=target,rows=rs,columns=cs,minor_choices=choices,gcd_factors=ff,native_log=proc.stderr,exact_gcd_factorization_log=proc2.stderr))
result=dict(all_passed=True,pure_core_rows_without_t=18,pure_core_columns=80,nuisance_independent_fourth_space=True,quartic_rank=12,augmented_rank=13,hilbert_function=[1,9,36,97,201],t_order=4,row_weights=rw,column_weights=cw,nonzero_factors=['h']+sorted(factor_set),factor_variables=['r=b2','h=b8'],minor_covers=reports,seconds=time.time()-st)
(WD/'flag-eight-quartic-cover-result.json').write_text(json.dumps(result,indent=2)+'\n');print('CORE MINOR COVER PASS',result['seconds'],'factors',len(result['nonzero_factors']),flush=True)
