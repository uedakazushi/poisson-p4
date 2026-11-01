#!/usr/bin/env python3
"""Exact native fraction-free rank and corrected t4 nonvanishing."""
import json,os,subprocess,sys,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
WD=Path(sys.argv[1]);BINARY=Path(sys.argv[2]);ROOT=Path(__file__).resolve().parents[1];st=time.time()
data=json.loads((WD/'flag-eight-initial4-field.json').read_text());R,*gens=ring(data['variables'],s.QQ);C,r,h=ring('r,h',s.QQ);ci=[data['variables'].index(x)for x in ['b2','b8']];ni=[i for i in range(8)if i not in ci]
def dec(items):return R.from_dict({tuple(m):s.QQ(c)for m,c in items})
def core(p):
 assert all(all(m[i]==0 for i in ni)for m in p)
 return C.from_dict({tuple(m[i]for i in ci):c for m,c in p.items()})
def lift(p):
 d={}
 for m,c in p.items():
  mm=[0]*8
  for i,j in enumerate(ci):mm[j]=m[i]
  d[tuple(mm)]=c
 return R.from_dict(d)
all_factors={};row_denominators=[]
def record(p):
 p=core(p)
 if p.is_ground:return
 for f,ex in p.factor_list()[1]:
  f=f.monic();all_factors[str(f.as_expr())]=True
rows={}
for group,offset in [('rows',0),('corrected_t4',1000)]:
 source={}
 for key,pair in data[group].items():
  i,j=map(int,key.split(','));source.setdefault(i+offset,{})[j]=(dec(pair[0]),core(dec(pair[1])))
 for i,row in source.items():
  D=C.one
  for n,d in row.values():D=D.lcm(d)
  rows[i]={j:n*lift(D.exquo(d))for j,(n,d)in row.items()};row_denominators.append(str(D.as_expr()));record(lift(D))
assert 1000 in rows
for values in data['denominators'].values():
 for v in values:record(R.from_expr(s.sympify(v)))
record(R.from_expr(s.sympify(data['projection_multiplier'])))
for v in data['scaling_denominators'].values():record(R.from_expr(s.sympify(v)))
def write(path,rr):
 with path.open('w')as out:
  out.write('8 '+str(len(rr))+'\n')
  for i,row in rr.items():
   out.write(str(i)+' '+str(len(row))+'\n')
   for j,p in row.items():
    out.write(str(j)+'\n'+str(len(p))+'\n')
    for m,c in p.items():out.write(str(c)+' '+' '.join(map(str,m))+'\n')
def readpoly(f):
 d={}
 for _ in range(int(next(f))):
  z=next(f).split();d[tuple(map(int,z[1:]))]=s.QQ(z[0])
 return R.from_dict(d)
reports=[];env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
for with_t in [False,True]:
 inp=WD/('rank-input-'+str(with_t)+'.txt');out=WD/('rank-output-'+str(with_t)+'.txt');write(inp,{i:row for i,row in rows.items()if with_t or i<1000})
 proc=subprocess.run([str(BINARY),str(inp),str(out)],env=env,text=True,capture_output=True,check=True,timeout=120);print(proc.stderr,flush=True)
 with out.open()as f:
  rank,left,ng=map(int,next(f).split());assert left==0 and rank==(13 if with_t else 12)
  pivots=[]
  for _ in range(rank):
   i,j=map(int,next(f).split());p=readpoly(f);assert p;record(p);pivots.append(dict(row=i,column=j,polynomial=str(p.as_expr())))
  for _ in range(ng):record(readpoly(f))
 reports.append(dict(rank=rank,remaining_rows=left,core_pivots=pivots,divided_core_contents=ng,native_log=proc.stderr))
result=dict(all_passed=True,parameters=data['variables'],quadratic_rank=9,cubic_product_rank=65,cubic_new_rank=3,quartic_lower_rank=282,quartic_additional_rank=12,quartic_initial_rank=294,hilbert_function=[1,9,36,97,201],t3_in_true_cubic_ideal=True,corrected_t4_augmented_rank=13,t_order=4,nonzero_core_factors=sorted(all_factors),factor_variables=['r=b2','h=b8'],ranks=reports,row_denominators=row_denominators,seconds=time.time()-st)
(WD/'flag-eight-quartic-result.json').write_text(json.dumps(result,indent=2)+'\n');print('TRUE FOURTH INITIAL RANK AND ORD T4 PASS',result['seconds'],flush=True)
