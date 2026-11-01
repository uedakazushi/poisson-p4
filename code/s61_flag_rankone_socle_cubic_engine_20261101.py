#!/usr/bin/env python3
"""Exact cubic classes at the rank-one socle core; internal isolated engine."""
import sys
from pathlib import Path
WD=Path(sys.argv[1])
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.rings import ring
start=time.time();cache=pickle.load(open(WD / 'flag-b1-universal-quadratic-cache.pkl','rb'));r0=-s.Rational(91,400)
variables=list(map(s.Symbol,cache['variables'][1:]));R,*_=ring(variables,s.QQ);zero=R.zero
def dec(p):
    ans={}
    for ex,c in p:
        key=tuple(ex[1:]);ans[key]=ans.get(key,s.Rational(0))+s.Rational(c)*r0**ex[0]
    return R.from_dict(ans)
mat=lambda key:[[dec(v) for v in row] for row in cache[key]]
A2,W2,T,U,J=map(mat,['A2','W2','T','U','J']);poly=lambda e:R.from_expr(s.cancel(s.sympify(e).subs(s.Symbol('b2'),r0)))
ud=poly(cache['ud']);rows=cache['rows'];fixed=cache['fixed'];chart=cache['chart']
mons2=list(cwr(range(9),2));mons3=list(cwr(range(9),3));ix3={m:i for i,m in enumerate(mons3)}
Q3=[[zero for _ in mons3] for _ in range(135)]
for row,f in enumerate(cache['F']):
    for mon,c in f:
        if fixed in mon:continue
        i,j=[chart.index(k) for k in mon]+[44]*(2-len(mon));co=poly(c)
        for a in range(9):
            for k,mm in enumerate(mons2):
                val=T[i][a]*W2[j][k]+T[j][a]*W2[i][k]
                if val:Q3[row][ix3[tuple(sorted((a,)+mm))]]+=co*val
print('raw3',time.time()-start,flush=True)
W3=[[zero for _ in mons3] for _ in range(45)]
for i,row in enumerate(U):
    for j,co in enumerate(row):
        if co:
            for k,v in enumerate(Q3[rows[j]]):
                if v:W3[i][k]-=co*v
A3=[[ud*v for v in row] for row in Q3]
for i,row in enumerate(J):
    for j,co in enumerate(row):
        if co:
            for k,v in enumerate(W3[j]):
                if v:A3[i][k]+=co*v
print('implicit3',time.time()-start,flush=True)
qq=json.load(open(WD / 'flag-b1-universal-quadratic.json'));chosen=qq['pivot_rows']
M=s.Matrix([[v.as_expr() for v in A2[i]] for i in chosen]);qr,piv=M.to_DM().rref();qr=qr.to_Matrix();iv=M[:,piv].inv();assert not set().union(*(v.free_symbols for v in iv))
B=s.Matrix([[A2[i][col].as_expr() for col in piv] for i in range(135)])*iv
C=[row[:] for row in A3]
for i in range(135):
    for j in range(9):
        co=poly(B[i,j])
        if co:
            for k,v in enumerate(A3[chosen[j]]):
                if v:C[i][k]-=co*v
mult=s.zeros(81,165)
for i in range(9):
    for a in range(9):
        for k,m in enumerate(mons2):mult[9*i+a,ix3[tuple(sorted((a,)+m))]]+=qr[i,k]
mr,mp=mult.to_DM().rref();mr=mr.to_Matrix()[:len(mp),:];assert len(mp)==64
for row in C:
    for i,p in enumerate(mp):
        v=row[p]
        if v:
            for j in range(165):
                if mr[i,j]:row[j]-=v*s.Rational(mr[i,j])
print('true cubics modulo64 products',time.time()-start,flush=True)
primitive=[];contents=[]
for row in C:
    nz=[v for v in row if v]
    if not nz:primitive.append(row);contents.append('0');continue
    g=nz[0]
    for v in nz[1:]:g=g.gcd(v)
    primitive.append([v.exquo(g) for v in row]);contents.append(str(g.as_expr()))
cc=s.Matrix([[v.as_expr() for v in row] for row in primitive]);sample={v:i+1 for i,v in enumerate(variables)};rr=list(cc.subs(sample).T.to_DM().rref()[1]);cr,cp=cc[rr,:].to_DM().rref();cr=cr.to_Matrix()[:len(rr),:]
assert (cc-cc[:,cp]*cr).applyfunc(s.cancel)==s.zeros(135,165)
t3=list(W3[44])
for i,pivot in enumerate(mp):
    value=t3[pivot]
    if value:
        for j in range(165):
            if mr[i,j]:t3[j]-=value*s.Rational(mr[i,j])
t3expr=s.Matrix([[v.as_expr() for v in t3]])
t3res=(t3expr-t3expr[:,cp]*cr).applyfunc(s.cancel)
t_order_at_least_four=(t3res==s.zeros(1,165))
print('t order at least four',t_order_at_least_four,flush=True)
out=dict(t_order_at_least_four=t_order_at_least_four,core=str(r0),quadratic_rank=9,quadratic_linear_multiple_rank=64,cubic_rank=len(rr),cubic_rref=[[[list(mons3[j]),str(v)] for j,v in enumerate(row) if v] for row in cr.tolist()],rref_symbols=list(map(str,sorted(set().union(*(v.free_symbols for v in cr)),key=str))),contents=contents,pivot_rows=rr,pivot_columns=list(cp),seconds=time.time()-start)
Path(WD / 'flag-rankone-socle-cubic.json').write_text(json.dumps(out,indent=2)+'\n');print('SPECIAL CORE CUBIC PASS',len(rr),out['rref_symbols'],time.time()-start,flush=True)
