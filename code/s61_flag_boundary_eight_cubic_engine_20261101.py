#!/usr/bin/env python3
"""Exact rational eight-variable flag calculation: cubic stage.

Internal engine of s61_flag_boundary_eight_cubic_20261101.py.
All temporary products are rebuilt in the supplied working directory.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]); WD.mkdir(parents=True,exist_ok=True)
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.rings import ring

start=time.time()
cache=pickle.load(open(WD / 'flag-b1-universal-quadratic-cache.pkl','rb'))
q=json.load(open(WD / 'flag-b1-universal-quadratic.json'))
variables=list(map(s.Symbol,cache['variables']))
R,*_=ring(variables,s.QQ); zero=R.zero
dec=lambda p:R.from_dict({tuple(e):s.Rational(c) for e,c in p})
mat=lambda key:[[dec(v) for v in row] for row in cache[key]]
A2,W2,T,U,J=map(mat,['A2','W2','T','U','J'])
enc=lambda M:[[[[list(e),str(c)] for e,c in v.items()] for v in row] for row in M]
poly=lambda e:R.from_expr(s.cancel(e))
ud=poly(s.sympify(cache['ud']))
rows=cache['rows']; chosen=q['pivot_rows']; piv=q['pivot_columns']
mons2=list(cwr(range(9),2)); active={4,5,7,8}; radical=[0,1,2,3,6]
qr=s.Matrix(q['quadratic_rref']).applyfunc(s.sympify)
qd=s.factor(s.lcm([s.denom(v) for v in qr])); QD=poly(qd)
replace={mons2[col]:-poly(qd*qr[i,44]) for i,col in enumerate(piv)}
replace[(8,8)]=QD
G=s.Matrix([[s.cancel(replace[tuple(sorted((i,j)))].as_expr()/qd) for j in sorted(active)] for i in sorted(active)])
assert G.subs({variables[0]:1}).rank()>=2
# The degree-two quotient has one active socle coordinate; rank >=2 of its
# multiplication form implies every all-active cubic vanishes over QQ(b2).
mons3=[m for m in cwr(range(9),3) if sum(i in active for i in m)<=1]
mons3 += [tuple(sorted((i,8,8))) for i in radical]
idx3={m:i for i,m in enumerate(mons3)}
assert len(mons3)==100
def project3(m):
    aa=tuple(i for i in m if i in active)
    if len(aa)>=3:return None,zero
    if len(aa)==2:
        rr=next(i for i in m if i not in active)
        return idx3[tuple(sorted((rr,8,8)))],replace[aa]
    return idx3[m],QD
Wred=[[zero for _ in mons2] for _ in W2]
for i,row in enumerate(W2):
    for k,v in enumerate(row):
        if not v:continue
        mon=mons2[k]
        if all(j in active for j in mon):Wred[i][44]+=replace[mon]*v
        else:Wred[i][k]+=QD*v
print('quadratic projection',time.time()-start,flush=True)
Q3=[[zero for _ in mons3] for _ in range(135)]
fixed=cache['fixed']; chart=cache['chart']
for row,f in enumerate(cache['F']):
    for mon,c in f:
        if fixed in mon:continue
        i,j=[chart.index(k) for k in mon]+[44]*(2-len(mon))
        co=poly(s.sympify(c))
        for a in range(9):
            for k,mm in enumerate(mons2):
                val=T[i][a]*Wred[j][k]+T[j][a]*Wred[i][k]
                if not val:continue
                col,mul=project3(tuple(sorted((a,)+mm)))
                if mul:Q3[row][col]+=co*mul*val
print('raw projected cubics',time.time()-start,flush=True)
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
print('all projected implicit cubics',time.time()-start,flush=True)
pickle.dump(enc(A3),open(WD / 'flag-b1-projected-cubic-intermediate.pkl','wb'))
M=s.Matrix([[A2[i][j].as_expr() for j in piv] for i in chosen])
iv=M.to_DM().to_field().inv().to_Matrix()
idn=s.factor(s.lcm([s.denom(v) for v in iv])); ID=poly(idn)
IV=[[poly(idn*v) for v in row] for row in iv.tolist()]
B=[[zero for _ in chosen] for _ in range(135)]
for i in range(135):
    for j in range(9):
        for k,col in enumerate(piv):B[i][j]+=A2[i][col]*IV[k][j]
print('quadratic row relations',time.time()-start,flush=True)
C=[[ID*v for v in row] for row in A3]
for i,row in enumerate(B):
    for j,co in enumerate(row):
        if co:
            for k,v in enumerate(A3[chosen[j]]):
                if v:C[i][k]-=co*v
print('true cubic classes',time.time()-start,flush=True)
# Common scalar is irrelevant over the function field; remove row contents
# before any row-space computation.
primitive=[]
for row in C:
    nonzero=[v for v in row if v]
    if not nonzero:primitive.append(row);continue
    g=nonzero[0]
    for v in nonzero[1:]:g=g.gcd(v)
    primitive.append([v.exquo(g) for v in row])
cc=s.Matrix([[v.as_expr() for v in row] for row in primitive])
sample={v:0 for v in variables};sample[variables[0]]=1
rr=list(cc.subs(sample).T.to_DM().rref()[1]); print('sample rank',len(rr),rr,flush=True)
cr,cp=cc[rr,:].to_DM().rref(); cr=cr.to_Matrix()[:len(rr),:]
assert (cc-cc[:,cp]*cr).applyfunc(s.cancel)==s.zeros(135,100)
symbols=sorted(set().union(*(v.free_symbols for v in cr)),key=str)
print('cubic rref symbols',symbols,'time',time.time()-start,flush=True)
pickle.dump(dict(C=enc(primitive),monomials=mons3,variables=cache['variables']),open(WD / 'flag-b1-universal-cubic-cache.pkl','wb'))
Path(WD / 'flag-b1-universal-cubic.json').write_text(json.dumps(dict(
    rank=len(rr),pivot_rows=rr,pivot_columns=list(cp),monomials=mons3,
    cubic_rref=[[str(v) for v in row] for row in cr.tolist()],
    rref_symbols=list(map(str,symbols)),quadratic_denominator=str(qd),
    row_relation_denominator=str(idn),seconds=time.time()-start),indent=2)+'\n')
print('UNIVERSAL CUBIC PASS',time.time()-start,flush=True)
