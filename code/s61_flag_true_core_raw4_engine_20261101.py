#!/usr/bin/env python3
"""Rebuild true rank-drop core raw four-jets with exact FLINT arithmetic.
Internal engine; supplied working directory and native executable are mandatory.
"""
import json,pickle,sys,tempfile,subprocess,os,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
ROOT=Path.cwd()
WD=Path(sys.argv[1]);BINARY=Path(sys.argv[2]);start=time.time();cache=pickle.load(open(WD/'flag-b1-universal-quadratic-cache.pkl','rb'));general=True;D=4
assert not s.sympify(cache['ud']).free_symbols and not s.sympify(cache['td']).free_symbols
oldvars=list(map(s.Symbol,cache['variables']));variables=oldvars if general else oldvars[1:];np=len(variables);nv=np+10;r0=-s.Rational(91,400)
label='flag-true-core-fourth'
def specialize(poly):
    if general:return [[tuple(e),s.Rational(c)] for e,c in poly]
    out={}
    for e,c in poly:
        key=tuple(e[1:]);out[key]=out.get(key,s.Rational(0))+s.Rational(c)*r0**e[0]
    return [(e,c) for e,c in out.items() if c]
def expr_poly(expr):
    expr=s.sympify(expr)
    if not general:expr=expr.subs(oldvars[0],r0)
    return [(m,c) for m,c in s.Poly(expr,*variables).terms() if c]
def writepoly(out,poly):
    out.write(str(len(poly))+'\n')
    for e,c in poly:out.write(str(c)+' '+' '.join(map(str,tuple(e)+(0,)*10))+'\n')
env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
with tempfile.TemporaryDirectory(prefix='root-parametric-jets-') as temp:
    inp=Path(temp)/'input.txt';output=Path(temp)/'output.txt'
    with inp.open('w') as out:
        out.write(f'135 {np} 10 {D} 35\n');out.write(' '.join(map(str,cache['rows']))+'\n');writepoly(out,expr_poly(cache['ud']))
        for key in ['J','T','U']:
            for row in cache[key]:
                for v in row:writepoly(out,specialize(v))
        for f in cache['F']:
            terms=[]
            for mon,c in f:
                if cache['fixed'] in mon:continue
                i,j=[cache['chart'].index(k) for k in mon]+[44]*(2-len(mon));terms.append((i,j,expr_poly(c)))
            out.write(str(len(terms))+'\n')
            for i,j,c in terms:out.write(f'{i} {j}\n');writepoly(out,c)
    print('native input bytes',inp.stat().st_size,'seconds',time.time()-start,flush=True)
    proc=subprocess.run([str(BINARY),str(inp),str(output)],env=env,check=True,timeout=240)
    print('native output bytes',output.stat().st_size,'seconds',time.time()-start,flush=True)
    matrices={}
    active={4,5,7,8}
    with output.open() as f:
        assert list(map(int,next(f).split()))==[135,np,10,D]
        for d in range(2,D+1):
            mons=list(cwr(range(10),d));idx={m:i for i,m in enumerate(mons)}
            columns=[j for j,m in enumerate(mons)if sum(i in active for i in m)<=2]if d==4 else list(range(len(mons)))
            ci={j:k for k,j in enumerate(columns)}
            for key,count in [('W',45),('A',135)]:
                M=[]
                for i in range(count):
                    keep=key=='A' or d==2 or i==44
                    row=[[] for m in columns]if keep else None
                    terms=int(next(f))
                    for k in range(terms):
                        data=next(f).split()
                        if not keep:continue
                        co=data[0];ex=list(map(int,data[1:]));assert len(ex)==nv
                        normal=tuple(j for j,e in enumerate(ex[np:]) for _ in range(e));assert len(normal)==d
                        j=idx[normal]
                        if j in ci:row[ci[j]].append([ex[:np],co])
                    if keep:M.append(row)
                if key=='W' and d>=3:matrices['T'+str(d)+'_surviving'if d==4 else 'T3']=M
                else:matrices[key+str(d)+'_surviving'if d==4 else key+str(d)]=M
            if d==4:matrices['A4_columns']=columns
    for key in ['W2','A2']:
        expected=cache[key]
        for i,row in enumerate(expected):
            for j,v in enumerate(row):
                lhs={tuple(e):s.Rational(c) for e,c in matrices[key][i][j]};rhs=dict(specialize(v));assert lhs==rhs,(key,i,j)
    print('ALL W2/A2 EXACT CROSS-CHECK PASS',time.time()-start,flush=True)
    matrices.update(variables=list(map(str,variables)),ud=str(s.sympify(cache['ud']) if general else s.sympify(cache['ud']).subs(oldvars[0],r0)),td=str(s.sympify(cache['td']) if general else s.sympify(cache['td']).subs(oldvars[0],r0)),seconds=time.time()-start)
    pickle.dump(matrices,open(WD/(label+'-flint-projection-input.pkl'),'wb'))
print('PARAMETRIC FOUR-JET CACHE READY',time.time()-start,flush=True)
