#!/usr/bin/env python3
"""Exact truncated R0 jets with linear/quadratic honest-second parameters."""
from fractions import Fraction as F
from ambient_B18_T_joint135_functional_tools_20261104_v63 import MON as MON25
MON=MON25+[(7,7,7,7)];INDEX={m:i for i,m in enumerate(MON)}
def put(v,k,a):
    z=v.get(k,F(0))+a
    if z:v[k]=z
    else:v.pop(k,None)
def pa(a,b,c=F(1)):
    z=dict(a)
    for k,v in b.items():put(z,k,c*v)
    return z
def pm(a,b):
    z={}
    for i,c in a.items():
        for j,d in b.items():
            mon=tuple(sorted(i+j));assert len(mon)<=2,('HONEST_PARAMETER_DEGREE',mon)
            put(z,mon,c*d)
    return z
def ja(a,b,c=F(1)):
    z={i:dict(p) for i,p in a.items()}
    for i,p in b.items():
        v=pa(z.get(i,{}),p,c)
        if v:z[i]=v
        else:z.pop(i,None)
    return z
def js(a,c):return {i:{m:v*c for m,v in p.items() if v*c} for i,p in a.items() if any(v*c for v in p.values())}
def jm(a,b):
    z={}
    for i,p in a.items():
        for j,q in b.items():
            k=INDEX.get(tuple(sorted(MON[i]+MON[j])))
            if k is not None:
                v=pa(z.get(k,{}),pm(p,q))
                if v:z[k]=v
                else:z.pop(k,None)
    return z
def ji(a):
    assert set(a.get(0,{}))=={()},('NONUNIT_SYMBOLIC_INVERSE',a.get(0));c=a[0][()];assert c
    h={i:p for i,p in a.items() if i};h=js(h,-1/c);v={0:{():F(1)}};t=v
    for _ in range(4):t=jm(t,h);v=ja(v,t)
    v=js(v,1/c);assert jm(a,v)=={0:{():F(1)}};return v
def scalar(a):return {0:{():F(a)}} if a else {}
def known(p):return {int(i):{():F(v)} for i,v in p if F(v)}
def count(a):return sum(map(len,a.values()))
def encode_poly(p):return [[list(k),str(v)] for k,v in sorted(p.items()) if v]
def encode_jet(v):return [[i,encode_poly(p)] for i,p in sorted(v.items()) if p]
