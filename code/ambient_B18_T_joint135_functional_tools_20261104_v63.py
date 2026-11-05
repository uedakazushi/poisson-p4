#!/usr/bin/env python3
"""Finite source algebra and actual proper-map evaluation for selected135."""
import collections,itertools
from fractions import Fraction as F
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_T_source2_cache_composition_20261104_v63 import MON as MON20
MON=MON20+[(7,7,7),(0,7,7),(1,7,7),(4,7,7),(5,7,7)]
INDEX={m:i for i,m in enumerate(MON)}
class Finite:
    def __init__(self,m=0,nus=()):
        self.m=m;self.b=1<<m;self.n=25*self.b;self.Z=(F(0),)*self.n;self.ONE=self.c(1);self.table={};self.nilcap=11
        def bp(a,b):
            z={}
            for i,c in a.items():
                for j,d in b.items():
                    k=INDEX.get(tuple(sorted(MON[i]+MON[j])))
                    if k is not None:z[k]=z.get(k,F(0))+c*d
            return {i:c for i,c in z.items() if c}
        for i in range(self.n):
            ai,ti=divmod(i,self.b)
            for j in range(self.n):
                aj,tj=divmod(j,self.b);v=bp({ai:F(1)},{aj:F(1)});mask=ti^tj;carry=ti&tj
                for k in range(m):
                    if carry&(1<<k):v=bp(v,{a:c for a,c in enumerate(nus[k]) if c})
                if v:self.table[i,j]=[(self.b*a+mask,c) for a,c in v.items()]
    def c(self,x):return (F(x),)+(F(0),)*(self.n-1)
    def add(self,*ps):
        v=[F(0)]*self.n
        for p in ps:
            for i,c in enumerate(p):
                if c:v[i]+=c
        return tuple(v)
    def neg(self,p):return tuple(-c for c in p)
    def scale(self,p,c):return tuple(c*a for a in p)
    def mul2(self,p,q):
        if p==self.Z or q==self.Z:return self.Z
        if p==self.ONE:return q
        if q==self.ONE:return p
        v=[F(0)]*self.n;qq=[(j,c) for j,c in enumerate(q) if c]
        for i,c in enumerate(p):
            if c:
                for j,d in qq:
                    for k,a in self.table.get((i,j),()):v[k]+=c*d*a
        return tuple(v)
    def mul(self,*ps):
        v=self.ONE
        for p in ps:v=self.mul2(v,p)
        return v
    def power(self,p,n):
        if n<0:return self.power(self.inv(p),-n)
        if n>=self.nilcap and not p[0]:return self.Z
        v=self.ONE
        while n:
            if n&1:v=self.mul2(v,p)
            n//=2
            if n:p=self.mul2(p,p)
        return v
    def inv(self,p):
        assert p[0];q=list(p);q[0]=F(0);q=self.scale(tuple(q),-1/p[0]);v=self.ONE;t=self.ONE
        for _ in range(self.nilcap-1):t=self.mul(t,q);v=self.add(v,t)
        v=self.scale(v,1/p[0]);assert self.mul(p,v)==self.ONE;return v
    def dual(self,p):
        v=[F(0)]*self.n
        for i,c in enumerate(p):v[self.b*i]=F(c)
        return tuple(v)
    def tau(self,k):
        v=[F(0)]*self.n;v[1<<k]=F(1);return tuple(v)
    def install(self):
        for k in ['c','add','neg','scale','mul','power','inv']:setattr(H,k,getattr(self,k))
        H.Z=self.Z;H.ONE=self.ONE;H.TA=self.tau(0) if self.m else self.Z;H.TB=self.tau(1) if self.m>1 else self.Z
def source_evaluate(nodes,roots,bind,a):
    ns=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in ns:continue
        ns.add(i)
        if i in bind:continue
        n=nodes[i]
        if n[0] in ['a','m','n','i']:todo.extend(n[1:])
        else:assert n[0]=='c',('UNCACHED_SOURCE_LEAF',i,n)
    uses=collections.Counter()
    for i in ns:
        if i not in bind and nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    uses.update(roots);v={};peak=0
    for i in sorted(ns):
        n=nodes[i];o=n[0]
        if i in bind:p=bind[i]
        elif o=='c':p=a.c(F(n[1],n[2]))
        elif o=='a':p=a.add(v[n[1]],v[n[2]])
        elif o=='n':p=a.neg(v[n[1]])
        elif o=='m':p=a.mul(v[n[1]],v[n[2]])
        elif o=='i':p=a.inv(v[n[1]])
        else:raise AssertionError(n)
        v[i]=p;peak=max(peak,len(v)*a.n);assert peak+30000<150000
        if i not in bind:
            for j in n[1:] if o in ['a','m','n','i'] else []:
                uses[j]-=1
                if not uses[j]:del v[j]
    return {i:v[i] for i in roots},len(ns),peak
def source_first_AD(nodes,roots,origin):
    ns=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in ns:continue
        ns.add(i)
        if i in origin:continue
        n=nodes[i]
        if n[0] in ['a','m','n','i']:todo.extend(n[1:])
        else:assert n[0]=='c'
    uses=collections.Counter()
    for i in ns:
        if i not in origin and nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    uses.update(roots);v={};g={};peak=0;live=0
    def add(a,b,c=F(1)):
        a=dict(a)
        for i,x in b.items():
            z=a.get(i,F(0))+c*x
            if z:a[i]=z
            else:a.pop(i,None)
        return a
    for i in sorted(ns):
        n=nodes[i];o=n[0]
        if i in origin:p=origin[i];q={i:F(1)}
        elif o=='c':p=F(n[1],n[2]);q={}
        elif o=='a':p=v[n[1]]+v[n[2]];q=add(g[n[1]],g[n[2]])
        elif o=='n':p=-v[n[1]];q={j:-c for j,c in g[n[1]].items()}
        elif o=='m':p=v[n[1]]*v[n[2]];q=add({j:c*v[n[2]] for j,c in g[n[1]].items() if c*v[n[2]]},g[n[2]],v[n[1]])
        elif o=='i':p=1/v[n[1]];q={j:-c*p*p for j,c in g[n[1]].items()}
        else:raise AssertionError(n)
        v[i]=p;g[i]=q;live+=1+len(q);peak=max(peak,live);assert peak+30000<150000
        if i not in origin:
            for j in n[1:] if o in ['a','m','n','i'] else []:
                uses[j]-=1
                if not uses[j]:live-=1+len(g[j]);del v[j];del g[j]
    return {i:(v[i],g[i]) for i in roots},peak
class Proper:
    def __init__(self,S,I,P,K,prep,GL):self.S=S;self.I=I;self.P=P;self.K=K;self.prep=prep;self.GL=GL
    def evaluate(self,a,phi,pairs,ev,od):
        a.install();S,I,P,K,prep,GL=self.S,self.I,self.P,self.K,self.prep,self.GL
        def frame(family,tauindex):
            f=S['families'][family];p=H.add(a.tau(tauindex),H.scale(a.dual(phi[f['p_quadratic_relation']['s']]),F(1,2)));G=H.ident(5)
            for name,(aa,bb) in f['coefficient_unknown_pairs'].items():
                if name.startswith('gamma') and name!='gamma0':
                    i,j=divmod(GL[int(name[5:])-1],5);G[i][j]=H.add(G[i][j],a.dual(phi[aa]),H.mul(a.dual(phi[bb]),p))
            core=H.ident(5)
            if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
            else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
            return H.mm(G,core)
        G=H.ident(5)
        for (i,j),x,y in zip(I['actual_group15_positions'],ev,od):G[i][j]=H.add(G[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
        HA=frame('CA',0);HB=frame('CB',1);FA=H.frame(H.mm(G,HA));FB=H.frame(H.mm(G,HB))
        PA=[H.add(a.dual(x),H.mul(a.dual(y),a.tau(0))) for x,y in pairs['CA']]
        rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];centers=list(map(F,P['actual_CA_source19_center']))
        def own(tr):
            p=H.physical_push(PA,tr,H.mi(tr));v=H.mm(H.constmat(P['CA_left_inverse19']),[[p[rows[j]]] for j in ar]);return [q[0] for q in v]
        ow=own(G);owc=own(H.mm(FA,G))
        constraints=[FA[1][0],FA[1][3],FA[1][4],FA[0][3],FA[0][4]]+[H.add(owc[j],H.c(-centers[j])) for j in P['CA_fiber10_fixed_coefficient_indices']]
        actual=[H.add(ow[j],H.c(-centers[j])) for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
        kap=[actual[j] for j in I['retained_source9_indices']];core=[kap[j] for j in [0,1,2,6]];du=H.add(a.ONE,H.neg(H.mul(kap[1],kap[1])));pu=H.add(H.c(5),H.mul(kap[1],kap[1]));rr=[]
        inds=[K['source9_variable_order'].index(n) for n in K['selected5_pivot_order']]
        for name in K['selected5_pivot_order']:
            p=K['five_actual_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
        gs=[H.add(kap[j],H.neg(v)) for j,v in zip(inds,rr)]+[H.evalpoly(p,core) for p in prep['generators3']]
        return constraints,gs,kap
