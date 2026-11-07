#!/usr/bin/env python3
"""Independent complete eight-family rank-two nonkernel plane review.

Fraction arithmetic, sparse QQ[tau,eta] Bareiss determinants, exact rational
functions and four quadratic number fields; no producer or CAS imports.
"""
import ast
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from math import isqrt
from pathlib import Path
import resource
import signal
import time

from classification_unimodular_logarithmic_review import RF, ud
import classification_degenerate_quadric_q4_nil31_I_review as q
import classification_plane_rank_zero_rank_two_review as previous
from classification_minimal_quartic_review import LaurentRing, determinant

ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_20261107.json'
OUTPUT='reports/CLASSIFICATION_PLANE_RANK_ZERO_NONKERNEL_REVIEW_20261107.json'
JOB='reports/CLASSIFICATION_PLANE_RANK_ZERO_NONKERNEL_REVIEW_JOB_20261107.json'
PAIRS=list(combinations(range(5),2));NAMES=('a','b','c','z','w');ZERO=(0,)*8

# Bivariate polynomial arithmetic uses dictionaries, with lexicographic exact
# division. Bareiss divisions are checked with remainder zero, never numerical.
def pa(*Ps):
    out={}
    for P in Ps:
        for e,c in P.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}
def ps(P,c):return {e:c*v for e,v in P.items() if c*v}
def pm(P,Q):
    out={}
    for e,c in P.items():
        for f,d in Q.items():
            g=(e[0]+f[0],e[1]+f[1]);out[g]=out.get(g,F(0))+c*d
    return {e:c for e,c in out.items() if c}
def pd(P,Q):
    assert Q
    remainder=dict(P);out={};f=max(Q);c=Q[f]
    while remainder:
        e=max(remainder);assert all(v>=w for v,w in zip(e,f)),('inexact bivariate division',e,f)
        g=(e[0]-f[0],e[1]-f[1]);d=remainder[e]/c;out[g]=out.get(g,F(0))+d
        remainder=pa(remainder,ps(pm({g:F(1)},Q),-d))
    assert pm(out,Q)==P
    return out
PONE={(0,0):F(1)}

class R2:
    def __init__(self,value=0,den=None):
        if isinstance(value,R2):self.p,self.q=value.p,value.q;return
        if isinstance(value,RF):
            self.p={(i,0):c for i,c in enumerate(value.p) if c};self.q={(i,0):c for i,c in enumerate(value.q) if c};return
        self.p=dict(value) if isinstance(value,dict) else ({(0,0):F(value)} if value else {})
        self.q=dict(den) if den is not None else dict(PONE)
        assert self.q
        if not self.p:self.q=dict(PONE)
        elif len(self.q)==1 and (0,0) in self.q:
            self.p=ps(self.p,1/self.q[(0,0)]);self.q=dict(PONE)
    def __add__(self,v):
        v=R2(v)
        if self.q==v.q:return R2(pa(self.p,v.p),self.q)
        return R2(pa(pm(self.p,v.q),pm(v.p,self.q)),pm(self.q,v.q))
    __radd__=__add__
    def __neg__(self):return R2(ps(self.p,-1),self.q)
    def __sub__(self,v):return self+-R2(v)
    def __rsub__(self,v):return R2(v)+-self
    def __mul__(self,v):
        v=R2(v);return R2(pm(self.p,v.p),pm(self.q,v.q))
    __rmul__=__mul__
    def __truediv__(self,v):
        v=R2(v);assert v.p;return R2(pm(self.p,v.q),pm(self.q,v.p))
    def __rtruediv__(self,v):return R2(v)/self
    def __bool__(self):return bool(self.p)
    def __eq__(self,v):
        v=R2(v);return pm(self.p,v.q)==pm(v.p,self.q)
    def __repr__(self):return repr((self.p,self.q))


def quadratic_field(d):
    assert d
    assert d<0 or isqrt(d.numerator)**2!=d.numerator or isqrt(d.denominator)**2!=d.denominator
    class Q2:
        def __init__(self,a=0,b=0):
            if isinstance(a,Q2):self.a,self.b=a.a,a.b;return
            self.a,self.b=F(a),F(b)
        def __add__(self,v):v=Q2(v);return Q2(self.a+v.a,self.b+v.b)
        __radd__=__add__
        def __neg__(self):return Q2(-self.a,-self.b)
        def __sub__(self,v):return self+-Q2(v)
        def __rsub__(self,v):return Q2(v)+-self
        def __mul__(self,v):v=Q2(v);return Q2(self.a*v.a+d*self.b*v.b,self.a*v.b+self.b*v.a)
        __rmul__=__mul__
        def __truediv__(self,v):
            v=Q2(v);norm=v.a*v.a-d*v.b*v.b;assert norm
            return self*Q2(v.a/norm,-v.b/norm)
        def __rtruediv__(self,v):return Q2(v)/self
        def __bool__(self):return bool(self.a or self.b)
        def __eq__(self,v):v=Q2(v);return self.a==v.a and self.b==v.b
        def __repr__(self):return str(self.a)+' + '+str(self.b)+'*alpha'
    return Q2


def read_scalar(text,env,field,root_polynomial=None):
    def read(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return field(n.value)
        if isinstance(n,ast.Name):return env[n.id]
        if isinstance(n,ast.UnaryOp):return -read(n.operand) if isinstance(n.op,ast.USub) else read(n.operand)
        if isinstance(n,ast.Call):
            assert isinstance(n.func,ast.Name) and n.func.id=='CRootOf' and len(n.args)==2
            assert root_polynomial is not None and n.args[1].value==0
            assert ast.dump(n.args[0])==ast.dump(ast.parse(root_polynomial,mode='eval').body)
            return env['alpha']
        if isinstance(n,ast.BinOp):
            a=read(n.left)
            if isinstance(n.op,ast.Pow):
                assert isinstance(n.right,ast.Constant) and n.right.value>=0
                out=field(1)
                for _ in range(n.right.value):out=out*a
                return out
            b=read(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
        raise ValueError(ast.dump(n))
    return read(ast.parse(str(text),mode='eval').body)


def parse(text,env,field):
    def atom(v=1,k=None):
        ex=[0]*8
        if k is not None:ex[k]=1
        return {tuple(ex):field(v)} if v else {}
    def read(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return atom(n.value)
        if isinstance(n,ast.Name):return atom(env[n.id]) if n.id in env else atom(1,NAMES.index(n.id))
        if isinstance(n,ast.UnaryOp):return q.scale(read(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
        if isinstance(n,ast.BinOp):
            a=read(n.left)
            if isinstance(n.op,ast.Pow):
                out=atom()
                for _ in range(n.right.value):out=q.mul(out,a)
                return out
            b=read(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            if isinstance(n.op,ast.Div):
                assert set(b)=={ZERO};return q.scale(a,field(1)/b[ZERO])
        raise ValueError(ast.dump(n))
    return read(ast.parse(str(text),mode='eval').body)


def det_generic(M):
    B=[list(r) for r in M];out=F(1)
    for j in range(len(B)):
        k=next((k for k in range(j,len(B)) if B[k][j]),None)
        if k is None:return out*0
        if k!=j:B[j],B[k]=B[k],B[j];out=-out
        a=B[j][j];out=out*a
        for k in range(j+1,len(B)):
            if not B[k][j]:continue
            c=B[k][j]/a
            for l in range(j+1,len(B)):B[k][l]=B[k][l]-c*B[j][l]
            B[k][j]=a*0
    return out


def det_bivariate(M):
    assert all(v.q==PONE for row in M for v in row)
    B=[[dict(v.p) for v in row] for row in M];previous=dict(PONE);sgn=1;n=len(B)
    if n==0:return R2(1)
    for j in range(n-1):
        k=next((k for k in range(j,n) if B[k][j]),None)
        if k is None:return R2()
        if k!=j:B[j],B[k]=B[k],B[j];sgn=-sgn
        pivot=B[j][j]
        for k in range(j+1,n):
            for l in range(j+1,n):B[k][l]=pd(pa(pm(B[k][l],pivot),ps(pm(B[k][j],B[j][l]),-1)),previous)
            B[k][j]={}
        previous=pivot
    return R2(ps(B[-1][-1],sgn))


def matrix(text,env,field,root_poly=None):return [[read_scalar(v,env,field,root_poly) for v in row] for row in text]
def encode(value):
    if isinstance(value,R2):return dict(numerator=[[list(e),str(c)] for e,c in sorted(value.p.items())],denominator=[[list(e),str(c)] for e,c in sorted(value.q.items())])
    if isinstance(value,RF):return dict(numerator=list(map(str,value.p)),denominator=list(map(str,value.q)))
    if hasattr(value,'a'):return [str(value.a),str(value.b)]
    return str(value)


def certificate(M,C,env,field,root_poly=None):
    n=len(M[0]);r=C['rank'];k=n-r
    assert C['nullity']==k
    rr,cc=C['minor_rows'],C['minor_columns']
    assert len(rr)==len(cc)==r and len(set(rr))==len(set(cc))==r
    sub=[[M[i][j] for j in cc] for i in rr]
    value=det_bivariate(sub) if field is R2 else det_generic(sub)
    assert value and value==read_scalar(C['minor_determinant'],env,field,root_poly)
    K=matrix(C['complete_kernel'],env,field,root_poly)
    assert len(K)==n and all(len(row)==k for row in K)
    for row in M:
        for col in zip(*K):
            out=field(0)
            for a,b in zip(row,col):
                if a and b:out=out+a*b
            assert not out
    # RREF kernel contains the full identity on its free-variable rows.
    identity_rows=[]
    for j in range(k):
        target=[field(int(l==j)) for l in range(k)]
        index=next((i for i,row in enumerate(K) if row==target),None)
        assert index is not None;identity_rows.append(index)
    assert len(set(identity_rows))==k
    return dict(rank=r,nullity=k,all_rows_kernel=True,independent_free_rows=identity_rows,
                minor_rows=rr,minor_columns=cc,nonzero_minor=encode(value)),value


def Padd(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{}) for P in Ps)) for ij in PAIRS})
def wedge(X,Y):return q.antisymmetric({ij:q.add(q.mul(X[ij[0]],Y[ij[1]]),q.scale(q.mul(X[ij[1]],Y[ij[0]]),-1)) for ij in PAIRS})


def model(name,env,field):
    p=lambda text:parse(text,env,field)
    En=[{},{},{},p('z'),p('w')]
    if name.startswith('ND'):
        V=[{},{},{},p('c*w'),p('b*z')]
        U=[{},p('1'),p('1') if name.startswith('NDSS') else {},{},{}]
        K=p('-(z**2-w**2)/2' if name.startswith('NDSS') else '-z**2/2')
        if name=='NDSS':T=[p('2*tau*a'),p('-tau*b'),p('-tau*c'),{},{}]
        elif name=='NDSSN':T=[p('b-c'),{},{},{},{}]
        elif name=='NDNSS':T=[p('2*a'),p('-b'),p('-c'),{},{}]
        else:assert name=='NDN';T=[p('c'),{},{},{},{}]
    else:
        V=[{},{},{},p('b*z+c*w'),p('-b*w')]
        U=[{},p('1') if name.startswith('BSS') else {},{} if name.startswith('BSS') else p('1'),{},{}]
        if name=='BSS':T=[p('(2*tau-2*eta)*a'),p('eta*b'),p('(eta-2*tau)*c'),{},{}]
        elif name=='BSSR':T=[p('-2*tau*a/3+c'),p('4*tau*b/3'),p('-2*tau*c/3'),{},{}]
        elif name=='BNS':T=[p('-2*a'),p('b'),p('c+2*tau*b'),{},{}]
        else:assert name=='BNR';T=[p('b'),{},p('2*tau*b'),{},{}]
        K=p('(1-2*tau)*z*w' if name.startswith('BSS') else '(1/2-tau)*w**2')
    return Padd(wedge(T,En),wedge(U,V),q.antisymmetric({(3,4):K}))


def shear(P,basis,field):
    x=[q.atom(1,i) for i in range(5)];columns=[]
    for normal in (3,4):
        for i in range(3):
            Y=[{} for _ in range(5)];Y[i]=q.atom(1,normal)
            orbit=q.antisymmetric({ij:q.add(
                *(q.scale(q.mul(Y[k],q.der(P.get(ij,{}),k)),-1) for k in range(5)),
                *(q.mul(P.get((k,ij[1]),{}),q.der(Y[ij[0]],k)) for k in range(5)),
                *(q.mul(P.get((ij[0],k),{}),q.der(Y[ij[1]],k)) for k in range(5))) for ij in PAIRS})
            first=[orbit.get((k,3),{}).get((0,0,0,2,0,0,0,0),field(0)) for k in range(3)]
            second=[orbit.get((k,4),{}).get((0,0,0,0,2,0,0,0),field(0)) for k in range(3)]
            V=[q.add(q.scale(x[3],first[k]),q.scale(x[4],second[k])) for k in range(3)]+[{},{}]
            canonical=Padd(orbit,wedge(x,V));coeffs=[]
            for B in basis:
                ij=next(ij for ij in PAIRS if B.get(ij));ex=next(iter(B[ij]));assert B[ij][ex]==1
                coeffs.append(canonical.get(ij,{}).get(ex,field(0)))
            assert canonical==q.antisymmetric({ij:q.add(*(q.scale(B.get(ij,{}),v) for B,v in zip(basis,coeffs))) for ij in PAIRS})
            columns.append(coeffs)
    return list(map(list,zip(*columns)))


def divisor(P):
    return q.add(*(q.scale(q.mul(q.atom(1,i),previous.pf4(P,[j for j in range(5) if j!=i])),(-1)**i) for i in range(5)))


def rf_at(value,t):
    value=RF(value)
    def evaluate(poly):
        out=t*0
        for c in reversed(poly):out=out*t+c
        return out
    return evaluate(value.p)/evaluate(value.q)


def r2_at(value,t,s):
    value=R2(value)
    def power(a,k):
        out=a*0+1
        for _ in range(k):out=out*a
        return out
    def evaluate(P):return sum((c*power(t,e[0])*power(s,e[1]) for e,c in P.items()),t*0)
    return evaluate(value.p)/evaluate(value.q)


def convert_matrix(M,fn):return [[fn(v) for v in row] for row in M]


def combined(Ms,C,env,field,affine,excluded=None,root_poly=None):
    certs={};product=field(1)
    for key,M in zip(('D1','D2','actual_shear'),Ms):
        cert,value=certificate(M,C[key],env,field,root_poly);certs[key]=cert;product=product*value
    d=certs['D1']['nullity']+certs['D2']['nullity'];bound=affine+6-1+d
    assert d==C['higher_fiber_upper_bound'] and bound==C['PGL_image_upper_bound']
    is_excluded=excluded is not None and not excluded
    if not is_excluded:assert bound<=23
    return dict(**certs,affine_symbol_dimension=affine,higher_fiber_upper_bound=d,
                projective_image_upper_bound=bound,excluded_higher_multiplicity=is_excluded),product


def factor_univariate(value,analysis):
    assert value.q==(F(1),)
    product=RF(1);factors=[]
    for text,multiplicity in analysis['all_irreducible_factors']:
        f=read_scalar(text,{'tau':previous.RHO},RF)
        assert f.q==(F(1),) and len(f.p) in (2,3) and multiplicity>0
        if len(f.p)==3:
            assert not f.p[1];d=-f.p[0]/f.p[2]
            assert d<0 or isqrt(d.numerator)**2!=d.numerator or isqrt(d.denominator)**2!=d.denominator
        for _ in range(multiplicity):product=product*f
        factors.append((text,f))
    ratio=value/product
    assert len(ratio.p)==1 and ratio.q==(F(1),) and ratio.p[0]
    assert value==read_scalar(analysis['exception_polynomial'],{'tau':previous.RHO},RF)
    return factors


def univariate(Ms,analysis,affine,excluded,position=None):
    env={'tau':previous.RHO}
    result,product=combined(Ms,analysis,env,RF,affine)
    factors=factor_univariate(product,analysis)
    assert len(factors)==len(analysis['all_exceptional_roots'])
    walls=[];covered=[]
    for text,f in factors:
        if len(f.p)==2:
            rho=-f.p[0]/f.p[1]
            matches=[w for w in analysis['all_exceptional_roots'] if 'minimal_polynomial' not in w and F(w['parameter'])==rho]
            assert len(matches)==1;wall=matches[0];field=RF;value=RF(rho);root_poly=None
        else:
            matches=[w for w in analysis['all_exceptional_roots'] if w.get('minimal_polynomial')==text]
            assert len(matches)==1;wall=matches[0];d=-f.p[0]/f.p[2]
            field=quadratic_field(d);value=field(0,1);root_poly=text
            assert wall['all_conjugate_roots_covered']
        specialized=[convert_matrix(M,lambda v:rf_at(v,value)) for M in Ms]
        ex=rf_at(excluded,value)
        cert,_=combined(specialized,wall,{'tau':value,'alpha':value},field,affine-1,ex,root_poly)
        assert cert['excluded_higher_multiplicity']==wall['excluded_multiplicity_higher']
        cert.update(minimal_polynomial=text,parameter=encode(value),all_roots_by_field_embeddings=len(f.p)==3)
        if position is not None:
            point=tuple(rf_at(v,value) for v in position)
            cert['point']=[encode(v) for v in point]
        walls.append(cert);covered.append(text)
    result.update(exceptional_factors=covered,all_exceptional_roots=walls)
    return result,product


def bivariate(Ms,analysis):
    env={'tau':R2({(1,0):F(1)}),'eta':R2({(0,1):F(1)})}
    result,product=combined(Ms,analysis,env,R2,12)
    factor_product=R2(1);factor_keys=[]
    for text,multiplicity in analysis['all_irreducible_factors']:
        f=read_scalar(text,env,R2);assert f.q==PONE and f.p and max(sum(e) for e in f.p)==1
        key=tuple(f.p.get(e,F(0)) for e in [(1,0),(0,1),(0,0)])
        scale=next(v for v in key if v);key=tuple(v/scale for v in key)
        assert key not in factor_keys;factor_keys.append(key)
        for _ in range(multiplicity):factor_product=factor_product*f
    assert product==read_scalar(analysis['exception_polynomial'],env,R2)
    ratio=pd(product.p,factor_product.p)
    assert len(ratio)==1 and (0,0) in ratio and ratio[(0,0)]
    lines=analysis['all_exceptional_lines'];assert len(lines)==len(factor_keys)==17
    assert {line['equation'] for line in lines}=={text for text,m in analysis['all_irreducible_factors']}
    reviewed=[];substitutions=[]
    for line in lines:
        sub={key:read_scalar(value,{'tau':previous.RHO},RF) for key,value in line['parameterization'].items()}
        assert set(sub)=={'tau','eta'}
        t,s=sub['tau'],sub['eta']
        assert all(v.q==(F(1),) and len(v.p)<=2 for v in (t,s))
        assert any(len(v.p)==2 and v.p[1] for v in (t,s))
        equation=read_scalar(line['equation'],env,R2)
        assert not r2_at(equation,t,s)
        ex=4*t-3*s;excluded=not ex
        assert excluded==line['excluded_multiplicity_higher']
        if excluded:
            reviewed.append(dict(equation=line['equation'],excluded_higher_multiplicity=True,
                                 complementary_m2_nonzero_row='BSSR whole family',all_strata_checked=True))
            substitutions.append((line,sub,None));continue
        specialized=[convert_matrix(M,lambda v:r2_at(v,t,s)) for M in Ms]
        cert,minor=univariate(specialized,line['analysis'],11,ex,position=(t,s))
        cert.update(equation=line['equation'],parameterization=line['parameterization'])
        reviewed.append(cert);substitutions.append((line,sub,minor))
    # Explicit check that every pairwise finite line intersection belongs to a
    # certified generic line patch or one of that line's retained root fields.
    intersections=[]
    for i,(line,sub,minor) in enumerate(substitutions):
        if minor is None:continue
        t,s=sub['tau'],sub['eta']
        for j,(other,_,_) in enumerate(substitutions):
            if j<=i:continue
            equation=r2_at(read_scalar(other['equation'],env,R2),t,s)
            assert equation.q==(F(1),) and len(equation.p)<=2
            if len(equation.p)<2:assert equation;continue
            rho=-equation.p[0]/equation.p[1]
            ex=rf_at(4*t-3*s,RF(rho))
            kind='excluded Delta0' if not ex else 'line maximal minor nonzero' if rf_at(minor,RF(rho)) else 'explicit retained root'
            if kind=='explicit retained root':
                assert any('minimal_polynomial' not in w and F(w['parameter'])==rho for w in line['analysis']['all_exceptional_roots'])
            intersections.append(dict(lines=[line['equation'],other['equation']],line_parameter=str(rho),coverage=kind))
    result.update(all_seventeen_linear_factors=[line['equation'] for line in lines],all_exceptional_lines=reviewed,
                  explicit_pairwise_intersections=intersections,excluded_line_count=1)
    return result


def normalization_identities():
    R=LaurentRing(('s','ell','tau','sigma','d'),NAMES,invert=('s','ell','sigma','d'))
    p=R.parse; x=list(map(R.variable,NAMES))
    def wed(X,Y):return {ij:R.add(R.mul(X[ij[0]],Y[ij[1]]),R.scale(R.mul(X[ij[1]],Y[ij[0]]),-1)) for ij in PAIRS}
    def make(kind,T):
        En=[{},{},{},p('z'),p('w')]
        if kind in ('NDSS','NDN'):
            U=[{},p('1'),p('1') if kind=='NDSS' else {},{},{}];V=[{},{},{},p('c*w'),p('b*z')]
            K=p('-(z**2-w**2)/2' if kind=='NDSS' else '-z**2/2')
        else:
            U=[{},p('1') if kind=='BSS' else {},{} if kind=='BSS' else p('1'),{},{}]
            V=[{},{},{},p('b*z+c*w'),p('-b*w')]
            K=p('(1-2*tau)*z*w' if kind=='BSS' else '(1/2-tau)*w**2')
        A=wed(list(map(p,T))+[{},{}],En);B=wed(U,V)
        return {ij:R.add(A[ij],B[ij],K if ij==(3,4) else {}) for ij in PAIRS}
    def subst(poly,images):
        out={}
        for ex,c in poly.items():
            value={ex[:R.n]+(0,)*5:c}
            for k,n in enumerate(ex[R.n:]):
                for _ in range(n):value=R.mul(value,images[k])
            out=R.add(out,value)
        return out
    def get(P,i,j):return P.get((i,j),{}) if i<j else R.scale(P.get((j,i),{}),-1) if i>j else {}
    specs=[
      ('NDSS remove row','NDSS',['2*s*a+ell*(b-c)','-s*b','-s*c'],['a+ell*(b-c)/(3*s)','b','c','z','w'],['a-ell*(b-c)/(3*s)','b','c','z','w'],'1',['2*s*a','-s*b','-s*c']),
      ('ND nil remove row','NDN',['2*s*a+ell*c','-s*b','-s*c'],['a+ell*c/(3*s)','b','c','z','w'],['a-ell*c/(3*s)','b','c','z','w'],'1',['2*s*a','-s*b','-s*c']),
      ('ND nil scale nonzero s','NDN',['2*s*a','-s*b','-s*c'],['a','b','sigma**2*c','sigma*z','w'],['a','b','c/sigma**2','z/sigma','w'],'sigma',['2*s*sigma*a','-s*sigma*b','-s*sigma*c']),
      ('BSS remove row Delta=d','BSS',['(2*tau-2*s)*a+ell*c','s*b','(s-2*tau)*c'],['a+ell*c/d','b','c','z','w'],['a-ell*c/d','b','c','z','w'],'1',['(2*tau-2*s)*a','s*b','(s-2*tau)*c']),
      ('BNS remove row','BN',['-2*s*a+ell*b','s*b','s*c+2*tau*b'],['a-ell*b/(3*s)','b','c','z','w'],['a+ell*b/(3*s)','b','c','z','w'],'1',['-2*s*a','s*b','s*c+2*tau*b']),
      ('BNS scale nonzero s','BN',['-2*s*a','s*b','s*c+2*tau*b'],['a','sigma*b','c','z/sigma','w'],['a','b/sigma','c','sigma*z','w'],'sigma',['-2*s*sigma*a','s*sigma*b','s*sigma*c+2*tau*b'])]
    results=[]
    for label,kind,T,G,inv,scalar,Tout in specs:
        P=make(kind,T);G=list(map(p,G));inv=list(map(p,inv));expected=make(kind,Tout)
        assert [subst(v,inv) for v in G]==x and [subst(v,G) for v in inv]==x
        L=[[R.derivative(v,k) for k in range(5)] for v in G]
        for ij in PAIRS:
            pushed=R.mul(p(scalar),subst(R.add(*(R.mul(L[ij[0]][k],get(P,k,l),L[ij[1]][l]) for k in range(5) for l in range(5))),inv))
            difference=R.add(pushed,R.scale(expected[ij],-1))
            if 'Delta=d' in label:
                # Evaluate d=4tau-3s by polynomial multiplication after clearing
                # the finitely many Laurent d powers. No inverse of a sum is used.
                minimum=min((ex[4] for ex in difference),default=0)
                out={}
                for ex,c in difference.items():
                    f=list(ex);f[4]=0;term={tuple(f):c}
                    for _ in range(ex[4]-minimum):term=R.mul(term,p('4*tau-3*s'))
                    out=R.add(out,term)
                assert not out
            else:assert not difference,(label,ij)
        results.append(dict(normalization=label,both_inverse_compositions=True,all_ten_coefficients=True,
                            map=[R.encode(v) for v in G],inverse=[R.encode(v) for v in inv],
                            tensor_scalar=scalar,determinant=R.encode(determinant(R,L))))
    # Every zero-diagonal/nonzero-row branch is normalized by scaling a alone.
    for kind,T,Tout in [('NDSS',['ell*(b-c)','0','0'],['b-c','0','0']),('NDN',['ell*c','0','0'],['c','0','0']),
                        ('BSS',['-2*tau*a/3+ell*c','4*tau*b/3','-2*tau*c/3'],['-2*tau*a/3+c','4*tau*b/3','-2*tau*c/3']),
                        ('BN',['ell*b','0','2*tau*b'],['b','0','2*tau*b'])]:
        P=make(kind,T);expected=make(kind,Tout);G=list(map(p,['a/ell','b','c','z','w']));inv=list(map(p,['ell*a','b','c','z','w']))
        L=[[R.derivative(v,k) for k in range(5)] for v in G]
        assert [subst(v,inv) for v in G]==x and [subst(v,G) for v in inv]==x
        for ij in PAIRS:
            pushed=subst(R.add(*(R.mul(L[ij[0]][k],get(P,k,l),L[ij[1]][l]) for k in range(5) for l in range(5))),inv)
            assert pushed==expected[ij]
        results.append(dict(normalization=kind+' nonzero row to1',all_ten_coefficients=True,both_inverse_compositions=True))
    return results


def main():
    started=time.monotonic();signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    output=ROOT/OUTPUT
    if output.exists():raise FileExistsError(output)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    packet=json.loads((ROOT/INPUT).read_text());assert packet['all_passed']
    assert hashlib.sha256(packet['producer_bytes'].encode()).hexdigest()==packet['producer_sha256']
    base1,base2,source_review=previous.source_bases(packet)
    normalizations=normalization_identities();results=[]
    dims={'NDSS':12,'NDSSN':11,'NDNSS':11,'NDN':10,'BSS':12,'BSSR':11,'BNS':11,'BNR':10}
    for c in packet['cases']:
        name=c['type'];field=R2 if name=='BSS' else RF
        env={'tau':R2({(1,0):F(1)}),'eta':R2({(0,1):F(1)})} if field is R2 else {'tau':previous.RHO}
        P=model(name,env,field)
        assert all(P.get(ij,{})==parse(c['normal_form_P0'][str(ij)],env,field) for ij in PAIRS)
        assert not any(q.jacobi(P).values()) and divisor(P)==parse(c['actual_D'],env,field)
        bases=[[q.antisymmetric({ij:{ex:field(v) for ex,v in B.get(ij,{}).items()} for ij in PAIRS}) for B in base] for base in (base1,base2)]
        D1=previous.operator(P,bases[0]);D2=previous.operator(P,bases[1]);S=shear(P,bases[0],field)
        Ms=[convert_matrix(M,field) for M in (D1,D2,S)]
        for M,key in zip(Ms,['complete_D1_350_by30','complete_D2_350_by9','complete_actual_shear_30_by6']):
            assert M==matrix(c[key],env,field)
        assert all(not sum((a*b for a,b in zip(row,col) if a and b),field(0)) for row in Ms[0] for col in zip(*Ms[2]))
        assert all(v.q==(PONE if field is R2 else (F(1),)) for M in Ms for row in M for v in row)
        if name=='BSS':analysis=bivariate(Ms,c['analysis'])
        elif c['parameters']:
            excluded=read_scalar(c['excluded_polynomial'],env,field)
            analysis,_=univariate(Ms,c['analysis'],dims[name],excluded)
        else:analysis,_=combined(Ms,c['analysis'],env,field,dims[name])
        analysis.update(type=name,full_350_by_30_and_9_reconstructed=True,actual_shear_30_by_6_reconstructed=True,
                        full_symbol_Jacobi_and_divisor=True)
        results.append(analysis)
        print(name,'PASS generic image',analysis['projective_image_upper_bound'],flush=True)
    old=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json').read_text())
    funcs=lambda text:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(text).body if isinstance(n,ast.FunctionDef)}
    old_funcs=funcs(old['producer_bytes']);current=funcs(packet['producer_bytes']);assert all(current.get(k)==v for k,v in old_funcs.items())
    result=dict(task_id='G',all_passed=True,scope='Entire mD2 rankzero-plane k1 rankM2/M(p)!=0 eight-family stratum has projective image dimension at most23; no E/fullG claim.',
        complements=['rankM<=1','k0','mD>=3','M(p)0 andrankM3 covered by separate independent reviews'],
        complete_source_review=source_review,universal_actual_normalizations=normalizations,cases=results,
        no_producer_or_CAS_import=True,old_AST_unchanged=True,inputs=job['exact_inputs_and_hashes'],
        producer_sha256=packet['producer_sha256'],producer_snapshot=packet['producer_bytes'],
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-started,
        maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUTPUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)

if __name__=='__main__':main()
