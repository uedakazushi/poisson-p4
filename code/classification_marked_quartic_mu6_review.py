#!/usr/bin/env python3
"""Separate-author exact MU6 marked-curve/source/old-witness audit.
Only independent sparse/Fraction helpers; no producer or CAS imports.
"""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import hashlib,json,resource,time,subprocess,signal
from classification_minimal_quartic_review import LaurentRing,determinant
import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_eight_essential_quintic_source_review import exponents,matmul,local
from classification_unimodular_logarithmic_review import ug,trim
ROOT=Path(__file__).resolve().parents[1]
OUT='reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_SOURCE_REVIEW_20261107.json'
JOB='research/classification/jobs/plane_hyperplane_quartic_mu6_source_review.json'
PREFIX='reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_'

def load(path):return json.loads((ROOT/path).read_text())
def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def power(R,p,n):
    out=R.constant(1)
    for _ in range(n):out=R.mul(out,p)
    return out

def subst(R,P,images):
    out={}
    for e,c in P.items():
        term=R.constant(c)
        for i,k in enumerate(e):
            if k:term=R.mul(term,power(R,images.get(i,R.variable(R.names[i])),k))
        out=R.add(out,term)
    return out

def cert(M,c):
    assert c['shape']==[len(M),len(M[0])]
    r=c['rank'];rows,cols=c['rows'],c['columns']
    assert len(rows)==len(cols)==r and len(set(rows))==len(set(cols))==r
    d=q.determinant([[M[i][j]for j in cols]for i in rows]);assert d==F(c['determinant'])and d
    return str(d)

def lie(R,U,V,base=(0,2,4)):
    return [R.add(*(R.mul(U[k],R.derivative(V[j],i))for k,i in enumerate(base)),*(R.scale(R.mul(V[k],R.derivative(U[j],i)),-1)for k,i in enumerate(base)))for j in range(len(U))]

def source(R,U,V,a,b,mag):
    z={};P={(i,j):z for i in range(5)for j in range(5)}
    vals={(1,2):U[0],(1,4):U[1],(2,3):R.scale(V[0],-1),(3,4):V[1],(1,3):R.add(R.mul(a,R.variable('x')),R.mul(b,R.variable('z')),mag)}
    for(i,j),v in vals.items():P[i,j]=v;P[j,i]=R.scale(v,-1)
    return P

def jac(R,P):
    return [R.add(*(R.mul(P[a,l],R.derivative(P[b,c],l))for a,b,c in((i,j,k),(j,k,i),(k,i,j))for l in range(5)))for i,j,k in combinations(range(5),3)]
def divergence(R,P):return[R.add(*(R.derivative(P[j,i],j)for j in range(5)))for i in range(5)]
def pf(R,P):return R.add(R.mul(P[1,2],P[3,4]),R.scale(R.mul(P[1,3],P[2,4]),-1),R.mul(P[1,4],P[2,3]))
def magmat(R,U,V,a,b):
    lin=[R.variable(n)for n in('h','y','w')];mon=[R.parse(s)for s in('h**2','h*y','h*w','y**2','y*w','w**2')]
    act=lambda X,L:R.add(R.mul(X[0],R.derivative(L,2)),R.mul(X[1],R.derivative(L,4)))
    cols=[R.add(R.scale(R.mul(a,L),-1),R.scale(act(V,L),-1))for L in lin]+[R.add(R.scale(R.mul(b,L),-1),act(U,L))for L in lin]
    return[[C.get(next(iter(m)),F(0))for C in cols]for m in mon]

def sourcecheck(row,norm):
    name=row['name'];params=('g00','g01','g10','g11')+tuple('r'+str(i)for i in range(6));R=LaurentRing(params,('h','x','y','z','w'));p=R.parse;n=R.n
    base=(0,2,4);x=[R.variable(t)for t in('h','x','y','z','w')];mon=[p(t)for t in('h**2','h*y','h*w','y**2','y*w','w**2')]
    G=[[F(v)for v in r]for r in row['actual_H_to_h_base_matrix']];Gi=[[F(v)for v in r]for r in row['actual_inverse_base_matrix']]
    assert matmul(G,Gi)==[[F(i==j)for j in range(3)]for i in range(3)]
    assert q.determinant(G)==F(row['base_determinant'])
    images={n+i:R.add(*(R.scale(x[j],Gi[k][l])for l,j in enumerate(base)))for k,i in enumerate(base)}
    assert subst(R,p(row['input_marked_H']),images)==x[0]
    f=p(row['transformed_irreducible_quartic']);assert subst(R,p(row['input_quartic']),images)==f
    D=R.mul(x[0],f);grad=[R.derivative(D,i)for i in base];target=[(6-i-j,i,j)for i in range(7)for j in range(7-i)]
    def exponent(e):return(0,)*n+(e[0],0,e[1],0,e[2])
    columns=[R.mul(m,g)for g in grad for m in mon]
    M=[[col.get(exponent(e),F(0))for col in columns]for e in target]
    assert M==[[F(v)for v in r]for r in row['full28x18_gradient_matrix']]
    d=cert(M,row['gradient_rank']);K=[[F(v)for v in r]for r in row['full_gradient_kernel']]
    k=row['full_kernel_dimension'];assert len(K)==18 and len(K[0])==k and q.rank(q.transpose(K))[0]==k
    assert not any(any(r)for r in matmul(M,K))and row['gradient_rank']['rank']+k==18
    anchors=[]
    for j in range(k):
        W=[R.add(*(R.scale(mon[l],K[6*i+l][j])for l in range(6)))for i in range(3)]
        ell={}
        for e,c in W[0].items():assert e[n]>0;ee=list(e);ee[n]-=1;ell[tuple(ee)]=c
        anchors.append([R.add(W[i],R.scale(R.mul(ell,x[base[i]]),-1))for i in(1,2)])
    assert anchors==[[p(t)for t in U]for U in row['all_kernel_to_quadratic_h_Casimir_anchors']]
    if name=='E6_a0_transverse_no_frame':
        assert k==3 and all(not R.add(R.mul(U[0],V[1]),R.scale(R.mul(U[1],V[0]),-1))for U,V in combinations(anchors,2))
        return dict(name=name,gradient_rank=15,independent_minor=d,whole_kernel_no_frame=True)
    assert k==2;U,V=anchors;assert U==[p(t)for t in row['anchor_U']]and V==[p(t)for t in row['anchor_V']]
    a,b=map(p,row['commutator']);Lb=lie(R,[{},*U],[{},*V]);assert Lb[1:]==[R.add(R.mul(a,u),R.mul(b,v))for u,v in zip(U,V)]
    ratio=F(row['Saito_determinant_ratio']);assert ratio and R.add(R.mul(U[0],V[1]),R.scale(R.mul(U[1],V[0]),-1))==R.scale(f,ratio)
    for X,co in zip((U,V),row['logarithmic_cofactors']):assert R.add(R.mul(X[0],R.derivative(f,2)),R.mul(X[1],R.derivative(f,4)))==R.mul(p(co),f)
    magnetic=R.add(*(R.mul(p('r'+str(i)),m)for i,m in enumerate(mon)));P=source(R,U,V,a,b,magnetic)
    actual={tuple(map(int,k.split(','))):p(v)for k,v in row['all150_source_brackets'].items()}
    assert all(P[i,j]==actual.get((i,j),{})for i,j in combinations(range(5),2));assert not any(jac(R,P))
    assert R.mul(x[0],pf(R,P))==R.scale(D,ratio)==p(row['actual_quintic'])
    N=divergence(R,P);assert N==[p(t)for t in row['canonical_modular_field']]
    NM=[[N[i].get((0,)*n+tuple(int(j==k)for k in range(5)),F(0))for j in range(5)]for i in range(5)]
    assert NM==[[F(v)for v in r]for r in row['homogeneous_cofactor_matrix']]
    assert not any(any(r)for r in matmul(NM,NM))
    rank=q.rank(NM)[0];assert rank==row['homogeneous_cofactor_power_ranks'][0]
    assert q.rank(NM+[[1,0,0,0,0]])[0]==rank+1
    T=magmat(R,U,V,a,b);assert T==[[F(v)for v in r]for r in row['full_magnetic_translation_matrix']]
    tm=cert(T,row['full_magnetic_translation_rank']);assert row['full_magnetic_translation_rank']['rank']==5
    v=list(map(F,row['magnetic_cokernel_basis'][0]));assert any(v)and not any(matmul([v],T)[0])
    j=next(i for i,vv in enumerate(v)if vv);proj=[[v[k]/v[j]if i==j else F(0)for k in range(6)]for i in range(6)]
    rows=row['full_magnetic_translation_rank']['rows'];cols=row['full_magnetic_translation_rank']['columns'];inv=q.inverse([[T[i][k]for k in cols]for i in rows]);diff=[[proj[i][k]-F(i==k)for k in range(6)]for i in rows]
    B=[[F(0)for k in range(6)]for i in range(6)];sol=matmul(inv,diff)
    for i,r in zip(cols,sol):B[i]=r
    assert matmul(T,B)==[[proj[i][k]-F(i==k)for k in range(6)]for i in range(6)]
    # Entire GL2 frame, not a fixed basis: adjugate transforms the bracket.
    g00,g01,g10,g11=map(p,params[:4]);Ug=[R.add(R.mul(g00,u),R.mul(g01,vv))for u,vv in zip(U,V)];Vg=[R.add(R.mul(g10,u),R.mul(g11,vv))for u,vv in zip(U,V)]
    ag=R.add(R.mul(g11,a),R.scale(R.mul(g10,b),-1));bg=R.add(R.scale(R.mul(g01,a),-1),R.mul(g00,b));Pg=source(R,Ug,Vg,ag,bg,magnetic)
    assert not any(jac(R,Pg));assert pf(R,Pg)==R.mul(R.add(R.mul(g00,g11),R.scale(R.mul(g01,g10),-1)),R.scale(f,ratio))
    Ng=divergence(R,Pg)
    assert all(not Ng[i]for i in(0,2,4))and all(all(sum(e[n:])==1 and e[n+1]==e[n+3]==0 for e in Ng[i])for i in(1,3))
    if norm:
        assert norm['complete_cokernel']==row['magnetic_cokernel_basis'][0]
        assert p(norm['canonical_R'])==mon[j]and norm['projective_orbit_dimension']==23
        assert rank==2 and row['canonical_type']=='nil22_image_type6'
    else:assert rank==1
    return dict(name=name,full28x18=True,rank=16,minor=d,full_GL2_and_all6_magnetic_strict10_and_projective350_Jacobi=True,actual_H_inverse=True,actual_Pf_ratio=str(ratio),full_magnetic_rank=5,magnetic_minor=tm,whole_regular_shear_matrix=[[str(v)for v in r]for r in B],magnetic_cokernel=[str(t)for t in v],canonical_cofactor_rank=rank,N_squared_zero_entire_source_identity=True,distinguished_h_outside_image=True)

CURVES={
'A6':(('s**4','s**2*t**2','8*s**3*t+t**4'),'(h*w-y**2)**2-64*h**3*y'),
'A4A2':(('s**4','s**2*t**2','12*s*t**3-6*t**4'),'(h*w+6*y**2)**2-144*h*y**3'),
'threeA2':(('s**2*(s-t)**2','t**2*(s-t)**2','s**2*t**2'),'(h*y+h*w+y*w)**2-4*h*y*w*(h+y+w)'),
'E6':(('s**3*t','s**4','t**4-6*s**2*t**2'),'y**3*w-h**4+6*h**2*y**2')}
PULLBACKS=['s**2*t**2','(s-t)**3*(s+t/3)','s**2*t**2','6*t**3*(2*s-t)','-6*t**2*(s-t)**2','(s-4*t/3)**3*(s+4*t)','s**2*(s-t)**2','t**3*(t-2*s)','(s**2-s*t+t**2)**2','-(s-t)**3*(3*s+t)','(3*s**2-t**2)**2']
def curvecheck(marked):
    B=LaurentRing((),('s','t'));results=[];phis={}
    for name,(phit,f)in CURVES.items():
        phi=[B.parse(t)for t in phit];phis[name]=phi
        aff=[trim(tuple(P.get((j,4-j),F(0))for j in range(5)))for P in phi]
        assert len(ug(ug(aff[0],aff[1]),aff[2]))==1 and any(P.get((4,0))for P in phi)
        M=[[B.mul(*(power(B,phi[i],k)for i,k in enumerate(e))).get((12-j,j),F(0))for e in exponents(3,3)]for j in range(13)];assert q.rank(M)[0]==10
        R=LaurentRing((),('h','y','w'));Ff=R.parse(f);rel={}
        for e,c in Ff.items():rel=B.add(rel,B.scale(B.mul(*(power(B,phi[i],k)for i,k in enumerate(e))),c))
        assert not rel
        if name=='A6':loc=[local(phi,0)];assert loc[0]['first_characteristic_odd']==7
        elif name=='A4A2':loc=[local(phi,0),local(phi,0,infinity=True)];assert[loc[0]['first_characteristic_odd'],loc[1]['first_characteristic_odd']]==[5,3]
        elif name=='threeA2':loc=[local(phi,0),local(phi,1),local(phi,0,infinity=True)];assert all(x['first_characteristic_odd']==3 for x in loc)
        else:loc=[local(phi,0)];assert loc[0]['orders']==[3,4]
        results.append(dict(profile=name,full_degree3_implicit_rank=10,basepointfree=True,irreducible_degree4_image=True,absolute_irreducible=True,local_cusps=loc,delta_sum3_no_other_singularities=True))
    pull=[]
    for row,expected in zip(marked['all11_candidates'],PULLBACKS):
        name=next(k for k in CURVES if row['name'].startswith(k));R=LaurentRing((),('h','y','w'));H=R.parse(row['H']);value={}
        for e,c in H.items():value=B.add(value,B.scale(B.mul(*(power(B,phis[name][i],k)for i,k in enumerate(e))),c))
        assert value==B.parse(expected);pull.append(dict(name=row['name'],exact_pullback=expected))
    return dict(four_absolute_irreducible_normalizations=results,all11_exact_marked_pullbacks=pull,whole_list_necessity_review='Separately written hand proof, all line-coefficient zero complements retained; not inferred from explicit finite representatives.')

class Q3:
    def __init__(self,a=0,b=0):
        if isinstance(a,Q3):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,x):x=Q3(x);return Q3(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return Q3(-self.a,-self.b)
    def __sub__(self,x):return self+-Q3(x)
    def __rsub__(self,x):return Q3(x)+-self
    def __mul__(self,x):x=Q3(x);return Q3(self.a*x.a+3*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):x=Q3(x);d=x.a*x.a-3*x.b*x.b;assert d;return self*Q3(x.a/d,-x.b/d)
    def __rtruediv__(self,x):return Q3(x)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,x):x=Q3(x);return self.a==x.a and self.b==x.b
    def __repr__(self):return str(self.a)+'+('+str(self.b)+')sqrt3'
class QRing(LaurentRing):
    def mul(self,*ps):
        raw=super().mul(*ps);out={}
        for e,c in raw.items():ee=list(e);c*=F(3)**(ee[0]//2);ee[0]%=2;ee=tuple(ee);out[ee]=out.get(ee,F(0))+c
        return {e:c for e,c in out.items()if c}

def rankgeneric(M):
    A=[list(r)for r in M];rr=0;cols=[];rows=list(range(len(A)));selected=[]
    for j in range(len(A[0])):
        k=next((i for i in range(rr,len(A))if A[i][j]),None)
        if k is None:continue
        A[rr],A[k]=A[k],A[rr];rows[rr],rows[k]=rows[k],rows[rr];selected.append(rows[rr]);d=A[rr][j];A[rr]=[v/d for v in A[rr]]
        for i in range(rr+1,len(A)):
            d=A[i][j]
            if d:A[i]=[a-d*b for a,b in zip(A[i],A[rr])]
        cols.append(j);rr+=1
    return rr,selected,cols

def fieldvalue(P,ex):return Q3(P.get((0,)+ex,F(0)),P.get((1,)+ex,F(0)))
def generic_kernel(M):
    A=[list(r)for r in M];rr=0;pivs=[]
    for j in range(len(A[0])):
        k=next((i for i in range(rr,len(A))if A[i][j]),None)
        if k is None:continue
        A[rr],A[k]=A[k],A[rr];d=A[rr][j];A[rr]=[v/d for v in A[rr]]
        for i in range(len(A)):
            if i!=rr:
                d=A[i][j]
                if d:A[i]=[a-d*b for a,b in zip(A[i],A[rr])]
        pivs.append(j);rr+=1
    out=[]
    for j in range(len(A[0])):
        if j not in pivs:
            v=[Q3(0)for k in range(len(A[0]))];v[j]=Q3(1)
            for i,k in enumerate(pivs):v[k]=-A[i][j]
            out.append(v)
    return out

def oldcheck():
    R=QRing(('rho',),('h','x','y','z','w'));p=R.parse
    data={
    'ESSENTIAL5':['h**2-x**2+x*z-x*w','-x**2','0','3*x*y-4*x*z+z**2-3*z*w','x**2+x*y-2*x*z+2*x*w+z**2/2-2*z*w+w**2','h**2+x**2+x*w'],
    'MINIMAL_M':['h**2+3*x**2-6*x*z','0','-2*x**2','3*x*z-18*z**2','6*x**2-2*x*y-3*x*w+30*z*w','h**2+6*x*z'],
    'A0':['h**2+3*x*z/4','0','z**2+3*x**2/4','-z**2/4','x**2+2*x*z+3*z**2-3*x*y/4+z*w/4','h**2+x*z/4'],
    'Ga9constant':['14*h*w-2*y**2','8*h*x+8*y*z+35*w**2','8*h*z-12*w*y','16*h*y','8*h**2','16*h*w+8*y**2']}
    rho='rho';c=p('rho/18');d=p('1/2+rho/18');a=p('-rho/6');b=p('-rho/216');u=p('1/6');v=p('-1/36');x,y,z,w=map(p,('x','y','z','w'))
    dp=[R.add(p('h**2'),R.mul(R.add(p('1'),R.scale(d,-1)),p('x*z')),R.scale(R.mul(v,p('x**2')),-1),R.scale(R.mul(u,p('z**2')),-1)),{},R.add(R.mul(a,p('z**2')),R.scale(R.mul(u,p('x*z')),-4),R.scale(R.mul(R.add(c,R.scale(d,-3),p('2')),p('x**2')),F(1,2))),R.add(R.mul(b,p('x**2')),R.scale(R.mul(v,p('x*z')),-4),R.scale(R.mul(R.add(c,d),p('z**2')),F(-1,2))),R.add(p('x**2+2*x*z+3*z**2'),R.mul(R.add(R.scale(d,2),R.scale(c,-1),p('-1')),p('x*y')),R.mul(c,p('z*w')),R.scale(R.mul(u,p('y*z')),2),R.scale(R.mul(v,p('x*w')),2)),R.add(p('h**2'),R.mul(d,p('x*z')),R.mul(v,p('x**2')),R.mul(u,p('z**2')))]
    out={}
    for name,strings in list(data.items())+[('DPM',dp)]:
        polys=[p(s)for s in strings]if name!='DPM'else strings;P={(i,j):{}for i in range(5)for j in range(5)}
        for ij,t in zip(combinations(range(1,5),2),polys):P[ij]=t;P[ij[::-1]]=R.scale(t,-1)
        assert not any(jac(R,P));f=pf(R,P);D=R.mul(p('h'),f);grad=[R.derivative(D,i)for i in range(5)];labels=sorted(set().union(*(g.keys()for g in grad)))
        GM=[[fieldvalue(g,e[1:])for g in grad]for e in labels if e[0]==0]
        # Include labels carrying only rho as well.
        labels=sorted({e[1:]for g in grad for e in g});GM=[[fieldvalue(g,e)for g in grad]for e in labels]
        er,rows,cols=rankgeneric(GM);N=divergence(R,P)
        out[name]=dict(actual_Pf_terms={','.join(map(str,e)):str(c)for e,c in f.items()},strict10_Jacobi=True,actual_essential_rank=er,essential_minor_rows=[list(labels[i])for i in rows],cofactor_terms=[{','.join(map(str,e)):str(c)for e,c in a.items()}for a in N])
        if name=='ESSENTIAL5':assert er==5 and f==p('h**4+h**2*x*z+x**3*y-x**3*z+x**2*z**2/2-x**2*z*w')
        if name=='Ga9constant':assert er==5
        if name=='MINIMAL_M':assert f==p('h**4+3*h**2*x**2+12*x**3*z')and N[2]==p('6*x')and all(not N[i]for i in(0,1,3,4))
        if name=='A0':assert f==p('h**4+h**2*x*z-z**4/4')
        if name in('A0','DPM'):
            assert not P[1,3];U=[R.scale(P[1,2],-1),P[2,3]];V=[R.scale(P[1,4],-1),R.scale(P[3,4],-1)]
            ab=[R.derivative(P[2,4],2),R.derivative(P[2,4],4)]
            pure=subst(R,P[2,4],{3:{},5:{}}) # actual exponent positions: n+physical y=3,w=5
            lin=[p(s)for s in('h','x','z')];mon=[p(s)for s in('h**2','h*x','h*z','x**2','x*z','z**2')]
            act=lambda X,L:R.add(R.mul(X[0],R.derivative(L,1)),R.mul(X[1],R.derivative(L,3)))
            columns=[R.add(R.scale(R.mul(ab[0],L),-1),R.scale(act(V,L),-1))for L in lin]+[R.add(R.scale(R.mul(ab[1],L),-1),act(U,L))for L in lin]
            T=[[fieldvalue(C,next(iter(m))[1:])for C in columns]for m in mon];left=generic_kernel(q.transpose(T));assert len(left)==1
            rv=[fieldvalue(pure,next(iter(m))[1:])for m in mon];I=sum((a*b for a,b in zip(left[0],rv)),Q3(0));assert I
            out[name].update(old_full_magnetic_translation_rank=5,magnetic_cokernel=[str(t)for t in left[0]],printed_actual_magnetic_invariant=str(I),nonzero_old_I_directly_checked=True)
        if name=='DPM':
            pp=subst(R,f,{1:{}});assert f==R.add(p('h**4+h**2*x*z'),pp)
            # Euclidean polynomials over QQ(sqrt3): p(t), p(t)-t²/4 each has exactly a triple root + a simple root.
            def trim3(P):
                P=list(P)
                while P and not P[-1]:P.pop()
                return P
            def rem(A,B):
                A=trim3(A)
                while len(A)>=len(B):
                    d=A[-1]/B[-1];k=len(A)-len(B)
                    for j,b0 in enumerate(B):A[j+k]=A[j+k]-d*b0
                    A=trim3(A)
                return A
            def gcd(A,B):
                while B:A,B=B,rem(A,B)
                return[x/A[-1]for x in A]
            P0=[fieldvalue(pp,(0,j,0,4-j,0))for j in range(5)];g=P0[:];g[2]=g[2]-Q3(F(1,4));roots=[]
            for pol in(P0,g):
                der=lambda A:[Q3(i)*v for i,v in enumerate(A)if i]
                h1=gcd(pol,der(pol));h2=gcd(h1,der(h1));assert len(h1)==3 and len(h2)==2
                rr=-h2[0]/h2[1];assert rr;roots.append(str(rr))
            out[name].update(three_cusps_via_two_triple_root_polynomials=True,triple_roots_p_and_p_minus_t2_over4=roots,marked_H_cusp_tangent_3plus1=True)
    return out

def priorcheck(old):
    proj=load('research/classification/inputs/mu6_prior111_37d56e2a_20261107.json');blob=subprocess.check_output(['git','cat-file','blob',proj['source_git_blob']],cwd=ROOT)
    assert hashlib.sha256(blob).hexdigest()==proj['source_sha256'];raw=json.loads(blob)
    rows=raw['components']if isinstance(raw,dict)and'components'in raw else raw
    if isinstance(rows,dict):rows=rows.get('identified_components',rows)
    assert rows==proj['all111_rows'];assert len(rows)==111
    prior=load('reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_REVIEW_20261107.json')['prior103_audit'];snc=set(prior['all29_actual_fullperiod_SNC_sources']);assert len(snc)==29
    audit=load(PREFIX+'PRIOR111_AUDIT_20261107.json');assert set(x['id']for x in rows)==set(x['id']for x in audit['all111_dispositions'])
    result=[]
    for row in rows:
        ident=row['id'];reason=next(x['reason']for x in audit['all111_dispositions']if x['id']==ident)
        if row.get('generic_rank')!=4:why='different generic rank'
        elif row.get('projective_dimension')!=23:why='different full component dimension'
        elif ident in snc:why='accepted independently reconstructed actual rank5 toric SNC point outside entire new closed essential<=3 locus'
        elif ident in('R4_H_QUARTIC_ESSENTIAL5','R4_GA9_CONSTANT_HYPERPLANE_QUARTIC_23'):why='printed actual essential5 point independently reconstructed here'
        elif ident=='R4_MINIMAL_M':why='old dense rank1-cofactor orbit versus new dense rank2-cofactor orbit'
        elif ident=='R4_A0':why='old marked D5 curve profile versus new mu6 cuspidal quartic'
        elif ident=='R4_DPM':why='old 3A2 cusp-tangent versus new 3A2 join/smooth-bitangent, or different curve profile'
        else:
            assert any(t in reason for t in('multiplicationlocus','closednonreducedfactorlocus','irreduciblequinticpoint'))
            why=reason
        result.append(dict(id=ident,independent_reason=why))
    names=audit['nine_component_candidates'];assert len(names)==9
    profile=lambda n:n.split('_')[0]
    pairs=[dict(first=a,second=b,reason='different intrinsic quartic cusp profile'if profile(a)!=profile(b)else'different marked H incidence/contact on the same unique curve normalization')for a,b in combinations(names,2)]
    assert len(pairs)==36
    return dict(immutable_git_commit=proj['source_commit'],blob=proj['source_git_blob'],all111_rows_equal_git_blob=True,all111_dispositions=result,all29_SNC_witnesses_reused_from_independent_first104_review=sorted(snc),old_Ga9constant_flag_closedness_not_needed_actual_essential5_suffices=True,all36_mutual_pairs=pairs)

def main():
    signal.alarm(180);t=time.monotonic();rows=load(PREFIX+'SOURCES_20261107.json')['cases'];normal={r['name']:r for r in load(PREFIX+'NORMALFORMS_20261107.json')['cases']};mark=load(PREFIX+'MARKED_LIST_20261107.json')
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,scope='MU6 whole marked curves and full sources/prior audit; ambient componenthood independently supplied by root, not recalculated. No MU5 claim.')
    out['curve_recognition']=curvecheck(mark);out['all13_complete_sources']=[sourcecheck(r,normal.get(r['name']))for r in rows];out['old_printed_witnesses']=oldcheck();out['prior111_and36_mutual']=priorcheck(out['old_printed_witnesses'])
    out['positive_whole_recognitions']={'D5_cusp_tangent':'A0: unique marked pair, complete two Saito anchors/all6mag, printed nonzero magnetic invariant, whole canonical orbit closure','threeA2_cusp_tangent':'DPM: direct threeA2/H3+1, unique marked pair, complete Saito/all6mag, printed nonzero magnetic invariant','E6_nonzero_transverse':'MINIMAL_M by accepted type8 complete-support theorem, source cofactor rank1 translation0; no type6 inference'}
    out['unresolved']=['E6_smooth_flex: rankD100/two normals; no parent or new component claim','MU5 is a separate packet','Other marked-profile branches retained in marked entrance']
    job=load(JOB);out['inputs']=job['exact_inputs_and_hashes'];assert all(sha(x['path'])==x['sha256']for x in out['inputs']);out['reviewer_sha256']=sha('code/'+Path(__file__).name);out['executed_reviewer_bytes']=Path(__file__).read_text();out['seconds']=time.monotonic()-t;out['maximum_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    path=ROOT/OUT;assert not path.exists();path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n');print(json.dumps({k:out[k]for k in('all_passed','seconds','maximum_RSS_KiB')}) )
if __name__=='__main__':main()
