#!/usr/bin/env python3
"""Independent whole mD=2 global-vector completion.

Only Fraction, sparse rational polynomials and existing independent reviewers.
No producer/CAS imports.  Full operators are reconstructed; rank lower bounds
and exhaustive maximal-minor walls are the weakest sufficient G certificate.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from functools import lru_cache
import ast, hashlib, json, resource, signal, time
from classification_unimodular_logarithmic_review import RF, ug, ud
import classification_degenerate_quadric_q4_nil31_I_review as q
import classification_plane_rank_zero_nonkernel_review as H
import classification_plane_scalar_review as B

ROOT=Path(__file__).resolve().parents[1]
PREF='reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_'
OUT='reports/CLASSIFICATION_PLANE_RANK_ZERO_M2_COMPLETION_REVIEW_20261107.json'
JOB='reports/CLASSIFICATION_PLANE_RANK_ZERO_M2_COMPLETION_REVIEW_JOB_20261107.json'
MAN='research/classification/inputs/plane_rank_zero_m2_full_completion_manifest_20261107.json'
FROZEN='research/classification/inputs/plane_rank_zero_m2_full_completion_frozen_jobs_20261107.json'
Z=(0,0,0);ONE={Z:F(1)}

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
            g=tuple(a+b for a,b in zip(e,f));out[g]=out.get(g,F(0))+c*d
    return {e:c for e,c in out.items() if c}
def pd(P,Q):
    assert Q
    rem=dict(P);out={};f=max(Q);c=Q[f]
    while rem:
        e=max(rem);assert all(a>=b for a,b in zip(e,f)),('nonexact polynomial division',e,f)
        g=tuple(a-b for a,b in zip(e,f));d=rem[e]/c
        out[g]=out.get(g,F(0))+d;rem=pa(rem,ps(pm({g:F(1)},Q),-d))
    return out

class R3:
    def __init__(self,value=0,den=None):
        if isinstance(value,R3):self.p,self.q=value.p,value.q;return
        self.p=dict(value) if isinstance(value,dict) else ({Z:F(value)} if value else {})
        self.q=dict(den) if den is not None else dict(ONE);assert self.q
        if not self.p:self.q=dict(ONE);return
        # Exact common monomial cancellation covers all graph denominators.
        ep=tuple(min(e[i] for e in self.p) for i in range(3));eq=tuple(min(e[i] for e in self.q) for i in range(3))
        g=tuple(min(a,b) for a,b in zip(ep,eq))
        if any(g):
            self.p={tuple(a-b for a,b in zip(e,g)):c for e,c in self.p.items()}
            self.q={tuple(a-b for a,b in zip(e,g)):c for e,c in self.q.items()}
        if len(self.q)==1 and Z in self.q:self.p=ps(self.p,1/self.q[Z]);self.q=dict(ONE)
    def __add__(self,v):
        v=R3(v)
        if self.q==v.q:return R3(pa(self.p,v.p),self.q)
        return R3(pa(pm(self.p,v.q),pm(v.p,self.q)),pm(self.q,v.q))
    __radd__=__add__
    def __neg__(self):return R3(ps(self.p,-1),self.q)
    def __sub__(self,v):return self+-R3(v)
    def __rsub__(self,v):return R3(v)+-self
    def __mul__(self,v):v=R3(v);return R3(pm(self.p,v.p),pm(self.q,v.q))
    __rmul__=__mul__
    def __truediv__(self,v):v=R3(v);assert v.p;return R3(pm(self.p,v.q),pm(self.q,v.p))
    def __rtruediv__(self,v):return R3(v)/self
    def __bool__(self):return bool(self.p)
    def __eq__(self,v):v=R3(v);return pm(self.p,v.q)==pm(v.p,self.q)
    def __repr__(self):return repr((self.p,self.q))

TAU=R3({(1,0,0):F(1)});ETA=R3({(0,1,0):F(1)});KAPPA=R3({(0,0,1):F(1)})
ENV={'tau':TAU,'eta':ETA,'kappa':KAPPA};U=H.previous.RHO
def scalar(s,env=ENV,field=R3):return H.read_scalar(s,env,field)
def subvalue(v,values,field):
    def one(P):
        out=field(0)
        for e,c in P.items():
            a=field(c)
            for x,k in zip(values,e):
                for _ in range(k):a=a*x
            out=out+a
        return out
    return one(v.p)/one(v.q)
def submatrix(M,values,field):
    # Sparse full350-row matrices repeat very few coefficient functions.  Cache
    # structurally equal fractions locally, including the overwhelmingly common0.
    cache={}
    def entry(v):
        key=(tuple(sorted(v.p.items())),tuple(sorted(v.q.items())))
        if key not in cache:cache[key]=subvalue(v,values,field)
        return cache[key]
    return [[entry(v)for v in row]for row in M]
def scaled(M,a):return [[a*v for v in row] for row in M]
def small_encode(v):
    if isinstance(v,R3):return dict(p=[[list(e),str(c)]for e,c in sorted(v.p.items())],q=[[list(e),str(c)]for e,c in sorted(v.q.items())])
    return H.encode(v)

def determinant(M):
    """Perfect matching/SCC blocks, then exact polynomial Bareiss per block."""
    if not M:return R3(1)
    if not isinstance(M[0][0],R3):
        if all(isinstance(v,RF)and len(v.p)<=1 and v.q==(F(1),)for row in M for v in row):
            return RF(H.det_generic([[v.p[0]if v.p else F(0)for v in row]for row in M]))
        return H.det_generic(M)
    n=len(M);assert all(len(row)==n for row in M)
    # A selected minor is polynomial. Exact division also clears graph fractions.
    A=[[pd(v.p,v.q) if v.q!=ONE else v.p for v in row]for row in M]
    match={}
    def augment(i,seen):
        for j in range(n):
            if not A[i][j] or j in seen:continue
            seen.add(j)
            if j not in match or augment(match[j],seen):match[j]=i;return True
        return False
    if not all(augment(i,set())for i in range(n)):return R3(0)
    perm=[next(j for j in match if match[j]==i)for i in range(n)]
    sign=(-1)**sum(perm[i]>perm[j]for i in range(n)for j in range(i+1,n))
    graph=[[j for j in range(n)if A[i][perm[j]]]for i in range(n)]
    counter=[0];idx={};low={};stack=[];active=set();parts=[]
    def visit(i):
        idx[i]=low[i]=counter[0];counter[0]+=1;stack.append(i);active.add(i)
        for j in graph[i]:
            if j not in idx:visit(j);low[i]=min(low[i],low[j])
            elif j in active:low[i]=min(low[i],idx[j])
        if low[i]==idx[i]:
            part=[]
            while True:
                j=stack.pop();active.remove(j);part.append(j)
                if j==i:break
            parts.append(sorted(part))
    for i in range(n):
        if i not in idx:visit(i)
    out={Z:F(sign)}
    for part in parts:
        C=[[dict(A[i][perm[j]])for j in part]for i in part];prev=ONE;sgn=1
        for j in range(len(C)-1):
            k=next((i for i in range(j,len(C))if C[i][j]),None)
            if k is None:return R3(0)
            if k!=j:C[j],C[k]=C[k],C[j];sgn=-sgn
            pivot=C[j][j]
            for k in range(j+1,len(C)):
                for l in range(j+1,len(C)):
                    C[k][l]=pd(pa(pm(C[k][l],pivot),ps(pm(C[k][j],C[j][l]),-1)),prev)
                C[k][j]={}
            prev=pivot
        out=pm(out,ps(C[-1][-1],sgn))
    return R3(out)

def minor(M,C,env,field):
    r=C['rank'];assert C['nullity']==len(M[0])-r
    rows,cols=C['minor_rows'],C['minor_columns'];assert len(rows)==len(cols)==r
    assert len(set(rows))==len(set(cols))==r
    v=determinant([[M[i][j]for j in cols]for i in rows])if r else field(1)
    assert v and v==scalar(C['minor_determinant'],env,field)
    return v

def triple(Ms,C,env,field,offset,allow_failure=False):
    minors=[minor(M,C[k],env,field)for M,k in zip(Ms,('D1','D2','actual_shear'))]
    d=39-C['D1']['rank']-C['D2']['rank'];bound=offset+d
    assert d==C['higher_fiber_upper_bound'] and bound==C['PGL_image_upper_bound']
    if not allow_failure:assert bound<=24
    return dict(ranks=[C[k]['rank']for k in ('D1','D2','actual_shear')],higher_fiber_upper=d,image_upper=bound),minors[0]*minors[1]*minors[2]

def factor_product(product,analysis,env,field):
    assert product==scalar(analysis['exception_polynomial'],env,field)
    fact=field(1);factors=[]
    for text,e in analysis['all_irreducible_factors']:
        v=scalar(text,env,field);assert v and e>=1
        factors.append((text,v))
        for _ in range(e):fact=fact*v
    if field is R3:
        a=pd(product.p,product.q);b=pd(fact.p,fact.q);quot=pd(a,b)
        assert set(quot)=={Z} and quot[Z]
    else:
        quot=product/fact;assert quot.q==(F(1),)and len(quot.p)==1 and quot.p[0]
    return factors

def model(name,env=ENV,field=R3):
    p=lambda text:H.parse(text,env,field)
    Hbase=[p('a'),p('b'),p('-2*c'),{},{}];N=[p('b'),{},{},{},{}]
    J=[p('b'),p('c'),{},{},{}];J2=[p('c'),{},{},{},{}]
    X=lambda T,S:H.Padd(q.antisymmetric({ij:q.add(q.scale(T.get(ij,{}),S))for ij in q.PAIRS}))
    tv=lambda coeff,V:[q.scale(v,scalar(coeff,env,field))for v in V]
    plus=lambda *Vs:[q.add(*(v[i]for v in Vs))for i in range(5)]
    if name=='SS_MIXH':T=plus(tv('tau',Hbase),N);V=Hbase
    elif name=='SS_MIXN':T=plus(Hbase,tv('tau',N));V=N
    elif name=='SS_NILREG':T=plus(tv('tau',J),J2);V=J
    elif name=='NIL_B0_MIX':T=plus(Hbase,tv('tau',N));V=plus(Hbase,N)
    elif name=='SS_MIX':T=plus(tv('tau',Hbase),tv('eta',N));V=plus(Hbase,N)
    else:
        assert name=='NIL_B0_CARTAN'
        T=[p('-(1+tau)*a/3'),p('(2-tau)*b/3'),p('(2*tau-1)*c/3'),{},{}]
        V=[p('-(1+eta)*a/3'),p('(2-eta)*b/3'),p('(2*eta-1)*c/3'),{},{}]
    normal=[{},{},{},p('z'),p('-w')]if name.startswith('SS')else[{},{},{},{},p('z')]
    En=[{},{},{},p('z'),p('w')]
    kval='kappa'if name in ('SS_MIX','NIL_B0_CARTAN')else'eta'
    K=p(kval+('*z*w'if name.startswith('SS')else'*z**2'))
    return H.Padd(H.wedge(T,En),H.wedge(V,normal),q.antisymmetric({(3,4):K}))

def bases():
    first=[]
    for normal in (3,4):
        for ij in combinations(range(3),2):
            for i in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,i)}))
    for normal,ex in [(4,(3,3)),(3,(3,4)),(4,(3,4)),(3,(4,4))]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*ex)}))
    second=[q.antisymmetric({ij:q.atom(1,*ex)})for ij in combinations(range(3),2)for ex in [(3,3),(3,4),(4,4)]]
    return first,second

def full_operators(c,first,second):
    P=model(c['type']);assert not any(q.jacobi(P).values())
    expected=q.antisymmetric({ast.literal_eval(k):H.parse(v,ENV,R3)for k,v in c['normal_form_P0'].items()})
    assert q.vec(P)==q.vec(expected) and H.divisor(P)==H.parse(c['actual_D'],ENV,R3)
    Ms=[H.convert_matrix(H.previous.operator(P,B.convert_basis(basis,R3)),R3)for basis in (first,second)]
    S=H.shear(P,first,R3);Ms.append(S)
    keys=('complete_D1_350_by30','complete_D2_350_by9','complete_actual_shear_30_by6')if 'complete_D1_350_by30'in c else('full_D1_350_by30','full_D2_350_by9','full_actual_shear_30_by6')
    for M,key in zip(Ms,keys):assert M==H.matrix(c[key],ENV,R3)
    assert all(not sum((v*w for v,w in zip(row,col)),R3(0))for row in Ms[0]for col in zip(*S))
    return P,Ms

def root_analysis(Ms,a,offset,domain,covered=None):
    # QQ(tau), all roots here are rational; equality of the whole factor product
    # proves completeness rather than relying on a list of sampled values.
    gen,prod=triple(Ms,a,{'tau':U},RF,offset,allow_failure=True)
    factors=factor_product(prod,a,{'tau':U},RF)
    if not covered:assert gen['image_upper']<=24
    roots=[];seen=set()
    for text,f in factors:
        assert f.q==(F(1),)and len(f.p)==2
        value=-f.p[0]/f.p[1];assert value not in seen;seen.add(value)
        matches=[r for r in a['all_exceptional_roots']if F(r['parameter'])==value];assert len(matches)==1
        r=matches[0];current=[H.convert_matrix(M,lambda v:at_rf(v,value))for M in Ms]
        res,_=triple(current,r,{'tau':RF(value)},RF,offset-1,allow_failure=True)
        outside=not at_rf(domain,value);assert outside==r['excluded_multiplicity_higher']
        elsewhere=covered is not None and not at_rf(covered,value)
        if not outside and not elsewhere:assert res['image_upper']<=24
        res.update(parameter=str(value),excluded_multiplicity_higher=outside,covered_by_nil_line=elsewhere)
        roots.append(res)
    assert len(roots)==len(a['all_exceptional_roots'])
    return dict(generic=gen,roots=roots,whole_factor_product_verified=True)

@lru_cache(maxsize=250000)
def at_rf(value,t):
    # Every final root in this packet is rational.  Horner in Fraction avoids
    # millions of redundant QQ(t)-gcd normalizations at constant specializations.
    t=F(t)
    def evaluate(P):
        out=F(0)
        for c in reversed(P):out=out*t+c
        return out
    return RF(evaluate(value.p)/evaluate(value.q))

def r3_linear_coefficients(v,n):
    assert v.q==ONE and all(sum(e)<=1 for e in v.p)
    return [[v.p.get(tuple(int(i==j)for i in range(3)),F(0))for j in range(n)],v.p.get(Z,F(0))]

def curve_analysis(Ms,wall,domain,covered=None):
    env1={'tau':U};maps={k:scalar(v,env1,RF)for k,v in wall['parameterization'].items()}
    t,s=maps['tau'],maps['eta'];f=scalar(wall['equation']);assert not subvalue(f,(t,s,RF(0)),RF)
    graphA=RF(1)
    if 'graph_A'in wall:
        A=scalar(wall['graph_A'],env1,RF);BB=scalar(wall['graph_B'],env1,RF)
        graphA=A
        assert A.q==(F(1),)and BB.q==(F(1),)and len(ug(A.p,BB.p))==1
        # The normalized parameter itself is one original coordinate, so graph
        # inversion is global wherever A!=0; gcd(A,B)=1 leaves no vertical points.
        assert t==U or s==U
        other=s if t==U else t;assert not A*other+BB
    else:
        co,constant=r3_linear_coefficients(f,2);assert any(co)
        assert all(v.q==(F(1),)and len(v.p)<=2 for v in (t,s)) and(t==U or s==U)
    pull=[submatrix(M,(t,s,RF(0)),RF)for M in Ms]
    dd=subvalue(domain,(t,s,RF(0)),RF);cc=subvalue(covered,(t,s,RF(0)),RF)if covered is not None else None
    flag=wall.get('excluded_multiplicity_higher',wall.get('excluded_by_graph_domain_or_higher_multiplicity',False))
    assert (not dd)==flag
    if not dd:return dict(equation=wall['equation'],excluded_whole_curve=True)
    if covered is not None and not cc:return dict(equation=wall['equation'],whole_curve_covered_by_nil_line=True)
    denominator=scalar(wall.get('rank_equivalent_denominator_clearing','1'),env1,RF);assert denominator
    current=[scaled(M,denominator)for M in pull]
    for M,key in zip(current,('full_cleared_D1_350_by30','full_cleared_D2_350_by9','full_cleared_actual_shear_30_by6')):
        if key in wall:assert M==H.matrix(wall[key],env1,RF)
    # denominatorzero is outside the true graph; include it in the excluded
    # function for root bookkeeping, exactly as in the producer contract.
    a=wall['analysis'];res=root_analysis(current,a,15,RF(dd.p)*denominator*graphA,cc)
    res['equation']=wall['equation'];return res

def plane_analysis(Ms,a,domain,covered=None,frontier=None):
    generic,prod=triple(Ms,a,ENV,R3,16,allow_failure=True)
    assert generic['image_upper']<=24
    factors=factor_product(prod,a,ENV,R3)
    walls=a.get('all_exceptional_lines',a.get('all_whole_curve_walls',[]));unsupported=a.get('unsupported_whole_curves',[])
    assert {x['equation']for x in walls}|{x['whole_curve_equation']for x in unsupported}=={t for t,f in factors}
    results=[curve_analysis(Ms,w,domain,covered)for w in walls]
    for x in unsupported:
        assert frontier=='hyperbola'and x['whole_curve_equation']=='eta*tau - tau + 3'
    return dict(generic=generic,whole_curves=results,retained_frontier=unsupported)

def uadd(*Ps):
    out={}
    for P in Ps:
        for e,c in P.items():out[e]=out.get(e,RF(0))+c
    return {e:c for e,c in out.items()if c}
def uscale(P,c):return {e:c*v for e,v in P.items()if c*v}
def umul(P,Q):
    out={}
    for e,c in P.items():
        for f,d in Q.items():
            g=tuple(a+b for a,b in zip(e,f));out[g]=out.get(g,RF(0))+c*d
    return {e:c for e,c in out.items()if c}
def upoly(text,n,env):
    zero=(0,)*n
    def read(a):
        if isinstance(a,ast.Constant):return {zero:RF(a.value)}if a.value else{}
        if isinstance(a,ast.Name):
            if a.id in env:return {zero:env[a.id]}if env[a.id]else{}
            assert a.id.startswith('u');j=int(a.id[1:]);assert j<n
            return {tuple(int(i==j)for i in range(n)):RF(1)}
        if isinstance(a,ast.UnaryOp):return uscale(read(a.operand),-1 if isinstance(a.op,ast.USub)else 1)
        x=read(a.left)
        if isinstance(a.op,ast.Pow):
            out={zero:RF(1)}
            for _ in range(a.right.value):out=umul(out,x)
            return out
        y=read(a.right)
        if isinstance(a.op,ast.Add):return uadd(x,y)
        if isinstance(a.op,ast.Sub):return uadd(x,uscale(y,-1))
        if isinstance(a.op,ast.Mult):return umul(x,y)
        if isinstance(a.op,ast.Div):assert set(y)=={zero};return uscale(x,RF(1)/y[zero])
        raise ValueError(ast.dump(a))
    return read(ast.parse(text,mode='eval').body)

def nil_line_one(c,Ms,first):
    value=c['parameter'];env={'tau':U if value is None else RF(F(value))}
    current=Ms if value is None else[H.convert_matrix(M,lambda v:H.rf_at(v,env['tau']))for M in Ms]
    D1,D2,S=current;certs={}
    for M,k in zip(current,('D1','D2','actual_shear')):certs[k],_=H.certificate(M,c[k],env,RF)
    rr=c['actual_shear']['minor_columns'];r=len(rr)
    slice_=H.matrix(c['full_firstweight_slice_30_by_n'],env,RF);n=len(slice_[0]);assert len(slice_)==30
    assert n+r==c['D1']['nullity']
    assert all(not sum((a*b for a,b in zip(row,col)),RF(0))for row in D1 for col in zip(*slice_))
    basis=[[S[i][j]for j in rr]+slice_[i]for i in range(30)]
    rows=c['shear_plus_slice_basis_rows'];square=[basis[i]for i in rows]
    assert square==H.matrix(c['shear_plus_slice_basis_minor'],env,RF)
    basisdet=H.det_generic(square);assert basisdet
    vectors=[]
    for col in zip(*slice_):
        vectors.append(q.antisymmetric({ij:q.add(*(q.scale(B.get(ij,{}),a)for B,a in zip(first,col)))for ij in q.PAIRS}))
    selfJ=[q.jvec(q.jacobi(P))for P in vectors];Q=[{}for _ in range(350)]
    for j in range(n):
        for k in range(j,n):
            vv=selfJ[j]if j==k else[a-b-c for a,b,c in zip(q.jvec(q.jacobi(H.Padd(vectors[j],vectors[k]))),selfJ[j],selfJ[k])]
            ex=tuple(int(i==j)+int(i==k)for i in range(n))
            for i,v in enumerate(vv):
                if v:Q[i][ex]=RF(v)
    assert Q==[upoly(t,n,env)for t in c['full_Q2_350']]
    rows,cols=c['D2']['minor_rows'],c['D2']['minor_columns']
    inv=H.matrix(c['D2_pivot_inverse'],env,RF)
    square=[[D2[i][j]for j in cols]for i in rows]
    assert q.matmul(square,inv)==[[RF(int(i==j))for j in range(len(cols))]for i in range(len(cols))]
    corrections=[uadd(*(uscale(Q[i],v)for i,v in zip(rows,row)))for row in inv]
    residual=[uadd(Q[i],*(uscale(P,-D2[i][j])for j,P in zip(cols,corrections)))for i in range(350)]
    assert residual==[upoly(t,n,env)for t in c['full_compatibility_350']]
    chain=c['exact_nonzero_or_gcd_chain'];forms=[]
    for item in chain:
        p=upoly(item['polynomial'],n,env);assert p==residual[item['row']];forms.append(p)
    # Every supplied chain starts with a monomial. Thus its common irreducible
    # factors can only be variables; minimum exponents prove gcd=1 exactly.
    assert len(forms[0])==1
    common=tuple(min(e[j]for P in forms for e in P)for j in range(n));assert not any(common)
    assert c['proven_compatibility_height_lower_bound']==2
    higher=r+n-2+c['D2']['nullity'];base=10 if value is None else 9;bound=base+5+higher
    assert higher==c['full_higher_dimension_upper']and bound==c['PGL_image_upper']and bound<=24
    extras=[]
    if value is None:
        br=c['binary_resultant_certificate'];p=upoly(br['first_binary_quadric'],n,env);qq=upoly(br['second_binary_quadric'],n,env)
        assert p in residual and qq in residual and n==2
        mons=((2,0),(1,1),(0,2));a,b,cc=[p.get(e,RF(0))for e in mons];d,e,f=[qq.get(e,RF(0))for e in mons]
        syl=[[a,b,cc,RF(0)],[RF(0),a,b,cc],[d,e,f,RF(0)],[RF(0),d,e,f]]
        assert syl==H.matrix(br['Sylvester_matrix'],env,RF)
        resultant=H.det_generic(syl);assert resultant and resultant==scalar(br['nonzero_binary_resultant'],env,RF)
        extras=[resultant,basisdet]+[v for row in slice_+inv for v in row]
        extras.extend(scalar(c[k]['minor_determinant'],env,RF)for k in ('D1','D2','actual_shear'))
    return dict(parameter=value,full_operator_kernel_and_slice_certificates=certs,shear_rank=r,slice_dimension=n,
                all350_true_cokernel_quadrics=True,monomial_first_gcd_certificate=chain,
                compatibility_height_at_least2=True,higher_fiber_upper=higher,image_upper=bound),extras

def nil_line(packet,Ms,first):
    data=[nil_line_one(packet['generic'],Ms,first)];roots=[F(c['parameter'])for c in packet['all_exceptional_points']]
    assert set(roots)=={F(-3),F(-1),F(-1,3),F(1,3),F(1),F(3)}and len(roots)==6
    seen=set()
    for v in data[0][1]:
        for polynomial in (v.p,v.q):
            rem=polynomial
            if not any(rem):continue
            for rho in roots:
                factor=(-rho,F(1))
                while len(rem)>1:
                    quotient,remainder=ud(rem,factor)
                    if any(remainder):break
                    rem=quotient;seen.add(rho)
            assert len(rem)==1,('uncovered generic denominator/minor/resultant root',rem)
    assert seen==set(roots)
    exception=scalar(packet['generic']['generic_exception_polynomial'],{'tau':U},RF)
    prod=RF(1)
    for rho in roots:prod=prod*(U-rho)
    assert exception==prod
    for c in packet['all_exceptional_points']:data.append(nil_line_one(c,Ms,first))
    return dict(generic=data[0][0],all_six_rank_basis_resultant_roots=[a for a,b in data[1:]],all_denominator_complements_covered=True)

def three_analysis(c,Ms):
    a=c['analysis'];res,prod=triple(Ms,a,ENV,R3,17)
    factors=factor_product(prod,a,ENV,R3)
    walls=a['all_whole_plane_walls'];unsupported=a['unsupported_whole_surfaces']
    assert {w['equation']for w in walls}|{w['whole_surface_equation']for w in unsupported}=={t for t,f in factors}
    results=[]
    domain=TAU-ETA
    for w in walls:
        f=scalar(w['equation']);co,const=r3_linear_coefficients(f,3);assert any(co)
        maps={k:scalar(v)for k,v in w['parameterization'].items()};vals=[maps[k]for k in ('tau','eta','kappa')]
        assert not subvalue(f,vals,R3)
        # Each plane parameterization is affine of full rank two; in these
        # canonical maps the two unsolved coordinates are literally tau,eta.
        assert vals.count(TAU)>=1 and vals.count(ETA)>=1
        assert all(v.q==ONE and all(e[2]==0 and sum(e)<=1 for e in v.p)for v in vals)
        dd=subvalue(domain,vals,R3);assert (not dd)==w['excluded_multiplicity_higher']
        if not dd:results.append(dict(equation=w['equation'],excluded_whole_plane=True));continue
        pulled=[submatrix(M,vals,R3)for M in Ms]
        for M,key in zip(pulled,('full_D1_350_by30','full_D2_350_by9','full_actual_shear_30_by6')):assert M==H.matrix(w[key],ENV,R3)
        rr=plane_analysis(pulled,w['analysis'],dd);rr['equation']=w['equation'];results.append(rr)
    if unsupported:
        assert c['type']=='NIL_B0_CARTAN'and len(unsupported)==2
        assert {u['whole_surface_equation']for u in unsupported}=={'eta*tau - 3*eta + 3*kappa*tau + tau','eta*tau + 2*eta + 3*kappa*tau - 3*kappa - 2*tau - 1'}
    return dict(generic=res,all_whole_plane_walls=results,retained_surfaces_closed_by_separate_receipt=unsupported)

def cartan_graph(packet,P,Ms):
    t,s,k=TAU,ETA,KAPPA;f1=s*t-3*s+3*k*t+t;f2=s*t+2*s+3*k*t-3*k-2*t-1
    assert scalar(packet['first_surface'])==f1 and scalar(packet['second_surface'])==f2
    assert subvalue(f1,(1-t,1-s,-k),R3)==f2
    # Literal actual base permutation and tensor -1 on the entire P0.
    sigma=[1,0,2,3,4];changed={}
    for i,j in q.PAIRS:
        old=P.get((sigma[i],sigma[j]),{})
        changed[i,j]={tuple(e[sigma[a]]for a in range(5))+e[5:]:-v for e,v in old.items()}
    target=model('NIL_B0_CARTAN',{'tau':1-t,'eta':1-s,'kappa':-k},R3)
    assert q.vec(q.antisymmetric(changed))==q.vec(target)
    graph=scalar(packet['first_surface_graph']);assert not subvalue(f1,(t,s,graph),R3)
    assert subvalue(f1,(R3(0),s,k),R3)==-3*s and subvalue(f2,(R3(1),s,k),R3)==3*(s-1)
    pulled=[scaled(submatrix(M,(t,s,graph),R3),t)for M in Ms]
    for M,key in zip(pulled,('full_cleared_D1_350_by30','full_cleared_D2_350_by9','full_cleared_actual_shear_30_by6')):assert M==H.matrix(packet[key],ENV,R3)
    result=plane_analysis(pulled,packet['analysis'],t*(t-s))
    result.update(actual_swap_tensor_inverse=True,both_vertical_complements_exactly_outside_m2=True)
    return result

def main():
    started=time.monotonic();signal.alarm(1200);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists()
    job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    manifest=json.loads((ROOT/MAN).read_text());frozen=json.loads((ROOT/FROZEN).read_text())
    reports={k:json.loads((ROOT/(PREF+k+'_20261107.json')).read_text())for k in ('TWOPARAMETER','HYPERBOLA','NIL_LINE','THREEPARAMETER','CARTAN_SURFACE')}
    source=reports['CARTAN_SURFACE']['producer_bytes'];assert hashlib.sha256(source.encode()).hexdigest()==frozen['runner_sha256']
    snapshots={x['path']:x for x in frozen['snapshots']}
    for row in manifest['files']:
        p=row['path']
        if p==frozen['runner_path']:assert row['sha256']==frozen['runner_sha256']
        elif p in snapshots:assert hashlib.sha256(snapshots[p]['bytes_utf8'].encode()).hexdigest()==row['sha256']==snapshots[p]['sha256']
        else:assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==row['sha256'],p
    first,second=bases();twos=[];raw={}
    for c in reports['TWOPARAMETER']['cases']:
        P,Ms=full_operators(c,first,second);raw[c['type']]=(P,Ms)
        domain=scalar(c['excluded_polynomial']);covered=ETA if c['type']=='SS_NILREG'else None
        r=plane_analysis(Ms,c['analysis'],domain,covered,'hyperbola'if c['type']=='NIL_B0_MIX'else None)
        r['type']=c['type'];twos.append(r);print(c['type'],'full operators and whole walls PASS',flush=True)
    original=raw['NIL_B0_MIX'][1]
    vals=(U,RF(1)-RF(3)/U,RF(0));hyper=[scaled(submatrix(M,vals,RF),U)for M in original]
    for M,key in zip(hyper,('full_cleared_D1_350_by30','full_cleared_D2_350_by9','full_cleared_actual_shear_30_by6')):assert M==H.matrix(reports['HYPERBOLA'][key],{'tau':U},RF)
    hyperres=root_analysis(hyper,reports['HYPERBOLA']['analysis'],15,U*(U-1))
    assert not subvalue(scalar('eta*tau-tau+3'),vals,RF)
    nilMs=[submatrix(M,(U,RF(0),RF(0)),RF)for M in raw['SS_NILREG'][1]]
    nilres=nil_line(reports['NIL_LINE'],nilMs,B.convert_basis(first,RF));print('whole nil K0 true compatibility PASS',flush=True)
    threes=[];three_raw={}
    for c in reports['THREEPARAMETER']['cases']:
        P,Ms=full_operators(c,first,second);three_raw[c['type']]=(P,Ms)
        r=three_analysis(c,Ms);r['type']=c['type'];threes.append(r);print(c['type'],'all whole surfaces/curves/points PASS',flush=True)
    cartan=cartan_graph(reports['CARTAN_SURFACE'],*three_raw['NIL_B0_CARTAN']);print('both Cartan surfaces whole graph PASS',flush=True)
    # Only audit AST compatibility of the immutable historical bytes; do not
    # bind or execute the mutable source or any mutable job paths.
    functions=lambda s:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(s).body if isinstance(n,ast.FunctionDef)}
    prior=None
    for k in ('TWOPARAMETER','HYPERBOLA','NIL_LINE','THREEPARAMETER','CARTAN_SURFACE'):
        d=reports[k];assert hashlib.sha256(d['producer_bytes'].encode()).hexdigest()==d['producer_sha256']
        now=functions(d['producer_bytes'])
        if prior is not None:assert all(now.get(n)==v for n,v in prior.items())
        prior=now
    result=dict(task_id='G',all_passed=True,
        theorem='Entire rank-zero-plane mD2 source is contained in the accepted MPS40 union or finitely many projective images of dimension<=24; no E/fullG claim.',
        independent_full_operators_and_rank_lower_bounds=True,no_CAS_or_producer_imports=True,
        two_parameter_whole_sources=twos,whole_hyperbola=hyperres,whole_nil_line=nilres,
        three_parameter_whole_sources=threes,both_Cartan_surfaces=cartan,
        retained_old_false_reports_immutable=True,mutable_source_and_jobs_not_bound=True,
        accepted_prior_low21_not_recomputed=True,
        complements=['mD>=3','pointwise-dependent T,V/global common-image high multiplicity','rank-two restriction on a surface','all other G/E obligations'],
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS',OUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)
if __name__=='__main__':main()

