#!/usr/bin/env python3
"""Independent Fraction verification of three whole q3 Anil3 orbits, two actual curves and faithful T39 membership.
No producer, producer engine, SymPy, NumPy or legacy arithmetic is imported.
"""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product, permutations
import hashlib,json,resource,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));QUADS=list(combinations(range(5),4))
N=15

def atom(c,*inds):
    e=[0]*N
    for i in inds:e[i]+=1
    return {tuple(e):F(c)} if c else {}

def add(*polys):
    out={}
    for p in polys:
        for e,c in p.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}

def scale(p,c):return {e:c*v for e,v in p.items() if c*v}

def mul(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            g=tuple(a+b for a,b in zip(e,f));out[g]=out.get(g,F(0))+c*d
    return {e:c for e,c in out.items() if c}

def der(p,k):
    out={}
    for e,c in p.items():
        if e[k]:
            f=list(e);f[k]-=1;out[tuple(f)]=c*e[k]
    return out

def parse(t):
    names=dict(zip(('h','x1','x2','x3','x4','t','L','n0','n1','n2','n3','n4','n5','n6','n7'),range(N)))
    def walk(p):
        if isinstance(p,ast.Constant) and isinstance(p.value,int):return atom(p.value)
        if isinstance(p,ast.Name):return atom(1,names[p.id])
        if isinstance(p,ast.UnaryOp):return scale(walk(p.operand),-1 if isinstance(p.op,ast.USub) else 1)
        if isinstance(p,ast.BinOp):
            if isinstance(p.op,ast.Add):return add(walk(p.left),walk(p.right))
            if isinstance(p.op,ast.Sub):return add(walk(p.left),scale(walk(p.right),-1))
            if isinstance(p.op,ast.Mult):return mul(walk(p.left),walk(p.right))
            if isinstance(p.op,ast.Div):
                den=walk(p.right);assert len(den)==1 and tuple([0]*N) in den
                return scale(walk(p.left),1/den[tuple([0]*N)])
            if isinstance(p.op,ast.Pow):
                out=atom(1)
                for _ in range(p.right.value):out=mul(out,walk(p.left))
                return out
        raise ValueError(ast.dump(p))
    return walk(ast.parse(t,mode='eval').body)

def antisymmetric(data):
    out={}
    for (i,j),v in data.items():out[i,j]=v;out[j,i]=scale(v,-1)
    return out

def decode(data):return antisymmetric({tuple(map(int,k.split(','))):parse(v) for k,v in data.items()})

def jacobi(P):
    J={}
    for i,j,k in combinations(range(5),3):
        J[i,j,k]=add(*(mul(P.get((a,l),{}),der(P.get((b,c),{}),l))
                    for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5)))
    return {q:add(*(scale(mul(atom(1,q[k]),J[q[:k]+q[k+1:]]),(-1)**k) for k in range(4))) for q in QUADS}

def monomials(d):
    out=[]
    for terms in combinations_with_replacement(range(5),d):
        e=[0]*5
        for k in terms:e[k]+=1
        out.append(tuple(e))
    return out

COLS=[(ij,e) for ij in PAIRS for e in monomials(2)]
ROWS=[(q,e) for q in QUADS for e in monomials(4)]

def vec(P):return [P.get(ij,{}).get(e+(0,)*(N-5),F(0)) for ij,e in COLS]
def jvec(J,ex=(0,0,0)):return [J[q].get(e+ex,F(0)) for q,e in ROWS]
def transpose(A):return [list(r) for r in zip(*A)]

def rank(A):
    B=[[F(v) for v in row] for row in A];r=0;piv=[]
    for j in range(len(B[0])):
        k=next((k for k in range(r,len(B)) if B[k][j]),None)
        if k is None:continue
        B[r],B[k]=B[k],B[r];d=B[r][j];B[r]=[v/d for v in B[r]]
        for k in range(r+1,len(B)):
            d=B[k][j]
            if d:B[k]=[u-d*v for u,v in zip(B[k],B[r])]
        piv.append(j);r+=1
    return r,piv

def determinant(A):
    B=[[F(v) for v in row] for row in A];out=F(1)
    for j in range(len(B)):
        k=next((k for k in range(j,len(B)) if B[k][j]),None)
        if k is None:return F(0)
        if k!=j:B[j],B[k]=B[k],B[j];out=-out
        d=B[j][j];out*=d
        for k in range(j+1,len(B)):
            a=B[k][j]/d
            if a:
                for l in range(j+1,len(B)):B[k][l]-=a*B[j][l]
                B[k][j]=0
    return out


def inverse(A):
    n=len(A);B=[[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k=next(k for k in range(j,n) if B[k][j]);B[j],B[k]=B[k],B[j]
        d=B[j][j];B[j]=[v/d for v in B[j]]
        for k in range(n):
            if k!=j:
                a=B[k][j];B[k]=[v-a*w for v,w in zip(B[k],B[j])]
    return [row[n:] for row in B]
def matmul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def substitute(p,lin):
    out={}
    for ex,c in p.items():
        v=atom(c)
        for i,k in enumerate(ex):
            for _ in range(k):v=mul(v,lin[i])
        out=add(out,v)
    return out

def evaluate(t,values):
    def f(n):
        if isinstance(n,ast.Constant):return F(n.value)
        if isinstance(n,ast.Name):return values[n.id]
        if isinstance(n,ast.UnaryOp):return -f(n.operand)
        a=f(n.left);b=f(n.right)
        if isinstance(n.op,ast.Add):return a+b
        if isinstance(n.op,ast.Sub):return a-b
        if isinstance(n.op,ast.Mult):return a*b
        if isinstance(n.op,ast.Div):return a/b
        if isinstance(n.op,ast.Pow):return a**int(b)
        raise ValueError(ast.dump(n))
    return f(ast.parse(t,mode='eval').body)


def decode_compact(d):
    return antisymmetric({tuple(map(int,k)):parse(v.replace('x0','h')) for k,v in d.items()})
def pdet(A):
    out={}
    for p in permutations(range(len(A))):
        v=atom((-1)**sum(p[i]>p[j] for i in range(len(A)) for j in range(i+1,len(A))))
        for i,j in enumerate(p):v=mul(v,A[i][j])
        out=add(out,v)
    return out

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_SUPPORT_20261107.json',
           'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_ORBITS_20261107.json',
           'reports/CLASSIFICATION_Q3_ANIL3_CURVES_20261107.json',
           'reports/CLASSIFICATION_Q3_ANIL3_UPPER_GRAPH_20261107.json',
           'legacy/v02/data/symplectic_slice_certificate.json','legacy/v02/data/symplectic_slice_equations.json',
           'legacy/v02/data/symplectic_slice_reduced_equations.json','reports/Q51_UNIVERSAL_SECTION_20261101.json']
    support,orbits,curves,graph,cert,model,reduced,toric=[json.loads((ROOT/p).read_text()) for p in paths]
    for name,rec in support['families'].items():
        base=decode_compact(rec['canonical_brackets']);target=decode_compact(rec['universal_brackets'])
        coords=[parse(v.replace('x0','h')) for v in rec['whole_family_coordinate_identity']]
        G=[[der(v,j) for j in range(5)] for v in coords]
        assert pdet(G)==parse(rec['coordinate_determinant'])
        lin=coords+[atom(1,j) for j in range(5,N)]
        for i,j in PAIRS:
            lhs=add(*(mul(mul(G[i][a],G[j][b]),base.get((a,b),{})) for a in range(5) for b in range(5)))
            assert lhs==substitute(target.get((i,j),{}),lin),(name,i,j)
        scaled=decode_compact(orbits['cases'][name]['canonical_brackets'])
        assert all(scaled.get(ij,{})==scale(base.get(ij,{}),F(orbits['cases'][name]['integer_scale'])) for ij in PAIRS)
    curve_checks={}
    for name,rec in curves['cases'].items():
        P=decode(rec['literal_actual_curve']);Pi=decode(rec['standard_unimodular_curve'])
        assert all(not v for v in jacobi(P).values())
        div=[add(*(der(P.get((i,j),{}),i) for i in range(5))) for j in range(5)]
        for i,j in PAIRS:
            expected=add(P.get((i,j),{}),scale(add(mul(atom(1,i),div[j]),scale(mul(atom(1,j),div[i]),-1)),F(-1,5)))
            assert expected==Pi.get((i,j),{})
        assert all(not add(*(der(Pi.get((i,j),{}),i) for i in range(5))) for j in range(5))
        for i,j,k in combinations(range(5),3):
            assert not add(*(mul(Pi.get((a,l),{}),der(Pi.get((b,c),{}),l)) for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5)))
        q=parse(rec['quadric']);c=parse(rec['actual_cubic']);D=mul(q,c)
        X=[parse(v) for v in rec['canonical_cofactor']]
        for j in range(5):assert add(*(mul(der(q,i),Pi.get((i,j),{})) for i in range(5)))==mul(q,X[j])
        A=[[der(v,j) for j in range(5)] for v in X]
        assert pdet([[add(atom(1,6) if i==j else {},scale(A[i][j],-1)) for j in range(5)] for i in range(5)])==parse(rec['cofactor_characteristic_polynomial'])
        H=[[der(der(q,i),j) for j in range(5)] for i in range(5)]
        assert all(not pdet([[H[i][j] for j in J] for i in I]) for I in QUADS for J in QUADS)
        assert pdet([r[:3] for r in H[:3]])==parse(rec['quadric_nonzero_3_minor'])
        for m in range(5):
            i,j,k,l=[v for v in range(5) if v!=m]
            pf=scale(add(mul(Pi.get((i,j),{}),Pi.get((k,l),{})),scale(mul(Pi.get((i,k),{}),Pi.get((j,l),{})),-1),mul(Pi.get((i,l),{}),Pi.get((j,k),{}))),(-1)**m)
            assert pf==der(D,m)
        lin=[atom(1,j) if j!=5 else {} for j in range(N)]
        original=decode_compact(orbits['cases'][rec['source_orbit']]['canonical_brackets'])
        assert all(substitute(P.get(ij,{}),lin)==original.get(ij,{}) for ij in PAIRS)
        c0=substitute(c,lin);assert c0==parse(rec['specialized_cubic'])
        # Direct primitive linear-polynomial tests certify absolute irreducibility.
        if name.endswith('_I'):
            assert der(c0,3)==parse('h**2/8')
            assert substitute(c0,[{}]+[atom(1,j) for j in range(1,N)])==parse('-(2*x1**3+5*x1**2*x4)/48')
        else:
            assert der(c0,4)==parse('-(10*h*x2+5*x1**2)/20736')
            assert substitute(c0,[atom(1,j) if j!=4 else {} for j in range(N)])==parse('-h**3/3456')
        curve_checks[name]={'all350_Jacobi':True,'canonical_gauge_and_strict_Jacobi':True,'all5_q_normal':True,'all5_Pfaffian_gradient':True,'all_t_quadric_rank':3,'characteristic_polynomial':rec['cofactor_characteristic_polynomial'],'endpoint_and_primitive_irreducible_cubic':True}
    P=decode_compact(orbits['cases']['upper']['canonical_brackets'])
    T=[[F(v) for v in row] for row in graph['translation']]
    B=[[F(v) for v in row] for row in graph['Darboux_frame']]
    A=[[F(i==j) for j in range(5)] for i in range(5)]
    A[0][1:]=list(map(F,graph['perspective']))
    M=matmul(A,matmul(B,T));Mi=inverse(M)
    lin=[add(*(atom(v,j) for j,v in enumerate(row))) for row in Mi]+[atom(1,i) for i in range(5,N)]
    pushed={ij:substitute(add(*(scale(P.get((a,b),{}),M[ij[0]][a]*M[ij[1]][b]) for a in range(5) for b in range(5))),lin) for ij in PAIRS}
    params=list(map(F,graph['all61_parameters']));assert len(params)==61
    def sparse_eval(data):
        out=F(0)
        for mon,num,den in data:
            v=F(num,den)
            for k in mon:v*=params[k]
            out+=v
        return out
    vv=vec({(1,2):atom(1,0,0),(3,4):atom(1,0,0)})
    for at,v in zip(cert['raw_free_indices'],params):vv[at]=v
    for e in model['eliminated_raw_variables']:vv[e['raw_index']]=sparse_eval(e['polynomial'])
    represented={ij:{} for ij in PAIRS}
    for ((ij,e),v) in zip(COLS,vv):
        if v:represented[ij][e+(0,)*(N-5)]=v
    Delta=antisymmetric({ij:add(pushed[ij],scale(represented[ij],-1)) for ij in PAIRS})
    for i,j,k in combinations(range(5),3):
        assert not add(mul(atom(1,i),Delta[j,k]),scale(mul(atom(1,j),Delta[i,k]),-1),mul(atom(1,k),Delta[i,j]))
    assert len(reduced['equations'])==134 and all(sparse_eval(e)==0 for e in reduced['equations'])
    vals={'u'+str(i+1):v for i,v in enumerate(params)}
    res=[params[k]-evaluate(q,vals) for k,q in zip(toric['weight2_u_indices_zero_based'],toric['quadratic_section_q'])]
    assert len(res)==45 and res==[0]*45

    result={'task_id':'E_reduced_nonnormal','all_passed':True,'method':'Independent standard-library Fraction sparse polynomial reconstruction; no producer/SymPy/S61 arithmetic imports',
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'three_universal_orbit_identities_all10_brackets':True,'three_coordinate_determinants_verified':True,
        'upper_full150_quotient_reconstruction':True,'upper_all134_equations':True,'upper_all45_T39_graph_equations':True,
        'curves':curve_checks,'scope':'All five q3 Anil3 cofactor types are covered by MPS, using the separately hand-reviewed complete-source reduction. Only actual families and orbits are needed, not the auxiliary weighted completed rings.',
        'seconds':time.monotonic()-started}
    (ROOT/'reports/CLASSIFICATION_Q3_ANIL3_SUPPORT_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'whole_sources_covered':3}))
if __name__=='__main__':main()
