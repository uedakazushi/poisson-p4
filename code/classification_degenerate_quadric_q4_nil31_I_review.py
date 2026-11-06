#!/usr/bin/env python3
"""Independent Fraction reconstruction of the full three-parameter smooth family and faithful T40 graph membership.
No producer, producer engine, SymPy, NumPy or legacy arithmetic is imported.
"""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
import hashlib,json,resource,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));QUADS=list(combinations(range(5),4))
N=8

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
    names=dict(zip(('h','x1','x2','x3','x4','tau0','tau1','tau2'),range(N)))
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

def vec(P):return [P.get(ij,{}).get(e+(0,0,0),F(0)) for ij,e in COLS]
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

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_WEIGHTED_20261107.json',
           'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_GRAPH_20261107.json',
           'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json',
           'legacy/v02/data/symplectic_slice_certificate.json','legacy/v02/data/symplectic_slice_equations.json',
           'legacy/v02/data/symplectic_slice_reduced_equations.json','reports/Q52_UNIVERSAL_SECTION_20261101.json']
    saved,graph,comp,cert,model,reduced,toric=[json.loads((ROOT/p).read_text()) for p in paths]
    P=decode(saved['canonical_brackets']);normals=[decode(r) for r in saved['normal_brackets']]
    case=next(c for c in comp['cases'] if c['source']=='nil31_long_I')
    for ij in PAIRS:
        assert P.get(ij,{})==parse(case['point_brackets'].get(str(ij),'0').replace('x0','h'))
    assert all(not q for q in jacobi(P).values())
    Dcols=[]
    for ij,e in COLS:
        v={key:dict(value) for key,value in P.items()};inc={e+(0,0,0):F(1)}
        v[ij]=add(v.get(ij,{}),inc);v[ij[::-1]]=scale(v[ij],-1);Dcols.append(jvec(jacobi(v)))
    D=transpose(Dcols);G=[]
    for a,b in product(range(5),repeat=2):
        V={};O={}
        for i,j in PAIRS:
            V[i,j]=add(atom(1,i,b) if j==a else {},atom(-1,j,b) if i==a else {})
            O[i,j]=add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},scale(mul(atom(1,b),der(P.get((i,j),{}),a)),-1))
        G.append(vec(V));G.append(vec(O))
    G.append(vec(P));G=transpose(G);rvecs=[vec(r) for r in normals]
    family=[row+[r[i] for r in rvecs] for i,row in enumerate(G)]
    assert rank(G)[0]==47 and rank(family)[0]==50
    assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in transpose(family))
    det=determinant([[D[i][j] for j in saved['implicit_coordinate_indices']] for i in saved['selected_equation_rows']])
    assert det and det==F(saved['implicit_minor_determinant'])
    assert len(saved['implicit_coordinate_indices'])==100
    C=decode(saved['full_curve_brackets'])
    expected=antisymmetric({ij:add(P.get(ij,{}),*(mul(atom(1,5+k),r.get(ij,{})) for k,r in enumerate(normals))) for ij in PAIRS})
    assert C==expected and all(not q for q in jacobi(C).values())
    # Independently reconstruct the recorded projective change, then compare
    # all ten E-wedge quotient relations with the immutable S61 reconstruction.
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
        if v:represented[ij][e+(0,0,0)]=v
    Delta=antisymmetric({ij:add(pushed[ij],scale(represented[ij],-1)) for ij in PAIRS})
    for i,j,k in combinations(range(5),3):
        assert not add(mul(atom(1,i),Delta[j,k]),scale(mul(atom(1,j),Delta[i,k]),-1),mul(atom(1,k),Delta[i,j]))
    assert len(reduced['equations'])==134 and all(sparse_eval(e)==0 for e in reduced['equations'])
    vals={'u'+str(i+1):v for i,v in enumerate(params)}
    res=[params[k]-evaluate(q,vals) for k,q in zip(toric['weight2_u_indices_zero_based'],toric['quadratic_section_q'])]
    assert len(res)==45 and res==[0]*45
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'method':'Independent standard-library Fraction reconstruction; no producer, SymPy or S61 arithmetic import',
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'full_D_shape':[350,150],'exact_D_rank':100,'orbit_vertical_scale_rank':47,'actual_three_parameter_family_rank':50,
        'nonzero_full_D_minor':str(det),'all150_times50_kernel_products_zero':True,'literal_three_parameter_family_full350_Jacobi_zero':True,
        'same_point_as_source_component_report':True,'PGL_frame_inverse_and_all_quotient_relations_reconstructed':True,
        'all134_original_reduced_equations_zero':True,'all45_T40_graph_equations_zero':True,
        'scope':'Witness is a smooth25 point of T40. Together with the separately proved complete-source one-orbit normalization, the entire q4 nil31-long I source lies in T40. The established universal T40 graph identification is an explicit dependency.',
        'seconds':time.monotonic()-started}
    (ROOT/'reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'full_D_rank':100,'family_rank':50,'T40_graph':True}))

if __name__=='__main__':main()
