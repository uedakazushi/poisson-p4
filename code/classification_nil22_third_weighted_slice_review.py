#!/usr/bin/env python3
"""Independent Fraction reconstruction of the full two-normal weighted slice certificate.
No producer, producer engine, SymPy, NumPy or legacy arithmetic is imported.
"""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
import hashlib,json,resource,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));QUADS=list(combinations(range(5),4))
N=7

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
    names=dict(zip(('h','x1','x2','x3','x4','tau0','tau1'),range(N)))
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

def vec(P):return [P.get(ij,{}).get(e+(0,0),F(0)) for ij,e in COLS]
def jvec(J,ex=(0,0)):return [J[q].get(e+ex,F(0)) for q,e in ROWS]
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

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    path=ROOT/'reports/CLASSIFICATION_NIL22_THIRD_WEIGHTED_SLICE_20261107.json'
    saved=json.loads(path.read_text());P=decode(saved['canonical_brackets']);normals=[decode(r) for r in saved['normal_brackets']]
    assert all(not q for q in jacobi(P).values())
    Dcols=[]
    for ij,e in COLS:
        v={key:dict(value) for key,value in P.items()}
        inc={e+(0,0):F(1)};v[ij]=add(v.get(ij,{}),inc);v[ij[::-1]]=scale(v[ij],-1)
        Dcols.append(jvec(jacobi(v)))
    D=transpose(Dcols)
    G=[]
    for a,b in product(range(5),repeat=2):
        V={};O={}
        for i,j in PAIRS:
            V[i,j]=add(atom(1,i,b) if j==a else {},atom(-1,j,b) if i==a else {})
            O[i,j]=add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},scale(mul(atom(1,b),der(P.get((i,j),{}),a)),-1))
        G.append(vec(V));G.append(vec(O))
    G.append(vec(P));G=transpose(G)
    rvecs=[vec(r) for r in normals];kernel=[row+[r[i] for r in rvecs] for i,row in enumerate(G)]
    assert all(sum(a*b for a,b in zip(drow,kcol))==0 for drow in D for kcol in transpose(kernel))
    assert rank(G)[0]==47 and rank(kernel)[0]==49
    piv=saved['actual_group_tangent_pivot_rows'];sl=saved['stable_slice_coordinate_indices']
    assert sorted(piv+sl)==list(range(150)) and len(sl)==103 and rank([G[i] for i in piv])[0]==47
    implicit=saved['implicit_coordinate_indices'];eq=saved['selected_equation_rows']
    assert len(implicit)==len(eq)==101 and set(implicit)<set(sl)
    det=determinant([[D[i][j] for j in implicit] for i in eq])
    assert det==F(saved['implicit_minor_determinant']) and det
    free=sorted(set(sl)-set(implicit));assert len(free)==2
    assert all(rvecs[a][i]==0 for a in range(2) for i in piv)
    assert [[rvecs[a][i] for i in free] for a in range(2)]==[[F(1),F(0)],[F(0),F(1)]]
    weights=saved['coordinate_weights']
    relative=[weights[i]+weights[j]-sum(a*b for a,b in zip(weights,e))-1 for (i,j),e in COLS]
    assert (min(relative),max(relative))==(-6,4)
    assert all({relative[i] for i,v in enumerate(r) if v}=={-1} for r in rvecs)
    assert {relative[i] for i,v in enumerate(vec(P)) if v}=={0}
    assert all(len({relative[i] for i,v in enumerate(col) if v})<=1 for col in transpose(G))
    curve=decode(saved['full_curve_brackets'])
    for at,((i,j),e) in enumerate(COLS):
        poly=curve.get((i,j),{})
        assert poly.get(e+(0,0),F(0))==vec(P)[at]
        assert poly.get(e+(1,0),F(0))==rvecs[0][at]
        assert poly.get(e+(0,1),F(0))==rvecs[1][at]
        for ex,v in poly.items():
            if ex[:5]!=e or sum(ex[5:])<2:continue
            assert at in implicit and relative[at]==-sum(ex[5:]) and sum(ex[5:])<=6
    J=jacobi(curve)
    residual=[]
    for q,e in ROWS:
        residual.append({(0,0,0,0,0)+ex[5:]:c for ex,c in J[q].items() if ex[:5]==e})
    assert all(not residual[i] for i in eq)
    assert all(residual[i]==parse(saved['full350_residual_polynomials'].get(str(i),'0')) for i in range(350))
    q1=residual[115];q2=residual[157];q4=residual[98]
    assert q1==parse('-5*tau0*(35*tau0+36*tau1)/1056')
    assert q2==parse('5*tau0*(7*tau0+12*tau1)/48')
    assert all(sum(e[5:])==4 for e in q4)
    a1=q1.get((0,0,0,0,0,2,0),F(0));b1=q1.get((0,0,0,0,0,1,1),F(0))
    a2=q2.get((0,0,0,0,0,2,0),F(0));b2=q2.get((0,0,0,0,0,1,1),F(0))
    coeffdet=a1*b2-a2*b1;assert coeffdet
    restricted={e:c for e,c in q4.items() if e[5]==0}
    assert restricted=={(0,0,0,0,0,0,4):F(-5*591951360,204578390016)}
    assert all(sum(e[5:])>=2 for p in residual for e in p)
    result={'task_id':'E_multiple_hyperplane','all_passed':True,
        'method':'Independent standard-library Fraction sparse reconstruction; no producer or producer-engine arithmetic import',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'full_D_shape':[350,150],'full_D_exact_rank':101,'orbit_vertical_scale_exact_rank':47,
        'full_kernel_basis_rank':49,'all_actual_direction_and_normal_products_zero':True,
        'stable_slice_dimension':103,'IFT_dimension':101,'IFT_minor_determinant':str(det),
        'normal_relative_weights':[-1,-1],'implicit_weight_bound':6,
        'full350_residuals_reconstructed':True,'all101_selected_equations_identically_zero':True,
        'all_saved_corrections_in_true_slice_and_correct_weight':True,
        'quadratic_rows':[115,157],'quadratic_coefficient_determinant':str(coeffdet),
        'quartic_row':98,'quartic_tau0_zero_coefficient':str(next(iter(restricted.values()))),
        'reduced_formal_support':'the origin in both true normal coordinates',
        'scope':'The two independent quadratic rows give tau0^2 and tau0*tau1, hence tau0 belongs to the radical; row98 modulo tau0 is a nonzero multiple of tau1^4. The origin is in the support. Exact nilpotence thickness not claimed.',
        'seconds':time.monotonic()-started}
    out=ROOT/'reports/CLASSIFICATION_NIL22_THIRD_WEIGHTED_SLICE_REVIEW_20261107.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'quadratic_det':str(coeffdet),'quartic':str(next(iter(restricted.values())))}))

if __name__=='__main__':main()
