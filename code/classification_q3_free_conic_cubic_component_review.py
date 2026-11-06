#!/usr/bin/env python3
"""Independent Fraction reconstruction of the three complete free conic-cubic source components.
No producer, producer engine, SymPy, NumPy or legacy arithmetic is imported.
"""
import ast
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
import hashlib,json,resource,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));QUADS=list(combinations(range(5),4))
N=6

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
    names=dict(zip(('a','b','c','r','z','tau') if N==6 else ('a','b','c','r','z')+tuple('m'+str(i) for i in range(6)),range(N)))
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

def decode(data):
    def ij(k):return ast.literal_eval(k) if k.startswith('(') else tuple(map(int,k))
    return antisymmetric({ij(k):parse(v) for k,v in data.items()})

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

def vec(P):return [P.get(ij,{}).get(e+(0,),F(0)) for ij,e in COLS]
def jvec(J,k=0):return [J[q].get(e+(k,),F(0)) for q,e in ROWS]
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

def substitute(p,values):
    out={}
    for e,c in p.items():
        term=atom(c)
        for i,k in enumerate(e):
            for _ in range(k):term=mul(term,values[i])
        out=add(out,term)
    return out

def full_source_check(src,canon):
    global N
    N=11
    P=decode(src['full_unimodular_brackets'])
    assert all(not v for v in jacobi(P).values())
    assert all(not add(*(der(P.get((i,j),{}),i) for i in range(5))) for j in range(5))
    base={ij:{e:c for e,c in poly.items() if not any(e[5:])} for ij,poly in P.items()}
    mu=parse(canon['magnetic_quotient_mu']);mon=parse(canon['canonical_complement_monomial'])
    source_mu={ij:dict(v) for ij,v in base.items()}
    source_mu[3,4]=add(source_mu[3,4],mul(mu,mon));source_mu[4,3]=scale(source_mu[3,4],-1)
    G=[[parse(v) for v in row] for row in canon['shear_G']]
    assert all(G[i][j]==(atom(1) if i==j else {}) for i in range(3) for j in range(5))
    assert all(G[i][j]==(atom(1) if i==j else {}) for i in [3,4] for j in [3,4])
    coords=[add(*(mul(G[i][j],atom(1,j)) for j in range(5))) for i in range(5)]
    values=coords+[atom(1,j) for j in range(5,11)]
    for i,j in PAIRS:
        actual=add(*(mul(mul(G[i][a],source_mu.get((a,b),{})),G[j][b]) for a,b in product(range(5),repeat=2)))
        assert actual==substitute(P.get((i,j),{}),values)
    U=[base[i,3] for i in range(3)];V=[base[i,4] for i in range(3)]
    du=add(*(der(U[i],i) for i in range(3)));dv=add(*(der(V[i],i) for i in range(3)))
    images=[add(V[i],mul(dv,atom(1,i))) for i in range(3)]+[scale(add(U[i],mul(du,atom(1,i))),-1) for i in range(3)]
    monomial_indices=[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
    monomials11=[next(iter(atom(1,*ij))) for ij in monomial_indices]
    M=[[p.get(e,F(0)) for p in images] for e in monomials11]
    assert M==[[F(v) for v in row] for row in canon['magnetic_translation_matrix']]
    assert rank(M)[0]==5
    complement=[mon.get(e,F(0)) for e in monomials11]
    assert rank([row+[v] for row,v in zip(M,complement)])[0]==6
    left=[mu.get(next(iter(atom(1,j))),F(0)) for j in range(5,11)]
    assert all(sum(a*b for a,b in zip(left,col))==0 for col in transpose(M))
    assert sum(a*b for a,b in zip(left,complement))==1
    canonical={ij:dict(v) for ij,v in base.items()}
    canonical[3,4]=add(canonical[3,4],mon);canonical[4,3]=scale(canonical[3,4],-1)
    canonical={ij:scale(v,F(canon['canonical_integer_scale'])) for ij,v in canonical.items()}
    for ij,poly in base.items():
        if ij[0]>ij[1]:continue
        # Mixed entries have no fiber variables; the vertical entry is fiber-linear.
        assert all(sum(e[3:5])==(1 if ij==(3,4) else 0) for e in poly)
    reduced={ij:{e[:5]+(0,):v for e,v in poly.items()} for ij,poly in canonical.items()}
    N=6
    assert {ij:v for ij,v in reduced.items() if v}==decode(canon['canonical_point_brackets'])
    return {'full_six_parameter_Jacobi_and_divergence_zero':True,'whole_source_PGL_shear_identity':True,
      'magnetic_translation_rank':5,'complement_augmented_rank':6,'mu_quotient_checked':True,
      'canonical_scale_identity_checked':True,'fiber_scaling_and_mu_zero_closure_by_verified_monomial_support':True}

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    paths={kind:ROOT/('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_'+kind+'_20261107.json') for kind in ['SOURCES','CANONICAL','NORMALS']}
    data={kind:{row['case']:row for row in json.loads(p.read_text())['cases']} for kind,p in paths.items()}
    records=[]
    for name in ['cusp_outside_1','node_smooth_contact','node_at_node']:
        begin=time.monotonic();canon=data['CANONICAL'][name];normal=data['NORMALS'][name]
        source_check=full_source_check(data['SOURCES'][name],canon)
        P=decode(canon['canonical_point_brackets']);R=decode(normal['literal_normals'][0])
        assert all(not v for v in jacobi(P).values())
        Dcols=[]
        for ij,e in COLS:
            v={key:dict(value) for key,value in P.items()}
            inc={e+(0,):F(1)};v[ij]=add(v.get(ij,{}),inc);v[ij[::-1]]=scale(v[ij],-1)
            Dcols.append(jvec(jacobi(v)))
        D=transpose(Dcols);VV=[];OO=[]
        for a,b in product(range(5),repeat=2):
            V={};O={}
            for i,j in PAIRS:
                V[i,j]=add(atom(1,i,b) if j==a else {},atom(-1,j,b) if i==a else {})
                O[i,j]=add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},scale(mul(atom(1,b),der(P.get((i,j),{}),a)),-1))
            VV.append(vec(V));OO.append(vec(O))
        V=transpose(VV);G=transpose(VV+OO+[vec(P)]);rv=vec(R);K=[row+[v] for row,v in zip(G,rv)]
        assert rank(V)[0]==24
        grank,gcols=rank(G);assert grank==48 and rank(K)[0]==49
        assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in transpose(K))
        rows=normal['implicit_D_minor_rows'];cols=normal['implicit_D_minor_columns']
        assert len(rows)==len(cols)==101
        det=determinant([[D[i][j] for j in cols] for i in rows])
        assert det and det==F(normal['implicit_D_minor_determinant'])
        witness=normal['left_null_witnesses'][0];left={int(i):F(c) for i,c in witness['functional'].items()}
        assert all(sum(c*D[i][j] for i,c in left.items())==0 for j in range(150))
        Q=jvec(jacobi(R));pair=sum(c*Q[i] for i,c in left.items())
        assert pair and pair==F(witness['pairings'][0])
        record={'case':name,'all_passed':True,**source_check,'point_brackets':canon['canonical_point_brackets'],
          'full_D_shape':[350,150],'full_D_exact_rank':101,'Euler_vertical_exact_rank':24,
          'actual_orbit_vertical_scale_exact_rank':48,'full_kernel_spanned_by_group_and_normal_rank':49,
          'all150_D_columns_annihilated_by_functional':True,'D_times_all_actual_group_and_normal_zero':True,
          'actual_group_basis_columns':gcols,'literal_normal':normal['literal_normals'][0],
          'derivative_minor_rows':rows,'derivative_minor_columns':cols,'derivative_minor_determinant':str(det),
          'left_null_functional':witness['functional'],'nonzero_quadratic_pairing':str(pair),
          'true_slice_dimension':102,'true_slice_basis':'the literal normal R and the101 standard coordinate vectors listed in derivative_minor_columns; this complements a48-dimensional actual group tangent because D kills group+R, which has rank49, and the101-column minor is invertible',
          'formal_reduced_transverse_support':'origin; the cokernel pairing gives nonzero coefficient of n² in a true one-normal IFT slice',
          'projective_orbit_closure_component_dimension':23,'component_implication':'True102slice plus101 IFT equations leaves one normal; its nonzero quadratic initial form forces zero-dimensional support. Krull gives local dimension at most48 before quotient, matching actual48-dimensional orbit. Quotient vertical24 and scale1 gives component23.',
          'seconds':time.monotonic()-begin}
        records.append(record);print(json.dumps({'case':name,'all_passed':True,'seconds':record['seconds'],'pairing':str(pair)}),flush=True)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'method':'Independent standard-library Fraction polynomial, Jacobi, derivative, orbit, determinant, and cokernel reconstruction; neither producer nor SymPy imported.',
      'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},
      'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':records,'seconds':time.monotonic()-started,
      'scope':'Independent arithmetic for three whole-source component23 claims. Finite geometric configuration exhaustion and component distinctness are separate proof obligations.'}
    out=ROOT/'reports/CLASSIFICATION_Q3_FREE_CONIC_CUBIC_COMPONENT_REVIEW_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

if __name__=='__main__':main()
