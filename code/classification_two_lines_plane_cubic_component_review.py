#!/usr/bin/env python3
"""Independent Fraction reconstruction of the five complete two-line plane-cubic source components.
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
    names=dict(zip(('h','x','y','z','w','tau') if N==6 else ('h','x','y','z','w')+tuple('r'+str(i) for i in range(6)),range(N)));names.update({'x1':1,'x2':2,'x3':3,'x4':4})
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
    def ij(k):return ast.literal_eval(k) if k.startswith('(') else tuple(map(int,k.split(','))) if ',' in k else tuple(map(int,k))
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

def rawjac(P):
    return {ij:add(*(mul(P.get((a,l),{}),der(P.get((b,c),{}),l)) for a,b,c in [ij,(ij[1],ij[2],ij[0]),(ij[2],ij[0],ij[1])] for l in range(5))) for ij in combinations(range(5),3)}

def full_source_check(source):
    global N
    N=11;P=decode(source['full_six_magnetic_source'])
    assert all(not v for v in rawjac(P).values())
    assert all(not P.get((0,j),{}) for j in range(5))
    pf=add(mul(P.get((1,2),{}),P.get((3,4),{})),scale(mul(P.get((1,3),{}),P.get((2,4),{})),-1),mul(P.get((1,4),{}),P.get((2,3),{})))
    actualD=mul(atom(1,0),pf)
    assert actualD==parse(source['actual_degeneracy_quintic'])
    assert actualD==scale(parse(source['plane_quintic']),-F(source['integer_scale'])**2)
    values=[atom(1,j) for j in range(5)]+[atom(F(source['parameter_values']['r'+str(i)])) for i in range(6)]
    point={ij:substitute(poly,values) for ij,poly in P.items()}
    directions=[{ij:substitute(der(poly,5+j),values) for ij,poly in P.items()} for j in range(6)]
    reduced=lambda A:{ij:{e[:5]+(0,):v for e,v in poly.items()} for ij,poly in A.items() if poly}
    point=reduced(point);directions=list(map(reduced,directions))
    N=6
    assert point==decode(source['point_brackets'])
    assert directions==[decode(d) for d in source['all6_parameter_directions']]
    assert len(directions)==6 and rank(transpose([vec(d) for d in directions]))[0]==6
    return point,directions,{'all_six_parameter_strict_Jacobi_zero':True,'actual_Pfaffian_for_all_parameters_checked':True,
      'actual_divisor_independent_of_all_magnetic_parameters':True,'all_six_true_derivative_directions_reconstructed':True}

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic();paths=[ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json',ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_NORMALS_20261107.json']
    sources={r['case']:r for r in json.loads(paths[0].read_text())['cases']}
    normals={r['case']:r for r in json.loads(paths[1].read_text())['cases'] if r['zero_dimensional_quadratic_support']}
    expected={'cusp_flex_tangent','cusp_two_tangents_incident','cusp_flex_cusp_secant','cusp_tangent_cusp_secant','cusp_cusp_tangent_ordinary_tangent'}
    assert set(normals)==expected;records=[]
    for name,normal in normals.items():
        begin=time.monotonic();source=sources[name];P,directions,source_check=full_source_check(source)
        assert all(not v for v in jacobi(P).values());assert len(normal['normals'])==1
        R=decode(normal['normals'][0]);Dcols=[]
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
        V=transpose(VV);G=transpose(VV+OO+[vec(P)]);Fcols=VV+OO+[vec(P)]+[vec(A) for A in directions];Family=transpose(Fcols)
        rv=vec(R);K=[row+[v] for row,v in zip(Family,rv)]
        assert rank(V)[0]==24 and rank(G)[0]==47
        frank,fcols=rank(Family);assert frank==48 and rank(K)[0]==49
        assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in transpose(K))
        cert=source['certificates']['D'];rows=cert['rows'];cols=cert['columns'];assert len(rows)==len(cols)==101
        det=determinant([[D[i][j] for j in cols] for i in rows]);assert det and det==F(cert['determinant'])
        witness=normal['one_normal_quadratic_cokernel_functional'];left={int(i):F(c) for i,c in witness['nonzero_coefficients'].items()}
        assert all(sum(c*D[i][j] for i,c in left.items())==0 for j in range(150))
        Q=jvec(jacobi(R));pair=sum(c*Q[i] for i,c in left.items());assert pair and pair==F(witness['pairing'])
        record={'case':name,'all_passed':True,**source_check,'point_brackets':source['point_brackets'],
          'full_D_shape':[350,150],'full_D_exact_rank':101,'Euler_vertical_exact_rank':24,
          'actual_group_vertical_scale_exact_rank':47,'full_actual_source_direction_count':57,'full_actual_source_exact_rank':48,
          'full_kernel_spanned_by_actual_source_and_normal_rank':49,'all150_D_columns_annihilated_by_functional':True,
          'D_times_all_actual_source_and_normal_zero':True,'actual_source_basis_columns':fcols,
          'literal_normal':normal['normals'][0],'derivative_minor_rows':rows,'derivative_minor_columns':cols,
          'derivative_minor_determinant':str(det),'left_null_functional':witness['nonzero_coefficients'],'nonzero_quadratic_pairing':str(pair),
          'true_slice_dimension':102,'true_slice_basis':'The literal R plus101 standard coordinate vectors indexed by the derivative minor. Their span complements the48-dimensional actual source tangent: D annihilates source+R of rank49 and has the given invertible101minor.',
          'formal_reduced_transverse_support':'origin from the nonzero coefficient of n² in the true one-normal IFT residual',
          'projective_full_source_closure_component_dimension':23,
          'component_implication':'Cut by48 linear functions nondegenerate on the actual source tangent and vanishing on R. The true slice has101 invertible derivative equations and nonzero one-normal quadratic residual, hence dimension0. Krull bounds ambient local dimension by48, matching the actual source image lower bound. The irreducible source closure is therefore a component. Local uniqueness and dense orbit are not inferred.',
          'seconds':time.monotonic()-begin}
        records.append(record);print(json.dumps({'case':name,'all_passed':True,'seconds':record['seconds'],'pairing':str(pair)}),flush=True)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'method':'Independent standard-library Fraction sparse polynomials, full Jacobi derivative, all actual source directions, determinants, and cokernel pairings; no producer arithmetic or SymPy import.',
      'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':records,'seconds':time.monotonic()-started,
      'scope':'Independent arithmetic and dimension implication for five full23-dimensional components. Geometric exhaustion and pairwise/known-component distinctness are separate proof obligations. No dense open orbit or unique ambient component is claimed.'}
    out=ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_COMPONENT_REVIEW_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

if __name__=='__main__':main()
