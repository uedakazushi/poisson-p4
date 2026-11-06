#!/usr/bin/env python3
"""Independent Fraction review of all31 degenerate-quadric cofactor sources.
Arithmetic primitives reused from the dependent-source independent reviewer;
no producer, SymPy, or producer Jacobi/rank helper is imported.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import ast,hashlib,json,resource,sys,time

resource.setrlimit(resource.RLIMIT_CPU,(115,116))
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
started=time.monotonic();ROOT=Path.cwd();N=28
pairs=list(combinations(range(5),2))
def term(a,*inds):
    m=[0]*N
    for i in inds:m[i]+=1
    return {tuple(m):F(a)} if a else {}
def add(*ps):
    r={}
    for p in ps:
        for m,a in p.items():r[m]=r.get(m,F(0))+a
    return {m:a for m,a in r.items() if a}
def scale(p,a):return {m:a*b for m,b in p.items() if a*b}
def mul(p,q):
    r={}
    for m,a in p.items():
        for n,b in q.items():
            k=tuple(x+y for x,y in zip(m,n));r[k]=r.get(k,F(0))+a*b
    return {m:a for m,a in r.items() if a}
def der(p,i):
    r={}
    for m,a in p.items():
        if m[i]:
            k=list(m);k[i]-=1;r[tuple(k)]=a*m[i]
    return r
def tensor(data):
    r={}
    for (i,j),p in data.items():r[i,j]=p;r[j,i]=scale(p,-1)
    return r
def projective_jacobi(A):
    triples={}
    for i,j,k in combinations(range(5),3):
        triples[i,j,k]=add(*(mul(A.get((a,l),{}),der(A.get((b,c),{}),l))
                for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5)))
    return {inds:add(*(scale(mul(term(1,inds[k]),triples[inds[:k]+inds[k+1:]]),(-1)**k)
              for k in range(4))) for inds in combinations(range(5),4)}
def monomials(degree):
    out=[]
    for factors in combinations_with_replacement(range(5),degree):
        m=[0]*5
        for i in factors:m[i]+=1
        out.append(tuple(m))
    return out
mon2=monomials(2);mon4=monomials(4)
column_basis=[(ij,m) for ij in pairs for m in mon2]
row_basis=[(inds,m) for inds in combinations(range(5),4) for m in mon4]
row_index={v:i for i,v in enumerate(row_basis)}
def bivector_vector(A):return [A.get(ij,{}).get(m+(0,)*(N-5),F(0)) for ij,m in column_basis]
def jacobi_vector(A):
    J=projective_jacobi(A)
    return [J[inds].get(m+(0,)*(N-5),F(0)) for inds,m in row_basis]
def parse_polynomial(text):
    names=dict(zip(tuple('x'+str(i) for i in range(5))+tuple('n'+str(i) for i in range(22)),range(27)))
    def walk(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int):return term(n.value)
        if isinstance(n,ast.Name) and n.id in names:return term(1,names[n.id])
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return scale(walk(n.operand),-1)
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.UAdd):return walk(n.operand)
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Add):return add(walk(n.left),walk(n.right))
            if isinstance(n.op,ast.Sub):return add(walk(n.left),scale(walk(n.right),-1))
            if isinstance(n.op,ast.Mult):return mul(walk(n.left),walk(n.right))
            if isinstance(n.op,ast.Pow) and isinstance(n.right,ast.Constant) and isinstance(n.right.value,int) and n.right.value>=0:
                r=term(1)
                for _ in range(n.right.value):r=mul(r,walk(n.left))
                return r
            if isinstance(n.op,ast.Div) and isinstance(n.right,ast.Constant) and isinstance(n.right.value,int):return scale(walk(n.left),F(1,n.right.value))
        raise ValueError(ast.dump(n))
    return walk(ast.parse(text,mode='eval').body)
def parse_tensor(data):
    return tensor({tuple(ast.literal_eval(ij)):parse_polynomial(p) for ij,p in data.items()})

def modular_minor(A,prime=1000003):
    original=[[int(a.numerator)*pow(int(a.denominator),-1,prime)%prime if isinstance(a,F) else a%prime
               for a in row] for row in A]
    B=[r[:] for r in original];nr=len(B);nc=len(B[0]);labels=list(range(nr));rows=[];cols=[];rank=0
    for j in range(nc):
        k=next((i for i in range(rank,nr) if B[i][j]),None)
        if k is None:continue
        B[rank],B[k]=B[k],B[rank];labels[rank],labels[k]=labels[k],labels[rank]
        rows.append(labels[rank]);cols.append(j)
        inv=pow(B[rank][j],-1,prime);B[rank]=[a*inv%prime for a in B[rank]]
        for i in range(rank+1,nr):
            a=B[i][j]
            if a:B[i]=[(x-a*y)%prime for x,y in zip(B[i],B[rank])]
        rank+=1
    minor=[[original[i][j] for j in cols] for i in rows];det=1
    for j in range(rank):
        k=next(i for i in range(j,rank) if minor[i][j])
        if k!=j:minor[k],minor[j]=minor[j],minor[k];det=-det
        a=minor[j][j];det=det*a%prime;inv=pow(a,-1,prime)
        for i in range(j+1,rank):
            b=minor[i][j]*inv%prime
            if b:minor[i]=[(x-b*y)%prime for x,y in zip(minor[i],minor[j])]
    assert det
    return {'prime':prime,'rank_mod_prime':rank,'minor_rows':rows,'minor_columns':cols,'minor_determinant_mod_prime':det}

def strict_jacobi(A):
    return [add(*(mul(A.get((a,l),{}),der(A.get((b,c),{}),l))
            for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5)))
            for i,j,k in combinations(range(5),3)]

def specialize(p, values):
    out={}
    for m,a in p.items():
        mm=list(m)
        for i,value in values.items():a*=value**mm[i];mm[i]=0
        out[tuple(mm)]=out.get(tuple(mm),F(0))+a
    return {m:a for m,a in out.items() if a}

def homogeneous_linear(p):
    out={}
    for m,a in p.items():
        assert sum(m)<=1
        mm=list(m);mm[0]+=1-sum(m);out[tuple(mm)]=a
    return out

def rational_relations(A):
    B=[[F(a) for a in row] for row in A];nr=len(B);nc=len(B[0]);pivots=[];rank=0
    for j in range(nc):
        k=next((i for i in range(rank,nr) if B[i][j]),None)
        if k is None:continue
        B[rank],B[k]=B[k],B[rank];a=B[rank][j];B[rank]=[v/a for v in B[rank]]
        for i in range(nr):
            a=B[i][j]
            if i!=rank and a:B[i]=[x-a*y for x,y in zip(B[i],B[rank])]
        pivots.append(j);rank+=1
    free=[j for j in range(nc) if j not in pivots];relations=[]
    for j in free:
        rel=[F(0)]*nc;rel[j]=1
        for i,k in enumerate(pivots):rel[k]=-B[i][j]
        assert all(sum(x*y for x,y in zip(row,rel))==0 for row in A)
        relations.append(rel)
    return {'rank':rank,'free_columns':free,'independent_kernel_relations':[[str(a) for a in rel] for rel in relations]}

def matrix_hash(A):
    # Match the frozen producer's SymPy-string matrix serialization exactly.
    return hashlib.sha256(json.dumps([[str(a) for a in row] for row in A]).encode()).hexdigest()

def compose(p,substitutions):
    out={}
    for m,a in p.items():
        value=term(a)
        for i,power in enumerate(m):
            factor=substitutions.get(i,term(1,i))
            for _ in range(power):value=mul(value,factor)
        out=add(out,value)
    return out

def homogenize(p,degree):
    out={}
    for m,a in p.items():
        assert m[0]==0 and sum(m[:5])<=degree
        mm=list(m);mm[0]=degree-sum(m[:5]);out[tuple(mm)]=a
    return out

def decode(data):return tensor({tuple(map(int,ij.split(','))):parse_polynomial(p) for ij,p in data.items()})

def pfaffian(A):
    return add(mul(A.get((1,2),{}),A.get((3,4),{})),scale(mul(A.get((1,3),{}),A.get((2,4),{})),-1),mul(A.get((1,4),{}),A.get((2,3),{})))

def five_pfaffians(A):
    out=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5) if a!=missing]
        out.append(scale(add(mul(A.get((i,j),{}),A.get((k,l),{})),
            scale(mul(A.get((i,k),{}),A.get((j,l),{})),-1),mul(A.get((i,l),{}),A.get((j,k),{}))),(-1)**missing))
    return out

def divide_exact(p,d):
    remainder=dict(p);quotient={};lead=max(d);leading=d[lead]
    while remainder:
        top=max(remainder);diff=tuple(a-b for a,b in zip(top,lead))
        assert min(diff)>=0
        t={diff:remainder[top]/leading};quotient=add(quotient,t)
        remainder=add(remainder,scale(mul(t,d),-1))
    assert mul(quotient,d)==p
    return quotient

def jacobi_coefficient_witnesses(P,targets):
    raw=strict_jacobi(P);slots=[]
    for triple,p in zip(combinations(range(5),3),raw):
        for xm in sorted(set(m[:5] for m in p)):
            coefficient={(0,)*5+m[5:]:a for m,a in p.items() if m[:5]==xm}
            slots.append((triple,xm,coefficient))
    witnesses=[]
    for text in targets:
        target=parse_polynomial(text);key=next(iter(target));found=None
        for triple,xm,p in slots:
            if key in p:
                ratio=p[key]/target[key]
                if p==scale(target,ratio):
                    found={'condition':text,'triple':list(triple),'coordinate_monomial':list(xm),'coefficient_ratio':str(ratio)};break
        assert found,text
        witnesses.append(found)
    return raw,witnesses

table_path=Path('code/data/classification_degenerate_quadric_cofactors.json')
table=json.loads(table_path.read_text());records=[];input_hashes={str(table_path):hashlib.sha256(table_path.read_bytes()).hexdigest()}
for spec in table['cases']:
    begun=time.monotonic();kind=spec['name']
    path=Path(f'reports/CLASSIFICATION_DEGENERATE_QUADRIC_{kind.upper()}_20261107.json')
    saved=json.loads(path.read_text());count=saved['source_dimension']
    input_hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    assert saved['fixed_specification']==spec
    assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['producer_sha256']
    q=parse_polynomial(spec['quadric']);V=[parse_polynomial(f) for f in spec['cofactor_components']]
    M=[[der(V[i],j).get((0,)*N,F(0)) for j in range(5)] for i in range(5)]
    assert sum(M[i][i] for i in range(5))==0
    assert not add(*(mul(V[i],der(q,i)) for i in range(5)))
    H=[[der(der(q,i),j).get((0,)*N,F(0)) for j in range(5)] for i in range(5)]
    qcert=modular_minor(H);assert qcert['rank_mod_prime']==spec['quadric_rank']
    P=tensor({tuple(map(int,key)):parse_polynomial(value) for key,value in saved['brackets'].items()})
    def linear_map(A):
        div=[add(*(der(A.get((i,j),{}),i) for i in range(5))) for j in range(5)]
        normal=[add(*(mul(der(q,i),A.get((i,j),{})) for i in range(5))) for j in range(5)]
        lie=[add(*(mul(V[k],der(A.get((i,j),{}),k)) for k in range(5)),
            *(scale(A.get((k,j),{}),-M[i][k]) for k in range(5)),
            *(scale(A.get((i,k),{}),-M[j][k]) for k in range(5))) for i,j in pairs]
        return div+normal+lie
    rhs=[{} for _ in range(5)]+[mul(q,v) for v in V]+[{} for _ in range(10)]
    assert linear_map(P)==rhs
    rows=[(j,m+(0,)*(N-5)) for j in range(20) for m in monomials(1 if j<5 else 3 if j<10 else 2)]
    columns=[]
    for ij,m in column_basis:
        vals=linear_map(tensor({ij:{m+(0,)*(N-5):F(1)}}))
        columns.append([vals[j].get(e,F(0)) for j,e in rows])
    full=list(map(list,zip(*columns)));cert=modular_minor(full)
    directions=[tensor({ij:der(P.get(ij,{}),5+j) for ij in pairs}) for j in range(count)]
    assert all(not p for A in directions for p in linear_map(A))
    dir_matrix=list(map(list,zip(*(bivector_vector(A) for A in directions))))
    dircert=modular_minor(dir_matrix)
    assert dircert['rank_mod_prime']==count and cert['rank_mod_prime']==150-count==saved['linear_rank']
    pfs=five_pfaffians(P)
    potential=scale(add(*(mul(term(1,i),pfs[i]) for i in range(5))),F(1,5))
    assert potential==parse_polynomial(saved['euler_contracted_pfaffian_not_yet_actual_divisor'])
    coefficients=set()
    for poly in strict_jacobi(P):
        for xm in set(m[:5] for m in poly):
            coefficient={(0,)*5+m[5:]:v for m,v in poly.items() if m[:5]==xm}
            coefficients.add(tuple(sorted(coefficient.items())))
    reported={tuple(sorted(parse_polynomial(p).items())) for p in saved['residual_jacobi_coefficients']}
    assert coefficients==reported,(kind,len(coefficients),len(reported))
    record={'case':kind,'quadric_rank':spec['quadric_rank'],'full_linear_matrix_shape':[len(full),150],
        'full_linear_rank_QQ_by_kernel_sandwich':150-count,'source_dimension':count,
        'source_all_linear_constraints_exact':True,'all_source_directions_independent':True,
        'linear_rank_lower_certificate':cert,'source_rank_lower_certificate':dircert,
        'Euler_contracted_Pfaffian_reproduced_before_Jacobi':True,
        'all_ten_Jacobi_polynomials_and_every_nonzero_coefficient_reproduced':True,
        'distinct_Jacobi_coefficients':len(coefficients),'seconds':time.monotonic()-begun}
    records.append(record)
    print(json.dumps({k:record[k] for k in ['case','source_dimension','distinct_Jacobi_coefficients','seconds']}),flush=True)
checker=Path(__file__).read_text()
result={'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ','method':'Independent Fraction sparse arithmetic, full350-slot linear constraints and full150 coefficient sources; no producer or SymPy import.',
    'input_sha256':input_hashes,'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
    'cases':records,'elapsed_seconds':time.monotonic()-started,
    'scope':'All31 full fixed-cofactor sources independently reconstructed. All residual Jacobi and Euler Pfaffian retained; this alone is not nonlinear support coverage.',
    'global_exhaustion_proved':False}
with Path(sys.argv[1]).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'all_passed':True,'cases':len(records),'elapsed_seconds':result['elapsed_seconds']},indent=2))

