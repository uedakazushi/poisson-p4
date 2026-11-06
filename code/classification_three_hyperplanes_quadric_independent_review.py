#!/usr/bin/env python3
"""Independent Fraction review of four full conic/triangle magnetic sources.
Arithmetic primitives reused from the dependent-source independent reviewer;
no producer, SymPy, or producer Jacobi/rank helper is imported.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import ast,hashlib,json,resource,sys,time

resource.setrlimit(resource.RLIMIT_CPU,(55,56))
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
started=time.monotonic();ROOT=Path.cwd();N=12
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
    names=dict(zip(('h','x','y','z','w')+tuple('r'+str(i) for i in range(6)),range(11)))
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

source_path='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_20261107.json'
normal_path='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_NORMALS_20261107.json'
source=json.loads((ROOT/source_path).read_text());saved_normals=json.loads((ROOT/normal_path).read_text())
for saved in source['source_snapshots']:
    assert hashlib.sha256(saved['source'].encode()).hexdigest()==saved['sha256']
    assert source['input_sha256'][saved['path']]==saved['sha256']
normal_cases={r['case']:r for r in saved_normals['cases']}
base_mons=[term(1,*inds) for inds in [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]]
records=[]
for case_index,saved in enumerate(source['cases']):
    case_start=time.monotonic();name=saved['case'];common=F(saved['source_integer_scale'])
    q=parse_polynomial(saved['quadric']);plane=mul(mul(term(1,0,1),term(1,2)),q)
    assert plane==parse_polynomial(saved['plane_quintic'])
    gradient=[der(plane,i) for i in range(3)]
    syzygies=[[parse_polynomial(p) for p in row] for row in saved['two_full_syzygy_vectors']]
    for row in syzygies:assert not add(*(mul(a,b) for a,b in zip(row,gradient)))
    syzcols=[mul(m,gradient[i]) for i in range(3) for m in base_mons]
    syzrows=sorted(set().union(*(p.keys() for p in syzcols)))
    syzmatrix=[[p.get(m,F(0)) for p in syzcols] for m in syzrows]
    syzcert=modular_minor(syzmatrix);assert syzcert['rank_mod_prime']==16
    syzvecs=[[p.get(next(iter(m)),F(0)) for p in row for m in base_mons] for row in syzygies]
    assert modular_minor(list(map(list,zip(*syzvecs))))['rank_mod_prime']==2
    # Project the two full homogeneous syzygies by subtracting their Euler components.
    projected=[]
    for row in syzygies:
        hh={}
        for m,a in row[0].items():
            assert m[0];mm=list(m);mm[0]-=1;hh[tuple(mm)]=a
        projected.append([add(row[i],scale(mul(term(1,i),hh),-1)) for i in [1,2]])
    U=list(map(parse_polynomial,saved['anchor_U']));V=list(map(parse_polynomial,saved['anchor_V']))
    assert U==projected[0]
    ratio=None
    for p,v in zip(projected[1],V):
        assert p.keys()==v.keys()
        for m,a in p.items():
            if ratio is None:ratio=v[m]/a
            assert v[m]==ratio*a
    assert ratio
    assert add(mul(U[0],V[1]),scale(mul(U[1],V[0]),-1))==mul(term(1,1,2),q)
    ch,dh=map(lambda p:homogeneous_linear(parse_polynomial(p)),saved['affine_anchor_bracket_coefficients'])
    bracket=[add(*(add(mul(U[j],der(V[i],j+1)),scale(mul(V[j],der(U[i],j+1)),-1)) for j in range(2))) for i in range(2)]
    assert bracket==[add(mul(ch,U[i]),mul(dh,V[i])) for i in range(2)]
    magnetic=add(*(mul(term(1,5+i),m) for i,m in enumerate(base_mons)))
    Pgen=tensor({(1,3):scale(U[0],common),(2,3):scale(U[1],common),
        (1,4):scale(V[0],common),(2,4):scale(V[1],common),
        (3,4):scale(add(scale(mul(ch,term(1,3)),-1),scale(mul(dh,term(1,4)),-1),magnetic),common)})
    assert Pgen==parse_tensor(saved['all_six_parameter_source'])
    assert all(not p for p in strict_jacobi(Pgen))
    pf=add(mul(Pgen.get((1,2),{}),Pgen.get((3,4),{})),scale(mul(Pgen.get((1,3),{}),Pgen.get((2,4),{})),-1),mul(Pgen.get((1,4),{}),Pgen.get((2,3),{})))
    assert pf==scale(mul(term(1,1,2),q),-common**2)
    values={5+i:F(saved['point_substitution']['r'+str(i)]) for i in range(6)}
    P=tensor({ij:specialize(Pgen.get(ij,{}),values) for ij in pairs})
    # Omit literal zero entries in tensor comparison.
    P={ij:p for ij,p in P.items() if p}
    assert P==parse_tensor(saved['point_brackets'])
    parameters=[parse_tensor(d) for d in saved['all6_actual_parameter_directions']]
    for i,R in enumerate(parameters):assert R==tensor({(3,4):scale(base_mons[i],common)})
    D=[[0]*150 for _ in row_basis]
    for col,(ij,m) in enumerate(column_basis):
        variation=tensor({ij:{m+(0,)*(N-5):F(1)}})
        mixed=tensor({kl:add(P.get(kl,{}),mul(term(1,11),variation.get(kl,{}))) for kl in pairs})
        for inds,p in projective_jacobi(mixed).items():
            for exp,a in p.items():
                if exp[11]==1:
                    assert a.denominator==1 and not any(exp[5:11])
                    D[row_index[inds,exp[:5]]][col]=int(a)
    vectors=[]
    for i,j in product(range(5),repeat=2):
        A=tensor({(k,l):add(term(int(l==i),k,j),term(-int(k==i),l,j)) for k,l in pairs})
        vectors.append(bivector_vector(A))
    for i,j in product(range(5),repeat=2):
        A=tensor({(k,l):add(scale(mul(term(1,j),der(P.get((k,l),{}),i)),-1),
            scale(P.get((j,l),{}),int(k==i)),scale(P.get((k,j),{}),int(l==i))) for k,l in pairs})
        vectors.append(bivector_vector(A))
    vectors.append(bivector_vector(P));vectors.extend(map(bivector_vector,parameters))
    source_matrix=list(map(list,zip(*vectors)))
    matrices={'full350_D':D,'vertical':[r[:25] for r in source_matrix],
        'orbit_vertical_scale':[r[:51] for r in source_matrix],'actual_full_source':source_matrix}
    for key,A in matrices.items():assert matrix_hash(A)==saved['matrix_sha256'][key]
    dcert=modular_minor(D)
    ranks={};upper_relations={};lower_certificates={}
    for key,A in matrices.items():
        if key=='full350_D':continue
        cert=modular_minor(A);upper=rational_relations(A)
        assert cert['rank_mod_prime']==upper['rank']==saved['ranks_QQ'][key]
        ranks[key]=upper['rank'];upper_relations[key]=upper;lower_certificates[key]=cert
    assert ranks['actual_full_source']==48
    normals=[parse_tensor(d) for d in normal_cases.get(name,{}).get('normals',[])]
    vectors.extend(map(bivector_vector,normals));kernel_matrix=list(map(list,zip(*vectors)))
    for row in D:
        nz=[(j,a) for j,a in enumerate(row) if a]
        assert all(sum(a*v[j] for j,a in nz)==0 for v in vectors)
    kcert=modular_minor(kernel_matrix)
    assert dcert['rank_mod_prime']==102-case_index==saved['ranks_QQ']['full350_D']
    assert dcert['rank_mod_prime']+kcert['rank_mod_prime']==150
    record={'case':name,'universal_all_six_parameter_source_reconstructed':True,
        'all_ten_strict_Jacobi_identically_zero':True,'universal_Pfaffian_exact':True,
        'full_gradient_syzygy_rank_QQ':16,'full_gradient_syzygy_dimension_QQ':2,
        'syzygy_rank_lower_certificate':syzcert,'all_four_frozen_matrix_hashes_reproduced':True,
        'full_D_rank_QQ_by_sandwich':dcert['rank_mod_prime'],'actual_source_rank_QQ':48,
        'full_kernel_dimension_QQ':kcert['rank_mod_prime'],'true_normal_dimension':len(normals),
        'D_times_all_actual_source_and_saved_normal_directions_exact_zero':True,
        'D_rank_lower_certificate':dcert,'kernel_rank_lower_certificate':kcert,
        'source_rank_lower_certificates':lower_certificates,'source_exact_upper_certificates':upper_relations,
        'projective_source_dimension_lower':23,'projective_tangent_dimension_upper':23+case_index,
        'smooth_component_dimension_if_irreducible_source_argument_accepted':23 if case_index==0 else None,
        'quadratic_obstruction_arithmetic_checked_here':False,'seconds':time.monotonic()-case_start}
    records.append(record)
    print(json.dumps({k:record[k] for k in ['case','full_D_rank_QQ_by_sandwich','actual_source_rank_QQ','true_normal_dimension','seconds']}),flush=True)

checker=Path(__file__).read_text()
result={'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ',
    'method':'Standalone Fraction sparse polynomial arithmetic, AST literal parsing and exact rational elimination. No producer or SymPy import.',
    'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [source_path,normal_path]},
    'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
    'cases':records,'scope':'Four complete anchor/magnetic sources and finite full-ambient dimensions. Component distinctness and nonzero-residue coverage are separate geometric arguments.',
    'global_exhaustion_proved':False,'elapsed_seconds':time.monotonic()-started}
with Path(sys.argv[1]).open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps({'all_passed':True,'elapsed_seconds':result['elapsed_seconds']},indent=2))
