#!/usr/bin/env python3
"""Independent Fraction review of two full nil(3,1) long-top sources.
Arithmetic primitives reused from the dependent-source independent reviewer;
no producer, SymPy, or producer Jacobi/rank helper is imported.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import ast,hashlib,json,resource,sys,time

resource.setrlimit(resource.RLIMIT_CPU,(55,56))
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
started=time.monotonic();ROOT=Path.cwd();N=17
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
    names=dict(zip(('h','x1','x2','x3','x4')+tuple('n'+str(i) for i in range(11)),range(16)))
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

data_path='code/data/classification_nil31_long_sources.json'
report_path='reports/CLASSIFICATION_NIL31_LONG_SOURCES_20261107.json'
data=json.loads((ROOT/data_path).read_text());report=json.loads((ROOT/report_path).read_text())
source_path=data['source_report'];source=json.loads((ROOT/source_path).read_text())
assert hashlib.sha256((ROOT/data_path).read_bytes()).hexdigest()==report['input_sha256']
assert hashlib.sha256((ROOT/source_path).read_bytes()).hexdigest()==report['source_report_sha256']
producer_path='code/classification_modular_source_components.py'
assert hashlib.sha256((ROOT/producer_path).read_bytes()).hexdigest()==report['producer_sha256']
records=[]
for name in ['r_two_thirds','r_three_halves']:
    case_start=time.monotonic()
    spec=next(c for c in data['branches'] if c['name']==name)
    saved=next(c for c in report['branches'] if c['name']==name)
    assert spec==saved['full_source_specification']
    chart=source['cases'][spec['chart_index']]
    assert chart['kind']=='31' and chart['constant_index']==2 and chart['full_free_parameter_count']==11
    sub={5+int(k[1:]):parse_polynomial(v) for k,v in spec['substitution'].items()}
    entries=[compose(parse_polynomial(p),sub) for p in chart['full_six_brackets']]
    assert all(not compose(parse_polynomial(p),sub) for p in chart['all_distinct_jacobi_coefficients'])
    affine=tensor(dict(zip(combinations(range(1,5),2),entries)))
    assert all(not p for p in strict_jacobi(affine))
    pf=pfaffian(affine)
    assert pf==compose(parse_polynomial(chart['actual_pfaffian']),sub)
    assert pf==parse_polynomial(saved['universal_pfaffian'])
    expected=('-n8*(x1-x2*x3+x3**3/3+n2+3*n3*(x4-n9)/5)' if name=='r_two_thirds'
              else 'n4*(x2-x3**2/2+n7)')
    assert pf==parse_polynomial(expected)
    common=F(saved['integer_scale']);assert common==(15 if name=='r_two_thirds' else 2)
    Pgen=tensor({ij:scale(homogenize(p,2),common) for ij,p in zip(combinations(range(1,5),2),entries)})
    assert Pgen==decode(saved['universal_brackets'])
    assert all(not p for p in strict_jacobi(Pgen))
    assert pfaffian(Pgen)==scale(homogenize(pf,4),common**2)
    modular=[add(*(der(Pgen.get((i,j),{}),i) for i in range(5))) for j in range(5)]
    assert modular==[{},term(common,2),term(common,3),term(common,0),{}]
    # Check the canonically associated unimodular lift, retaining its projective gauge.
    Pi=tensor({(i,j):add(Pgen.get((i,j),{}),scale(add(mul(term(1,i),modular[j]),scale(mul(term(1,j),modular[i]),-1)),-F(1,5))) for i,j in pairs})
    assert all(not add(*(der(Pi.get((i,j),{}),i) for i in range(5))) for j in range(5))
    assert all(not p for p in projective_jacobi(Pi).values())
    params=[5+int(k[1:]) for k in spec['free_parameters']]
    assert len(params)==5
    values={5+int(k[1:]):F(v) for k,v in spec['witness'].items()}
    assert set(values)==set(params) and all(a==1 for a in values.values())
    P=tensor({ij:specialize(Pgen.get(ij,{}),values) for ij in pairs})
    P={ij:p for ij,p in P.items() if p};assert P==decode(saved['witness_brackets'])
    parameter_directions=[tensor({ij:specialize(der(Pgen.get(ij,{}),p),values) for ij in pairs}) for p in params]
    D=[[0]*150 for _ in row_basis]
    for col,(ij,m) in enumerate(column_basis):
        variation=tensor({ij:{m+(0,)*(N-5):F(1)}})
        mixed=tensor({kl:add(P.get(kl,{}),mul(term(1,16),variation.get(kl,{}))) for kl in pairs})
        for inds,p in projective_jacobi(mixed).items():
            for exp,a in p.items():
                if exp[16]==1:
                    assert a.denominator==1 and not any(exp[5:16])
                    D[row_index[inds,exp[:5]]][col]=int(a)
    vectors=[]
    for i,j in product(range(5),repeat=2):
        A=tensor({(k,l):add(term(int(l==i),k,j),term(-int(k==i),l,j)) for k,l in pairs})
        vectors.append(bivector_vector(A))
    for i,j in product(range(5),repeat=2):
        A=tensor({(k,l):add(scale(mul(term(1,j),der(P.get((k,l),{}),i)),-1),
            scale(P.get((j,l),{}),int(k==i)),scale(P.get((k,j),{}),int(l==i))) for k,l in pairs})
        vectors.append(bivector_vector(A))
    vectors.append(bivector_vector(P));vectors.extend(map(bivector_vector,parameter_directions))
    Fmatrix=list(map(list,zip(*vectors)))
    for row in D:
        nz=[(j,a) for j,a in enumerate(row) if a]
        assert all(sum(a*v[j] for j,a in nz)==0 for v in vectors)
    matrices={'Jacobi':D,'vertical':[r[:25] for r in Fmatrix],
        'orbit_scale':[r[:51] for r in Fmatrix],'family':Fmatrix}
    certificates={k:modular_minor(A) for k,A in matrices.items()}
    lower={k:c['rank_mod_prime'] for k,c in certificates.items()}
    assert lower=={'Jacobi':103,'vertical':24,'orbit_scale':47,'family':47}==saved['matrix_ranks_mod_prime']
    # Exact DF=0 and complementary nonzero minors prove all ranks over QQ.
    assert lower['Jacobi']+lower['family']==150
    assert all(sum(Fmatrix[i][5*j+j] for j in range(5))==0 for i in range(150))
    record={'name':name,'localization':spec['localization'],'all_five_parameters_retained':spec['free_parameters'],
        'source_substitution_exact':True,'all_six_saved_source_equations_zero':True,
        'universal_strict_Jacobi_exact_zero':True,'universal_Pfaffian_exact':expected,
        'unscaled_affine_modular_field':'x2*d/dx1 + x3*d/dx2 + d/dx3',
        'canonical_unimodular_lift_divergence_zero':True,'canonical_unimodular_lift_projective_Jacobi_zero':True,
        'integer_scale':int(common),'witness_and_all_actual_derivatives_reconstructed':True,
        'all_350_Jacobi_rows_and_150_coefficients_retained':True,
        'exact_D_times_all_56_actual_columns_zero':True,'matrix_ranks_QQ_by_exact_kernel_sandwich':lower,
        'rank_lower_certificates':certificates,'matrix_sha256':{k:matrix_hash(A) for k,A in matrices.items()},
        'projective_tangent_dimension_QQ':22,'projective_actual_family_dimension_lower':22,
        'projective_orbit_dimension_QQ':22,'smooth_full_component_dimension':22,
        'source_open_irreducible':True,'whole_90_variable_linear_reduction_reproved':False,
        'distinctness_from_other_components_checked_here':False,'seconds':time.monotonic()-case_start}
    records.append(record)
    print(json.dumps({k:record[k] for k in ['name','matrix_ranks_QQ_by_exact_kernel_sandwich','smooth_full_component_dimension','seconds']}),flush=True)

checker=Path(__file__).read_text()
paths=[data_path,report_path,source_path,producer_path]
result={'task_id':'E_multiple_hyperplane','all_passed':True,'field':'QQ',
    'method':'Independent standard-library Fraction sparse arithmetic and AST parser; no producer, SymPy, NumPy, legacy Jacobi or rank-helper import.',
    'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
    'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
    'branches':records,'elapsed_seconds':time.monotonic()-started,
    'scope':'Only nil(3,1) long-top r=2/3 and3/2 actual five-parameter sources. Universal identities, modular and unimodular gauges, full ambient dimensions and source componenthood. Completeness of the90-variable reduction and geometric distinctness are separate dependencies; r=-1 is not retested.',
    'global_exhaustion_proved':False}
with Path(sys.argv[1]).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'all_passed':True,'elapsed_seconds':result['elapsed_seconds']},indent=2))

