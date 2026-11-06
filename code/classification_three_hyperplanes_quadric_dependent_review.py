from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import ast,hashlib,json,resource,sys,time

resource.setrlimit(resource.RLIMIT_CPU,(55,56))
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
started=time.monotonic();ROOT=Path.cwd();N=6
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
def bivector_vector(A):return [A.get(ij,{}).get(m+(0,),F(0)) for ij,m in column_basis]
def jacobi_vector(A):
    J=projective_jacobi(A)
    return [J[inds].get(m+(0,),F(0)) for inds,m in row_basis]
def parse_polynomial(text):
    names=dict(zip(('h','x','y','z','w'),range(5)))
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

source_path='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_20261107.json'
quadrics_path='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_NORMALS_20261107.json'
source=json.loads((ROOT/source_path).read_text());quadrics=json.loads((ROOT/quadrics_path).read_text())
P=parse_tensor(source['point_brackets'])
parameters=[parse_tensor(d) for d in source['all8_actual_parameter_directions']]
normals=[parse_tensor(d) for d in quadrics['normals']]
assert len(parameters)==8 and len(normals)==2
expected_parameters=[tensor({(2,4):term(1,0,2)}),tensor({(2,3):term(1,2,2),(2,4):term(1,1,2),(3,4):term(1,2,4)})]
expected_parameters += [tensor({(3,4):term(-1,*factors)}) for factors in [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]]
assert parameters==expected_parameters
f=add(term(1,1,1),term(-1,0,1))
q=add(f,term(1,0,2),term(1,1,2))
assert P==tensor({(1,3):f,(2,3):add(term(1,2,2),term(2,1,2),term(-1,0,2)),(2,4):q,
 (3,4):add(term(1,2,4),*(term(-a,*factors) for a,factors in zip((1,2,3,5,7,11),[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)])))})
pf=add(mul(P.get((1,2),{}),P.get((3,4),{})),scale(mul(P.get((1,3),{}),P.get((2,4),{})),-1),mul(P.get((1,4),{}),P.get((2,3),{})))
assert pf==scale(mul(f,q),-1)
assert all(not p for p in projective_jacobi(P).values())
D=[[0]*150 for _ in row_basis]
for col,(ij,m) in enumerate(column_basis):
    variation=tensor({ij:{m+(0,):F(1)}})
    mixed=tensor({kl:add(P.get(kl,{}),mul(term(1,5),variation.get(kl,{}))) for kl in pairs})
    for inds,p in projective_jacobi(mixed).items():
        for exp,a in p.items():
            if exp[5]==1:
                assert a.denominator==1
                D[row_index[inds,exp[:5]]][col]=int(a)
vectors=[];direction_names=[]
for i,j in product(range(5),repeat=2):
    V=tensor({(k,l):add(term(int(l==i),k,j),term(-int(k==i),l,j)) for k,l in pairs})
    vectors.append(bivector_vector(V));direction_names.append(f'vertical_{i}_{j}')
for i,j in product(range(5),repeat=2):
    action=tensor({(k,l):add(scale(mul(term(1,j),der(P.get((k,l),{}),i)),-1),
                                scale(P.get((j,l),{}),int(k==i)),
                                scale(P.get((k,j),{}),int(l==i))) for k,l in pairs})
    vectors.append(bivector_vector(action));direction_names.append(f'linear_action_{i}_{j}')
vectors.append(bivector_vector(P));direction_names.append('scale_P')
vectors.extend(map(bivector_vector,parameters));direction_names.extend('parameter_'+str(i) for i in range(8))
source_matrix=[list(row) for row in zip(*vectors)]
vectors.extend(map(bivector_vector,normals));direction_names.extend('normal_'+str(i) for i in range(2))
kernel_matrix=[list(row) for row in zip(*vectors)]
for row in D:
    nonzero=[(j,a) for j,a in enumerate(row) if a]
    assert all(sum(a*v[j] for j,a in nonzero)==0 for v in vectors)

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
d_certificate=modular_minor(D);source_certificate=modular_minor(source_matrix);kernel_certificate=modular_minor(kernel_matrix)
assert (d_certificate['rank_mod_prime'],source_certificate['rank_mod_prime'],kernel_certificate['rank_mod_prime'])==(99,49,51)

# Exact upper certificate for the 59-column genuine-source tangent matrix.
B=[[F(a) for a in row] for row in source_matrix];pivots=[];rank=0
for j in range(59):
    k=next((i for i in range(rank,150) if B[i][j]),None)
    if k is None:continue
    B[rank],B[k]=B[k],B[rank];a=B[rank][j];B[rank]=[v/a for v in B[rank]]
    for i in range(150):
        a=B[i][j]
        if i!=rank and a:B[i]=[x-a*y for x,y in zip(B[i],B[rank])]
    pivots.append(j);rank+=1
assert rank==49
free=[j for j in range(59) if j not in pivots];relations=[]
for j in free:
    relation=[F(0)]*59;relation[j]=F(1)
    for i,k in enumerate(pivots):relation[k]=-B[i][j]
    assert all(sum(x*y for x,y in zip(row,relation))==0 for row in source_matrix)
    relations.append(relation)
assert len(relations)==10

normal_pairs=list(combinations_with_replacement(range(2),2))
squares=[jacobi_vector(R) for R in normals];qcols=[]
for i,j in normal_pairs:
    if i==j:qcols.append(squares[i])
    else:
        value=jacobi_vector(tensor({kl:add(normals[i].get(kl,{}),normals[j].get(kl,{})) for kl in pairs}))
        qcols.append([a-b-c for a,b,c in zip(value,squares[i],squares[j])])
Q=[list(row) for row in zip(*qcols)]
functionals=[];pairings=[]
for saved in quadrics['left_cokernel_functionals']:
    functional={int(i):F(a) for i,a in saved['nonzero_coefficients'].items()}
    assert all(sum(a*D[i][j] for i,a in functional.items())==0 for j in range(150))
    pairing=[sum(a*Q[i][j] for i,a in functional.items()) for j in range(3)]
    assert pairing==list(map(F,saved['quadratic_coefficients']))
    functionals.append(functional);pairings.append(pairing)
assert len(pairings)==3
p=pairings
determinant=p[0][0]*(p[1][1]*p[2][2]-p[1][2]*p[2][1])-p[0][1]*(p[1][0]*p[2][2]-p[1][2]*p[2][0])+p[0][2]*(p[1][0]*p[2][1]-p[1][1]*p[2][0])
assert determinant
checker_source=Path(__file__).read_text()
result={
 'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ',
 'method':'Independent Fraction sparse polynomial arithmetic, AST literal parser and rational elimination; no SymPy or producer import.',
 'independent_checker_source':checker_source,'source_sha256':hashlib.sha256(checker_source.encode()).hexdigest(),
 'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [source_path,quadrics_path]},
 'literal_input_snapshot':{'point':source['point_brackets'],'eight_actual_parameter_directions':source['all8_actual_parameter_directions'],'two_normals':quadrics['normals']},
 'point_and_all_eight_directions_match_displayed_source_exactly':True,
 'Pfaffian_equal_minus_fQ_exactly':True,
 'Jacobi_row_basis':[[list(inds),list(m)] for inds,m in row_basis],
 'bivector_column_basis':[[list(ij),list(m)] for ij,m in column_basis],
 'full_D_shape':[350,150],'full_D_sha256':hashlib.sha256(json.dumps(D,separators=(',',':')).encode()).hexdigest(),
 'D_times_all_source_and_normal_directions_exact_zero':True,'direction_names':direction_names,
 'full_D_rank_lower_certificate':d_certificate,'source_rank_lower_certificate':source_certificate,
 'source_plus_normals_rank_certificate':kernel_certificate,
 'source_QQ_rank_upper_certificate':{'columns':59,'rank':49,'free_columns':free,'ten_independent_kernel_relations':[[str(a) for a in rel] for rel in relations]},
 'full_D_rank_QQ_by_sandwich':99,'full_kernel_dimension_QQ':51,'actual_source_tangent_rank_QQ':49,'normal_quotient_dimension':2,
 'normal_monomial_order':['n0^2','n0*n1','n1^2'],
 'three_quadratic_vectors_sparse':[{str(i):str(a) for i,a in enumerate(col) if a} for col in qcols],
 'quadratic_obstruction_rank_QQ':3,
 'three_left_null_functionals':[{str(i):str(a) for i,a in functional.items()} for functional in functionals],
 'all_three_functionals_annihilate_all_150_columns_exactly':True,
 'three_pairings':[[str(a) for a in row] for row in pairings],
 'pairing_matrix_determinant':str(determinant),
 'genuine_quadratic_initial_rowspace':['n0^2','n0*n1','n1^2'],
 'quadratic_support_zero_dimensional':True,
 'scope':'Independent exact arithmetic for the true two-normal slice. Geometric Krull/source-image and distinctness arguments are separate dependencies.',
 'full_initial_ideal_or_nil_thickness_claimed':False,'global_exhaustion_proved':False,
 'elapsed_seconds':time.monotonic()-started,
}
output=Path(sys.argv[1])
with output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps({k:result[k] for k in ['all_passed','full_D_rank_QQ_by_sandwich','actual_source_tangent_rank_QQ','normal_quotient_dimension','quadratic_obstruction_rank_QQ','pairing_matrix_determinant','elapsed_seconds']},indent=2))
