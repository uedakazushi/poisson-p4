#!/usr/bin/env python3
"""Independent Fraction review of the three complete MPS deformation families.
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
    names=dict(zip(('h','x','y','z','w','t','alpha','beta','gamma'),range(9)))
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

def determinant(A):
    B=[[F(a) for a in row] for row in A];value=F(1)
    for i in range(len(B)):
        j=next((j for j in range(i,len(B)) if B[j][i]),None)
        if j is None:return F(0)
        if j!=i:B[i],B[j]=B[j],B[i];value=-value
        pivot=B[i][i];value*=pivot
        for j in range(i+1,len(B)):
            ratio=B[j][i]/pivot
            B[j]=[a-ratio*b for a,b in zip(B[j],B[i])]
    return value

paths=['reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_M1_MPS_20261107.json',
       'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_M23_MPS_20261107.json',
       'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_20261107.json']
one=json.loads((ROOT/paths[0]).read_text());two=json.loads((ROOT/paths[1]).read_text());old=json.loads((ROOT/paths[2]).read_text())
for report in [one,two]:
    assert hashlib.sha256(report['producer_source'].encode()).hexdigest() in report['input_sha256'].values()
records=[];base_mons=[term(1,*f) for f in [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]]
def base_vector(p):return [p.get(next(iter(m)),F(0)) for m in base_mons]
for case_index,saved in enumerate([one,*two['cases']],1):
    begun=time.monotonic();name=saved['case'];P=parse_tensor(saved['full_polynomial_family'])
    assert all(sum(m[:5])==2 for p in P.values() for m in p)
    assert all(not p for p in strict_jacobi(P))
    Y=list(map(parse_polynomial,saved['coordinate_change']))
    assert Y[:3]==[term(1,i) for i in range(3)]
    T=[[der(y,i) for i in range(5)] for y in Y]
    det=add(mul(T[3][3],T[4][4]),scale(mul(T[3][4],T[4][3]),-1))
    assert det==term(1,5,5)
    t=term(1,5);one_tgamma=add(term(1),term(1,5,8))
    if case_index==1:
        Q=add(scale(mul(Y[3],Y[4]),F(1,2)),term(-1,0,1));k=-1
        p=term(1,5,6);r=add(term(F(1,2)),term(1,5,7))
        w1=add(term(1),scale(p,-1));canon=add(term(F(1,2),3,4),term(-1,0,1))
    else:
        k=-2 if case_index==2 else 1;sign=1 if case_index==2 else -1
        Q=add(mul(Y[3],Y[4]),term(sign,0,0));canon=add(term(1,3,4),term(sign,0,0))
        p=add(term(-1 if case_index==2 else 1),term(1,5,6))
        r=add(term(2 if case_index==2 else -1),term(1,5,7));w1=scale(p,-1)
    D1=[{},Y[1],{},mul(p,Y[3]),mul(w1,Y[4])]
    D2=[{},{},Y[2],mul(r,Y[3]),scale(mul(r,Y[4]),-1)]
    canonical=tensor({(i,j):mul(t,add(mul(D1[i],D2[j]),scale(mul(D2[i],D1[j]),-1))) for i,j in pairs})
    canonical[3,4]=add(canonical[3,4],scale(mul(mul(t,one_tgamma),Q),k));canonical[4,3]=scale(canonical[3,4],-1)
    pushed=tensor({(i,j):add(*(mul(mul(T[i][a],P.get((a,b),{})),T[j][b]) for a in range(5) for b in range(5))) for i,j in pairs})
    assert pushed==canonical
    Hess=[[der(der(canon,i),j).get((0,)*N,F(0)) for j in range(5)] for i in range(5)]
    assert rational_relations(Hess)['rank']==(4 if case_index==1 else 3)
    pf=add(mul(P.get((1,2),{}),P.get((3,4),{})),scale(mul(P.get((1,3),{}),P.get((2,4),{})),-1),mul(P.get((1,4),{}),P.get((2,3),{})))
    assert pf==scale(mul(mul(term(1,1,2),one_tgamma),Q),k)
    assert mul(term(1,0),pf)==parse_polynomial(saved['actual_quintic'])
    assert P.get((1,2))==term(1,5,1,2)
    limit={ij:specialize(p,{5:F(0)}) for ij,p in P.items()}
    limit={ij:p for ij,p in limit.items() if p};assert limit==parse_tensor(saved['limit_source'])
    assert specialize(Q,{5:F(0)})==parse_polynomial(saved['limit_quadric'])
    U=[limit.get((1,3),{}),limit.get((2,3),{})];V=[limit.get((1,4),{}),limit.get((2,4),{})]
    c=scale(der(limit[(3,4)],3),-1);d=scale(der(limit[(3,4)],4),-1)
    assert [add(*(add(mul(U[j],der(V[i],j+1)),scale(mul(V[j],der(U[i],j+1)),-1)) for j in range(2))) for i in range(2)]==[add(mul(c,U[i]),mul(d,V[i])) for i in range(2)]
    magnetic=[der(limit[(3,4)],j) for j in [6,7,8]]
    assert magnetic==list(map(parse_polynomial,saved['limit_magnetic_terms']))
    assert limit[(3,4)]==add(scale(mul(c,term(1,3)),-1),scale(mul(d,term(1,4)),-1),*(mul(term(1,j),p) for j,p in zip([6,7,8],magnetic)))
    cob=[add(*(mul(V[i],der(term(1,j),i+1)) for i in range(2)),mul(c,term(1,j))) for j in range(3)]
    cob += [add(mul(d,term(1,j)),*(scale(mul(U[i],der(term(1,j),i+1)),-1) for i in range(2))) for j in range(3)]
    assert cob==list(map(parse_polynomial,saved['all_six_translation_coboundaries']))
    cobmatrix=list(map(list,zip(*map(base_vector,cob))))
    assert rational_relations(cobmatrix)['rank']==3
    six=list(map(parse_polynomial,saved['full_magnetic_span_six_generators']))
    assert six[:3]==magnetic
    enlarged=list(map(list,zip(*map(base_vector,cob+six[3:]))))
    assert rational_relations(enlarged)['rank']==3
    matrix=list(map(list,zip(*map(base_vector,six))))
    assert [[str(v) for v in row] for row in matrix]==saved['full_magnetic_span_matrix']
    value=determinant(matrix);assert value==F(saved['full_magnetic_span_determinant'])
    # Relate the new limit anchors to the previously certified complete source by an exact constant GL2.
    old_case=next(cc for cc in old['cases'] if cc['case']==name)
    U0=list(map(parse_polynomial,old_case['anchor_U']));V0=list(map(parse_polynomial,old_case['anchor_V']))
    u0=base_vector(U0[0])+base_vector(U0[1]);v0=base_vector(V0[0])+base_vector(V0[1])
    i,j=next((i,j) for i,j in combinations(range(12),2) if u0[i]*v0[j]-u0[j]*v0[i])
    pivot=u0[i]*v0[j]-u0[j]*v0[i];coeffs=[]
    for W in [U,V]:
        w=base_vector(W[0])+base_vector(W[1]);a=(w[i]*v0[j]-w[j]*v0[i])/pivot;b=(u0[i]*w[j]-u0[j]*w[i])/pivot
        assert w==[a*u+b*v for u,v in zip(u0,v0)];coeffs.append([a,b])
    anchor_det=coeffs[0][0]*coeffs[1][1]-coeffs[0][1]*coeffs[1][0]
    assert anchor_det==-k
    record={'case':name,'all_parameters':['t','alpha','beta','gamma'],
        'full_polynomial_family_strict_Jacobi_identically_zero':True,
        'all_coefficients_quadratic_and_regular_at_t_zero':True,
        'coordinate_determinant':'t^2','canonical_pushforward_identity_exact':True,
        'canonical_quadric_rank':4 if case_index==1 else 3,
        'nonzero_residue':'a=t','rank_four_open':'t*(1+t*gamma)!=0',
        'universal_actual_quintic_exact':True,'full_t_zero_limit_exact':True,
        'limit_anchor_GL2_columns_old_basis':[[str(a) for a in row] for row in coeffs],
        'limit_anchor_GL2_determinant':str(anchor_det),
        'all_six_translation_coboundaries_exact':True,'translation_image_rank_QQ':3,
        'magnetic_span_matrix':[[str(a) for a in row] for row in matrix],
        'magnetic_span_determinant':str(value),'all_six_magnetic_coefficients_reached':True,
        'seconds':time.monotonic()-begun}
    records.append(record);print(json.dumps({'case':name,'all_passed':True,'magnetic_determinant':str(value),'seconds':record['seconds']}),flush=True)

checker=Path(__file__).read_text()
result={'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ',
    'method':'Independent standard-library Fraction sparse arithmetic; no producer or SymPy import.',
    'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
    'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
    'cases':records,'elapsed_seconds':time.monotonic()-started,
    'scope':'Universal actual MPS deformations of all three remaining independent-three-hyperplane conic magnetic sources, including canonical nonzero-residue recognition and every six-coefficient magnetic limit. Uses the previously proved all-support nonzero-residue MPS theorem.',
    'global_exhaustion_proved':False}
with Path(sys.argv[1]).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'all_passed':True,'elapsed_seconds':result['elapsed_seconds']},indent=2))

