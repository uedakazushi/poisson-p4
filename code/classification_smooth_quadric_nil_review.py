#!/usr/bin/env python3
"""Independent Fraction review of all three nilpotent smooth-quadric cofactor sources.
Arithmetic primitives reused from the dependent-source independent reviewer;
no producer, SymPy, or producer Jacobi/rank helper is imported.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import ast,hashlib,json,resource,sys,time

resource.setrlimit(resource.RLIMIT_CPU,(55,56))
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
started=time.monotonic();ROOT=Path.cwd();N=21
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
    names=dict(zip(tuple('x'+str(i) for i in range(5))+tuple('n'+str(i) for i in range(15)),range(20)))
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

paths=['reports/CLASSIFICATION_SMOOTH_QUADRIC_'+kind+'_LINEAR_20261107.json' for kind in ['NIL5','NIL311','NIL221']]
records=[]
for path,count in zip(paths,[5,11,15]):
    begun=time.monotonic();saved=json.loads((ROOT/path).read_text());kind=saved['case']
    assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['producer_sha256']
    q=parse_polynomial(saved['q']);P=tensor({tuple(map(int,key)):parse_polynomial(value) for key,value in saved['brackets'].items()})
    M=[[0]*5 for _ in range(5)]
    if kind=='nil5':
        for i in range(4):M[i+1][i]=i+1
    elif kind=='nil311':M[1][0]=1;M[2][1]=2
    else:assert kind=='nil221';M[1][0]=1;M[3][2]=1
    V=[add(*(term(M[i][j],j) for j in range(5))) for i in range(5)]
    assert not add(*(mul(V[i],der(q,i)) for i in range(5)))
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
    record={'case':kind,'full_linear_matrix_shape':[len(full),150],
        'full_linear_rank_QQ_by_kernel_sandwich':150-count,'source_dimension':count,
        'source_all_linear_constraints_exact':True,'all_source_directions_independent':True,
        'linear_rank_lower_certificate':cert,'source_rank_lower_certificate':dircert,
        'Euler_contracted_Pfaffian_reproduced_before_Jacobi':True}
    if kind=='nil5':
        targets=['(n4-12)*(n4-4)','n2*(3*n4+4)',
            '16*n2**2+9*n3*n4+12*n3','15*n0*n4+20*n0+8*n1*n2+16*n3**2']
        raw,witnesses=jacobi_coefficient_witnesses(P,targets)
        checks=[]
        for value in [4,12]:
            sub={5:F(0),7:F(0),8:F(0),9:F(value)}
            assert all(not specialize(p,sub) for p in raw)
            pp=tensor({ij:specialize(P.get(ij,{}),sub) for ij in pairs})
            ff=five_pfaffians(pp);actual=specialize(potential,sub)
            target=scale(mul(term(1,0),mul(q,q)),-F(1,2)) if value==4 else {}
            assert actual==target
            assert all(der(actual,i)==ff[i] for i in range(5))
            checks.append({'n4':value,'n0_n2_n3_zero':True,'n1_free':True,
                'all_Jacobi_zero':True,'actual_divisor':'-x0*q^2/2' if value==4 else '0',
                'all_five_gradient_identities':True})
        record.update(necessary_Jacobi_literal_coefficients=witnesses,all_support_branches=checks,
            irreducible_cubic_factor_impossible=True)
    elif kind=='nil311':
        quotient=divide_exact(potential,mul(term(1,0),q))
        assert all(sum(m[:5])==2 for m in quotient)
        record.update(Euler_Pfaffian_divisible_by_x0_q_identically=True,
            remaining_factor_coordinate_degree=2,irreducible_cubic_factor_impossible=True,
            residual_Jacobi_support_not_solved=True)
    else:
        polynomial=parse_polynomial('(n8*x0**2+2*n12*x0*x2-n1*x2**2)*x4/2+(n5*x0**3+3*n7*x0**2*x2+3*n10*x0*x2**2-n0*x2**3)/3')
        assert P[(0,2)]==der(polynomial,4) and P[(0,4)]==scale(der(polynomial,2),-1) and P[(2,4)]==der(polynomial,0)
        raw,witnesses=jacobi_coefficient_witnesses(P,['n1*n8+n12**2','n0*n8+n1*n7+2*n10*n12','n1*n5-n10*n8+2*n12*n7'])
        zero_params={5+i:F(0) for i in [0,1,5,7,8,10,12]}
        Pzero=tensor({ij:specialize(P.get(ij,{}),zero_params) for ij in pairs})
        assert all(not p for p in five_pfaffians(Pzero))
        # All invariant cubics, independently retaining the full35-monomial space.
        cubic_mons=monomials(3);Ncols=[]
        for m in cubic_mons:
            p={m+(0,)*(N-5):F(1)}
            Np=add(*(mul(V[i],der(p,i)) for i in range(5)))
            Ncols.append([Np.get(mm+(0,)*(N-5),F(0)) for mm in cubic_mons])
        invmatrix=list(map(list,zip(*Ncols)));invcert=modular_minor(invmatrix);assert invcert['rank_mod_prime']==22
        base_cubics=[term(1,*inds) for inds in combinations_with_replacement([0,2,4],3)]
        cubics=base_cubics+[mul(q,term(1,i)) for i in [0,2,4]]
        assert all(not add(*(mul(V[i],der(p,i)) for i in range(5))) for p in cubics)
        cubmatrix=[[p.get(mm+(0,)*(N-5),F(0)) for p in cubics] for mm in cubic_mons]
        cubcert=modular_minor(cubmatrix);assert cubcert['rank_mod_prime']==13
        rank_one=specialize(polynomial,{5+i:F(0) for i in [0,1,10,12]})
        quotient=divide_exact(rank_one,term(1,0,0));assert all(sum(m[:5])==1 for m in quotient)
        zero_A=specialize(polynomial,{5+i:F(0) for i in [1,8,12]})
        assert all(m[4]==0 and not m[1] and not m[3] for m in zero_A)
        record.update(base_Jacobian_cubic_exact=True,
            base_Jacobian_cubic='(n8*u^2+2*n12*u*v-n1*v^2)*w/2+(n5*u^3+3*n7*u^2*v+3*n10*u*v^2-n0*v^3)/3',
            necessary_Jacobi_literal_coefficients=witnesses,
            F_zero_implies_all_five_Pfaffian_minors_zero=True,
            all_invariant_cubics_kernel_dimension=13,invariant_cubic_rank_certificate=invcert,
            invariant_cubic_basis_certificate=cubcert,
            all_invariant_cubics='C3(u,v,w)+q*ell(u,v,w)',
            binary_A_discriminant_zero_necessary=True,
            normalized_rank_one_F_divisible_by_u_squared=True,
            A_zero_F_binary_cubic=True,
            irreducible_cubic_factor_impossible_given_Casimir_argument=True,
            full_fifteen_parameter_Jacobi_support_not_solved=True)
    record['seconds']=time.monotonic()-begun;records.append(record)
    print(json.dumps({'case':kind,'linear_rank':150-count,'source_dimension':count,'seconds':record['seconds']}),flush=True)

checker=Path(__file__).read_text()
result={'task_id':'E_reduced_nonnormal','all_passed':True,'field':'QQ',
    'method':'Independent standard-library Fraction sparse arithmetic; no producer or SymPy import.',
    'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
    'source_sha256':hashlib.sha256(checker.encode()).hexdigest(),'independent_checker_source':checker,
    'cases':records,'elapsed_seconds':time.monotonic()-started,
    'scope':'All three full150 linear nilpotent cofactor sources. Only necessary support conditions and Pfaffian factors needed to exclude smoothq times irreducible cubic; the geometric nil221 Casimir implication is reviewed separately.',
    'global_exhaustion_proved':False}
with Path(sys.argv[1]).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'all_passed':True,'elapsed_seconds':result['elapsed_seconds']},indent=2))

