#!/usr/bin/env python3
"""Independent full quadratic_A source, symplectic arc and actual SNC proof.

Standard-library sparse Laurent polynomials over the independently implemented
Fraction(lambda) field. No producer, SymPy or CAS imports; no negative search
receipt is needed. All six magnetic variables remain symbolic.
"""
import ast
from fractions import Fraction as F
from itertools import combinations,permutations
import hashlib,json,resource,signal,time
from pathlib import Path
from classification_unimodular_logarithmic_review import RF,rank as field_rank

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SOURCE_20261107.json'
ARC='reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SPECTRAL_ARC_20261107.json'
ORIGINAL='reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
OUTPUT='reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_REVIEW_20261107.json'
JOB='reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_REVIEW_JOB_20261107.json'
N=14;ZERO=(0,)*N;LAMBDA=RF((F(0),F(1)));PAIRS=list(combinations(range(4),2))
NAMES=('h','u1','u2','v1','v2','epsilon','I')+tuple('m'+str(i) for i in range(1,8))


def a(c=1,*inds):
    e=[0]*N
    for k in inds:e[k]+=1
    return {tuple(e):RF(c)} if c else {}
def add(*Ps):
    out={}
    for P in Ps:
        for e,c in P.items():out[e]=out.get(e,RF())+c
    return {e:c for e,c in out.items() if c}
def scale(P,c):return {e:c*v for e,v in P.items() if c*v}
def mul(*Ps):
    out=a()
    for P in Ps:
        new={}
        for e,c in out.items():
            for f,d in P.items():
                g=tuple(v+w for v,w in zip(e,f));new[g]=new.get(g,RF())+c*d
        out={e:c for e,c in new.items() if c}
    return out

def inv(P):
    assert len(P)==1
    e,c=next(iter(P.items()));assert all(not v or i in (5,6) for i,v in enumerate(e))
    return {tuple(-v for v in e):RF(1)/c}
def der(P,k):
    out={}
    for e,c in P.items():
        if e[k]:f=list(e);f[k]-=1;out[tuple(f)]=c*e[k]
    return out

def parse(text,env=None):
    names={name:a(1,i) for i,name in enumerate(NAMES)}
    names.update({'lambda':a(LAMBDA),'lam':a(LAMBDA),'alpha':a(LAMBDA-F(1,2)),
                  'beta':a(-LAMBDA-F(1,2)),'m0':add(a(1,6),scale(a(1,10),-1))})
    if env:names.update(env)
    def read(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return a(n.value)
        if isinstance(n,ast.Name):return names[n.id]
        if isinstance(n,ast.UnaryOp):return scale(read(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
        if isinstance(n,ast.BinOp):
            left=read(n.left)
            if isinstance(n.op,ast.Pow):
                power=ast.literal_eval(n.right);assert isinstance(power,int)
                if power<0:left=inv(left)
                return mul(*([left]*abs(power)))
            right=read(n.right)
            if isinstance(n.op,ast.Add):return add(left,right)
            if isinstance(n.op,ast.Sub):return add(left,scale(right,-1))
            if isinstance(n.op,ast.Mult):return mul(left,right)
            if isinstance(n.op,ast.Div):return mul(left,inv(right))
        raise ValueError(ast.dump(n))
    return read(ast.parse(str(text).replace('lambda','lam'),mode='eval').body)


def subst(P,images):
    out={}
    for ex,c in P.items():
        rest=list(ex)
        for k in images:rest[k]=0
        value={tuple(rest):c}
        for k,image in images.items():
            exponent=ex[k];base=inv(image) if exponent<0 else image
            for _ in range(abs(exponent)):value=mul(value,base)
        out=add(out,value)
    return out

def skew(data):
    out={}
    for (i,j),P in data.items():
        if P:out[i,j]=P;out[j,i]=scale(P,-1)
    return out

def plus(*Ps):return skew({ij:add(*(P.get(ij,{}) for P in Ps)) for ij in PAIRS})
def strict(P):
    return {ijk:add(*(mul(P.get((i,l),{}),der(P.get((j,k),{}),l+1))
                    for i,j,k in (ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1])) for l in range(4))) for ijk in combinations(range(4),3)}
def pf(P):return add(mul(P.get((0,1),{}),P.get((2,3),{})),scale(mul(P.get((0,2),{}),P.get((1,3),{})),-1),mul(P.get((0,3),{}),P.get((1,2),{})))
def decode(data,env=None):return skew({tuple(map(int,k.split(','))):parse(value,env) for k,value in data.items()})
def encode(P):return [[list(e),{'p':list(map(str,c.p)),'q':list(map(str,c.q))}] for e,c in sorted(P.items())]


def push(P,maps,back,scalar=None):
    M=[[der(v,j+1) for j in range(4)] for v in maps]
    return skew({ij:mul(scalar or a(),subst(add(*(mul(M[ij[0]][k],P.get((k,l),{}),M[ij[1]][l]) for k in range(4) for l in range(4))),dict(zip(range(1,5),back)))) for ij in PAIRS})
def determinant(M):
    out={}
    for perm in permutations(range(len(M))):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))
        out=add(out,scale(mul(*(M[i][j] for i,j in enumerate(perm))),sign))
    return out

def homogenize(P):
    out={}
    for ij in PAIRS:
        poly={}
        for ex,c in P.get(ij,{}).items():
            degree=sum(ex[1:5]);assert ex[0]==0 and degree<=2
            f=list(ex);f[0]=2-degree;poly[tuple(f)]=c
        out[ij[0]+1,ij[1]+1]=poly
    return skew(out)

def projective(P):
    strict5={ijk:add(*(mul(P.get((i,l),{}),der(P.get((j,k),{}),l))
                       for i,j,k in (ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1])) for l in range(5))) for ijk in combinations(range(5),3)}
    return {quad:add(*(scale(mul(a(1,quad[k]),strict5[quad[:k]+quad[k+1:]]),(-1)**k) for k in range(4))) for quad in combinations(range(5),4)}


def original_source():
    env={'q1':a(1,1),'q2':a(1,2),'x':a(1,3),'y':a(1,4)}
    return decode({'0,1':'m0+m1*x+m2*y+m3*x**2+m4*x*y+m5*y**2+beta*y*q1-alpha*x*q2',
                   '0,2':'alpha*x**2','0,3':'1+beta*x*y','1,2':'1+alpha*x*y','1,3':'beta*y**2'},env),env


def generic_parent(A,B,c):
    u1,u2,v1,v2=[a(1,k) for k in range(1,5)]
    return skew({(0,1):scale(mul(c,u1,u2),-1),(0,2):add(a(),mul(A,u1,v1),mul(B,u2,v2)),
                 (0,3):mul(c,u1,v2),(1,2):mul(add(B,scale(c,-1)),u2,v1),
                 (1,3):add(a(),mul(B,u2,v2)),(2,3):mul(add(B,scale(c,-1)),v1,v2)})


def snc(P,A,B):
    H=homogenize(P);weights=(3,0,1,0,1);limit={};all_weights=[]
    for ij in combinations(range(5),2):
        out={}
        for ex,c in H.get(ij,{}).items():
            weight=sum(v*w for v,w in zip(ex[:5],weights))-weights[ij[0]]-weights[ij[1]]
            assert weight>=0;all_weights.append(weight)
            if weight==0:out[ex]=c
        limit[ij]=out
    limit=skew(limit)
    for ij in combinations(range(5),2):
        for ex in limit.get(ij,{}):assert ex[ij[0]]==ex[ij[1]]==1
    Q=skew({ij:limit.get((ij[0]+1,ij[1]+1),{}) for ij in PAIRS})
    assert pf(Q)==scale(mul(A,B,a(1,1,2,3,4)),-1)
    assert not any(projective(limit).values())
    return dict(weights=list(weights),all_term_weights=sorted(set(all_weights)),limit_divisor='-A*B*h*u1*u2*v1*v2',all_period_c_free=True)


def five_source(packet):
    env=dict(zip(('a','f','p','q','r'),[a(1,k) for k in range(6,11)]))
    p=lambda t:parse(t,env)
    J=skew({(0,2):a(),(1,3):a()})
    Q=decode({'0,1':'p*u1*u2','0,2':'a*u1*v1+(q+r)*u2*v2','0,3':'q*u1*v2',
              '1,2':'(p+q+r)*u2*v1','1,3':'(p+q)*u1*v1+f*u2*v2','2,3':'r*v1*v2'},env)
    assert Q==decode(packet['five_parameter_Q'],env)
    JP=strict(plus(J,Q));JQ=strict(Q);assert JP==JQ
    coefficients=[];labels=[]
    for tri,value in JP.items():
        mons=sorted({ex[1:5] for ex in value},reverse=True)
        for mon in mons:
            coefficient={}
            for ex,c in value.items():
                if ex[1:5]==mon:
                    f=list(ex);f[1:5]=[0]*4;coefficient[tuple(f)]=c
            coefficients.append(coefficient);labels.append(dict(triple=list(tri),exponents=list(mon)))
    assert labels==packet['Jacobi_coefficient_labels']
    assert coefficients==[p(t) for t in packet['complete_Jacobi_coefficients']]
    generators=[p('(q+r)*(p+q)'),p('(q+r)*(q+r-f)'),p('(p+q)*(p+q-a)')]
    assert generators==[p(t) for t in packet['three_ideal_generators']]
    for row,value in zip(packet['Jacobi_to_generators'],coefficients):assert add(*(scale(g,F(c)) for c,g in zip(row,generators)))==value
    for row,value in zip(packet['generators_to_Jacobi'],generators):assert add(*(scale(g,F(c)) for c,g in zip(row,coefficients)))==value
    U=p('u1*v1');V=p('u2*v2')
    assert pf(plus(J,Q))==p('-1-(a+p+q)*u1*v1-(q+r+f)*u2*v2-a*(p+q)*(u1*v1)**2-a*f*u1*u2*v1*v2-(q+r)*f*(u2*v2)**2')
    return dict(displayed_five_parameter_source_all_Jacobi_coefficients=True,both_ideal_span_directions=True,
                exact_generators=['b*e','b*(b-f)','e*(e-a)'],whole_support_case_split=['b=e=0','e=0,f=b','b=0,a=e'],
                discovery_classification_not_needed_for_literal_curve=True)


def main():
    started=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    dest=ROOT/OUTPUT
    if dest.exists():raise FileExistsError(dest)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    source=json.loads((ROOT/SOURCE).read_text());arc=json.loads((ROOT/ARC).read_text())
    proof5=five_source(source)
    original,oldenv=original_source();assert not any(strict(original).values())
    raw=json.loads((ROOT/ORIGINAL).read_text())['cases'][7];assert raw['name']=='nil22_image_quadratic_A'
    data=skew({(i,j):parse(raw['full_six_magnetic_source'].get(str(i+1)+','+str(j+1),'0'),dict(oldenv,h=a(),**{f'x{k}':a(1,k) for k in range(1,5)})) for i,j in PAIRS})
    assert data==original
    norm=[parse(t,oldenv) for t in source['full_normalization_map']]
    inverse=[parse(t,oldenv) for t in source['full_normalization_inverse']]
    explicit_inverse=[parse(t,oldenv) for t in ['I*q1-m2/beta-m5*y/(2*beta)','I*q2+m1/alpha+m3*x/(2*alpha)-m4*y','x','y']]
    assert inverse==explicit_inverse
    X=[a(1,k) for k in range(1,5)]
    assert [subst(v,dict(zip(range(1,5),inverse))) for v in norm]==X
    assert [subst(v,dict(zip(range(1,5),norm))) for v in inverse]==X
    Nd=determinant([[der(v,j) for j in range(1,5)] for v in norm])
    assert Nd==parse(source['normalization_determinant'])==parse('I**(-2)')
    normalized=push(original,norm,inverse,a(1,6))
    target_old=decode({'0,1':'1+beta*y*q1-alpha*x*q2','0,2':'alpha*x**2','0,3':'1+beta*x*y','1,2':'1+alpha*x*y','1,3':'beta*y**2'},oldenv)
    assert normalized==target_old
    C=[parse(t,oldenv) for t in source['symplectic_standardization']]
    Cback=[parse(t) for t in ['u1','v1-v2','u2','v1']]
    target=push(target_old,C,Cback)
    assert target==decode(source['target_standard']) and not any(strict(target).values())
    J=decode({'0,2':'1','1,3':'1'});assert J==decode(source['symplectic_J'])
    # Entire actual three-parameter parent and its independent SNC degeneration.
    A,B,c=[a(1,k) for k in (6,7,8)];parent=generic_parent(A,B,c)
    assert not any(strict(parent).values())
    factors=[add(a(),mul(B,a(1,2,4))),add(a(),mul(B,a(1,2,4)),mul(A,a(1,1,3)))]
    assert pf(parent)==scale(mul(*factors),-1)
    snc_proof=snc(parent,A,B)
    # Exact source B of the accepted pair-zero theorem, with no period wall.
    order=(0,2,1,3);pairzero=skew({ij:parent[order[ij[0]],order[ij[1]]] for ij in PAIRS})
    expected=decode({'0,1':'1+A*u1*u2+B*v1*v2','2,3':'1+B*v1*v2','0,2':'-c*u1*v1','0,3':'c*u1*v2','1,2':'(c-B)*u2*v1','1,3':'(B-c)*u2*v2'},dict(A=A,B=B,c=c))
    # Reorder variables too: new (x,y,z,w)=(old u1,v1,u2,v2).
    reordered=skew({ij:subst(v,{1:a(1,1),2:a(1,3),3:a(1,2),4:a(1,4)}) for ij,v in pairzero.items() if ij in PAIRS})
    assert reordered==expected
    # Literal positive arc, independently transformed from actual parent.
    eps=a(1,5);alpha=parse('alpha');beta=parse('beta')
    actualA=parse('-epsilon*(2*epsilon+1)');actualB=parse('-epsilon*(epsilon+1)');actualc=mul(eps,alpha)
    trueparent=generic_parent(actualA,actualB,actualc)
    T=[parse(t) for t in ['u1+u2/epsilon','u2','v1','v2-v1/epsilon']]
    Tback=[parse(t) for t in ['u1-u2/epsilon','u2','v1','v2+v1/epsilon']]
    assert [subst(v,dict(zip(range(1,5),Tback))) for v in T]==X
    assert [subst(v,dict(zip(range(1,5),T))) for v in Tback]==X
    assert push(J,T,Tback)==J and determinant([[der(v,k) for k in range(1,5)] for v in T])==a()
    family=push(trueparent,T,Tback)
    literal=decode({'0,1':'alpha*u2**2-epsilon*alpha*u1*u2',
        '0,2':'1-epsilon*(2*epsilon+1)*u1*v1-epsilon*(epsilon+1)*u2*v2+beta*u2*v1',
        '0,3':'epsilon*alpha*u1*v2+(2*epsilon-beta)*u1*v1-alpha*u2*v2-u2*v1',
        '1,2':'epsilon*(beta-epsilon)*u2*v1','1,3':'1-epsilon*(epsilon+1)*u2*v2+alpha*u2*v1',
        '2,3':'(beta-epsilon)*v1*(v1+epsilon*v2)'})
    row=arc['curves'][1];assert row['name']=='rank3_rank5'
    assert family==literal==decode(row['actual_regular_affine_family'])
    assert all(ex[5]>=0 and coeff.q==(F(1),) for ij in PAIRS for ex,coeff in family.get(ij,{}).items())
    assert skew({ij:subst(P,{5:{}}) for ij,P in family.items() if ij in PAIRS})==target
    assert not any(strict(family).values())
    homogeneous=homogenize(family);assert homogeneous==decode(row['actual_regular_homogeneous_family']) and not any(projective(homogeneous).values())
    assert pf(family)==parse(row['actual_affine_Pfaffian'])
    # Pole bound before endpoint matching: the spectral variables have pole2,
    # the genuine single period has pole1; higher coefficients vanish at0.
    pole_parent=generic_parent(scale(A,-1),scale(B,-1),c)
    pole_family=push(pole_parent,T,Tback);poles={6:0,7:0,8:0}
    for ij in PAIRS:
        for ex,coefficient in pole_family.get(ij,{}).items():
            assert sum(ex[6:9]) in (0,1) and not any(ex[9:])
            for k in poles:
                if ex[k]:poles[k]=min(poles[k],ex[5])
    assert poles[6]>=-2 and poles[7]>=-2 and poles[8]>=-1
    solution=[RF(1),RF(2),RF(0),LAMBDA-F(1,2)]
    assert [parse(row['spectral_period_solution'][name]) for name in ['r0','r1','c0','c1']]==[a(v) for v in solution]
    M=[[parse(v).get(ZERO,RF()) for v in r] for r in row['matrix']];rhs=[parse(v).get(ZERO,RF()) for v in row['rhs']]
    assert field_rank(M)==4 and all(sum((v*w for v,w in zip(r,solution)),RF())==b for r,b in zip(M,rhs))
    binding={'a':actualA,'f':actualB,'p':scale(actualc,-1),'q':actualc,'r':add(actualB,scale(actualc,-1))}
    assert {k:parse(v) for k,v in row['true_source_periods'].items()}==binding
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,
        theorem='Every point of the entire seven-parameter quadratic_A source and its PGL saturation lies in the accepted finite closed MPS40 union. No new component and no local ideal claim.',
        all_six_magnetic_symbolic_source=True,normalizer_both_inverse_compositions=True,normalizer_determinant='I^-2',
        target_dense_open='alpha*beta*I !=0',all_target_complements_returned_by_irreducible_A7_closure=True,
        standardization_all_six_brackets=True,three_parameter_actual_parent_all_Jacobi=True,
        parent_Pfaffian='-(1+B*u2*v2)*(1+A*u1*v1+B*u2*v2)',actual_SNC_proof=snc_proof,
        accepted_pair_zero_B_source_all_parameter_identity=True,
        symplectic_matrix_and_inverse=True,all150_homogeneous_family_coefficients=True,
        all4_strict_and350_projective_Jacobi=True,all_lambda_endpoint=True,
        actual_parent_binding={k:encode(v) for k,v in binding.items()},literal_regular_family={str(ij):encode(family.get(ij,{})) for ij in PAIRS},
        spectral_and_period_pole_lower_bounds={str(k):v for k,v in poles.items()},four_unknown_actual_endpoint_system_rank=4,
        generic_epsilon_localizations=['epsilon','epsilon+1','2epsilon+1'],no_lambda_wall_for_literal_curve=True,
        displayed_five_source_review=proof5,negative_receipts_not_used_or_modified=True,
        no_producer_or_CAS_imports=True,inputs=job['exact_inputs_and_hashes'],
        source_snapshot=source['producer_bytes'],arc_snapshot=arc['producer_bytes'],
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUTPUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)

if __name__=='__main__':main()
