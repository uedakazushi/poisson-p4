#!/usr/bin/env python3
"""Independent complete rank-two/M(p)=0 plane first-jet review.

Only prior independent Fraction and sparse arithmetic is imported. Producer
snapshots are data, never imported. All 350 rows, full 30/9 sources and all rho
walls are reconstructed, including six actual shears. Laurent identities check
universal invertible normalizations rather than isolated parameter samples.
"""
import ast
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

from classification_unimodular_logarithmic_review import RF, det, rank, ud
import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_minimal_quartic_review import LaurentRing, determinant

ROOT = Path(__file__).resolve().parents[1]
INPUT = 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json'
OUTPUT = 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_REVIEW_20261107.json'
JOB = 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_REVIEW_JOB_20261107.json'
NAMES = ('a','b','c','z','w')
PAIRS = list(combinations(range(5),2))
ZERO = (0,)*8
RHO = RF((F(0),F(1)))


def atom(value=1, index=None, n=8):
    ex = [0]*n
    if index is not None: ex[index] = 1
    return {tuple(ex):RF(value)} if value else {}


def parse(text, names=NAMES, n=8):
    zero = (0,)*n
    def read(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:
            return atom(node.value,n=n)
        if isinstance(node,ast.Name):
            return atom(RHO,n=n) if node.id == 'tau' else atom(1,names.index(node.id),n)
        if isinstance(node,ast.UnaryOp):
            return q.scale(read(node.operand),-1 if isinstance(node.op,ast.USub) else 1)
        if isinstance(node,ast.BinOp):
            left = read(node.left)
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant) and node.right.value >= 0
                value = atom(1,n=n)
                for _ in range(node.right.value): value = q.mul(value,left)
                return value
            right = read(node.right)
            if isinstance(node.op,ast.Add): return q.add(left,right)
            if isinstance(node.op,ast.Sub): return q.add(left,q.scale(right,-1))
            if isinstance(node.op,ast.Mult): return q.mul(left,right)
            if isinstance(node.op,ast.Div):
                assert set(right) == {zero}
                return q.scale(left,RF(1)/right[zero])
        raise ValueError(ast.dump(node))
    return read(ast.parse(str(text),mode='eval').body)


def scalar(text):
    value = parse(text)
    assert not value or set(value) == {ZERO}
    return value.get(ZERO,RF())


def matrix(data): return [[scalar(value) for value in row] for row in data]
def encode(v): return {'numerator':list(map(str,RF(v).p)),'denominator':list(map(str,RF(v).q))}
def mat_encode(M): return [[encode(v) for v in row] for row in M]
def wedge(X,Y):
    return q.antisymmetric({(i,j):q.add(q.mul(X[i],Y[j]),q.scale(q.mul(X[j],Y[i]),-1)) for i,j in PAIRS})
def padd(*Ps): return q.antisymmetric({ij:q.add(*(P.get(ij,{}) for P in Ps)) for ij in PAIRS})
def pf4(P,indices):
    a,b,c,d=indices
    return q.add(q.mul(P.get((a,b),{}),P.get((c,d),{})),
                 q.scale(q.mul(P.get((a,c),{}),P.get((b,d),{})),-1),
                 q.mul(P.get((a,d),{}),P.get((b,c),{})))
def divisor(P):
    return q.add(*(q.scale(q.mul(atom(1,i),pf4(P,[j for j in range(5) if j!=i])),(-1)**i) for i in range(5)))


def source_bases(packet):
    base = []
    parameters = []
    for prefix,normal in [('l',3),('n',4)]:
        for i,j in combinations(range(3),2):
            for k in range(3):
                ex=[0]*8;ex[normal]+=1;ex[k]+=1
                base.append(q.antisymmetric({(i,j):{tuple(ex):RF(1)}}))
                parameters.append(prefix+str(len([v for v in parameters if v.startswith(prefix)])))
    for prefix,normal,j in [('d',(3,3),4),('b',(3,4),3),('e',(3,4),4),('c',(4,4),3)]:
        for i in range(3):
            ex=[0]*8
            for k in normal: ex[k]+=1
            base.append(q.antisymmetric({(i,j):{tuple(ex):RF(1)}}));parameters.append(prefix+str(i))
    second=[]
    for i,j in combinations(range(3),2):
        for normal in [(3,3),(3,4),(4,4)]:
            ex=[0]*8
            for k in normal:ex[k]+=1
            second.append(q.antisymmetric({(i,j):{tuple(ex):RF(1)}}))
    assert len(base)==30 and len(second)==9
    # Read the entire symbolic sources, with each parameter independent.
    for key,pars,basis in [('complete_source1',parameters,base),('complete_source2',['r'+str(i) for i in range(9)],second)]:
        names=NAMES+tuple(pars);n=len(names)
        for ij in PAIRS:
            source=parse(packet[key][str(ij)],names,n)
            extracted=[{} for _ in pars]
            for ex,value in source.items():
                assert sum(ex[5:])==1 and sum(ex[:5])==2
                k=next(k for k,v in enumerate(ex[5:]) if v)
                extracted[k][ex[:5]+(0,0,0)]=value
            assert extracted==[P.get(ij,{}) for P in basis]
    # Entire weight-one raw space is 36, modulo the six Euler gauges.
    raw=[]
    for ij,e in q.COLS:
        normal_count=e[3]+e[4]
        if normal_count-sum(k>=3 for k in ij)==1:
            raw.append(q.antisymmetric({ij:{e+(0,0,0):RF(1)}}))
    assert len(raw)==36
    X=[atom(1,k) for k in range(5)]
    gauges=[]
    for normal in (3,4):
        for i in range(3):
            V=[{} for _ in range(5)];V[i]=atom(1,normal)
            gauges.append(wedge(X,V))
    columns=[q.vec(P) for P in base+gauges]
    assert rank(list(map(list,zip(*columns))))==36
    assert all(not any(v for (ij,e),v in zip(q.COLS,q.vec(G)) if e[3]+e[4]-sum(k>=3 for k in ij)!=1) for G in gauges)
    # No linear Euler gauge has weight two. All nine raw directions are kept.
    raw2=[(ij,e) for ij,e in q.COLS if e[3]+e[4]-sum(k>=3 for k in ij)==2]
    assert len(raw2)==9 and rank(list(map(list,zip(*(q.vec(P) for P in second)))))==9
    return base,second,dict(weight_one_raw=36,Euler_gauge=6,canonical_weight_one=30,canonical_weight_two=9,
                            all_symbolic_source_parameters_retained=True)


def model(kind):
    x=[atom(1,i) for i in range(5)];a,b,c,z,w=x
    En=[{},{},{},z,w];U=[atom(),{},{},{},{}]
    if kind=='ND':
        T=[{},q.scale(b,2),q.scale(c,-2),{},{}]
        V=[{},{},{},q.mul(c,w),q.mul(b,z)];K=q.scale(q.mul(z,w),-2)
    else:
        V=[{},{},{},q.add(q.mul(b,z),q.mul(c,w)),q.scale(q.mul(b,w),-1)]
        if kind=='BSS_rho':T=[q.add(q.scale(a,2),q.scale(b,RHO)),q.scale(b,2),q.scale(c,-4),{},{}];K=q.scale(q.mul(z,w),-6)
        elif kind=='BN1':T=[c,{},q.scale(b,2),{},{}];K=q.scale(q.mul(w,w),-1)
        else:
            assert kind=='BN0';T=[{},{},q.scale(b,2),{},{}];K=q.scale(q.mul(w,w),-1)
    return padd(wedge(T,En),wedge(U,V),q.antisymmetric({(3,4):K}))


def operator(P,basis):
    assert not any(q.jacobi(P).values())
    cols=[]
    for B in basis:
        assert not any(q.jacobi(B).values())
        cols.append(q.jvec(q.jacobi(padd(P,B))))
    return list(map(list,zip(*cols)))


def shear(P,basis):
    x=[atom(1,i) for i in range(5)]; columns=[]
    # Actual derivative of the pushforward under x_i -> x_i + eps*z/w.
    for normal in (3,4):
        for i in range(3):
            Y=[{} for _ in range(5)];Y[i]=atom(1,normal)
            orbit=q.antisymmetric({ij:q.add(
                *(q.scale(q.mul(Y[k],q.der(P.get(ij,{}),k)),-1) for k in range(5)),
                *(q.mul(P.get((k,ij[1]),{}),q.der(Y[ij[0]],k)) for k in range(5)),
                *(q.mul(P.get((ij[0],k),{}),q.der(Y[ij[1]],k)) for k in range(5))) for ij in PAIRS})
            first=[orbit.get((k,3),{}).get((0,0,0,2,0,0,0,0),RF()) for k in range(3)]
            second=[orbit.get((k,4),{}).get((0,0,0,0,2,0,0,0),RF()) for k in range(3)]
            V=[q.add(q.scale(x[3],first[k]),q.scale(x[4],second[k])) for k in range(3)]+[{},{}]
            canonical=padd(orbit,wedge(x,V))
            coefficients=[]
            for B in basis:
                ij=next(ij for ij in PAIRS if B.get(ij))
                ex=next(iter(B[ij]));assert B[ij][ex]==1
                coefficients.append(canonical.get(ij,{}).get(ex,RF()))
            reconstructed=q.antisymmetric({ij:q.add(*(q.scale(B.get(ij,{}),c) for c,B in zip(coefficients,basis))) for ij in PAIRS})
            assert reconstructed==canonical
            columns.append(coefficients)
    return list(map(list,zip(*columns)))


def certificate(M,C):
    r=C['rank'];n=len(M[0]);k=n-r
    assert C['nullity']==k
    rr,cc=C['minor_rows'],C['minor_columns']
    assert len(rr)==len(cc)==r and len(set(rr))==len(set(cc))==r
    minor=det([[M[i][j] for j in cc] for i in rr])
    assert minor and minor==scalar(C['minor_determinant'])
    K=matrix(C['complete_kernel'])
    assert len(K)==n and all(len(row)==k for row in K)
    assert all(not sum((a*b for a,b in zip(row,col) if a and b),RF()) for row in M for col in zip(*K))
    assert not k or rank(K)==k
    return dict(shape=[len(M),n],rank=r,nullity=k,every_row_kernel_checked=True,
                independent_kernel=True,nonzero_minor=encode(minor),minor_rows=rr,minor_columns=cc)


def evaluate(v,rho):
    v=RF(v)
    def at(p):
        result=F(0)
        for c in reversed(p):result=result*rho+c
        return result
    return RF(at(v.p)/at(v.q))


def specialization(M,rho):return [[evaluate(v,rho) for v in row] for row in M]


def normalizations():
    R=LaurentRing(('alpha','beta','rb','rc','epsilon'),NAMES,invert=('alpha','beta','epsilon'))
    p=R.parse; x=[R.variable(v) for v in NAMES]
    def wed(X,Y):return {(i,j):R.add(R.mul(X[i],Y[j]),R.scale(R.mul(X[j],Y[i]),-1)) for i,j in PAIRS}
    def get(P,i,j):return P.get((i,j),{}) if i<j else R.scale(P.get((j,i),{}),-1) if i>j else {}
    def mod(kind):
        U=[p('1'),{},{},{},{}];En=[{},{},{},p('z'),p('w')]
        V=[{},{},{},p('c*w' if kind=='ND' else 'b*z+c*w'),p('b*z' if kind=='ND' else '-b*w')]
        if kind=='ND':T=[p('rb*b+rc*c'),p('2*alpha*b'),p('-2*alpha*c'),{},{}];K=p('-2*alpha*z*w')
        elif kind=='BSS':T=[p('2*alpha*a/3+rb*b+rc*c'),p('2*alpha*b/3'),p('2*beta*b-4*alpha*c/3'),{},{}];K=p('-2*alpha*z*w-beta*w**2')
        elif kind=='BN':T=[p('rb*b+rc*c'),{},p('2*beta*b'),{},{}];K=p('-beta*w**2')
        elif kind=='BN_eps':T=[p('epsilon*c'),{},p('2*b'),{},{}];K=p('-w**2')
        else:raise ValueError(kind)
        A=wed(T,En);B=wed(U,V)
        return {ij:R.add(A[ij],B[ij],K if ij==(3,4) else {}) for ij in PAIRS}
    def subst(poly,images):
        out={}
        for ex,coeff in poly.items():
            value={ex[:R.n]+(0,)*5:coeff}
            for k,n in enumerate(ex[R.n:]):
                for _ in range(n):value=R.mul(value,images[k])
            out=R.add(out,value)
        return out
    records=[]
    data=[
        ('ND','ND', ['a-rb*b/(2*alpha)+rc*c/(2*alpha)','b/alpha','c/alpha','z','w'],
         ['a+rb*b/2-rc*c/2','alpha*b','alpha*c','z','w'],'1/alpha','0'),
        ('BSS','BSS_rho', ['a+rc*c/(2*alpha)-rc*beta*b/(2*alpha**2)','3*b/alpha','3*c/alpha-3*beta*b/alpha**2','z+beta*w/(2*alpha)','w'],
         ['a-rc*c/6','alpha*b/3','alpha*c/3+beta*b/3','z-beta*w/(2*alpha)','w'],'3/alpha','rb+beta*rc/alpha'),
        ('BN','BN_eps', ['a-rb*c/(2*beta)','b','c/beta','z/beta','w'],
         ['a+rb*c/2','b','beta*c','beta*z','w'],'1','beta*rc'),
        ('BN_eps','BN1', ['a','b/epsilon','c','epsilon*z','w'],
         ['a','epsilon*b','c','z/epsilon','w'],'1/epsilon','1')]
    for source,target,G,inv,sigma,residual in data:
        G=list(map(p,G));inv=list(map(p,inv));P=mod(source)
        assert [subst(v,inv) for v in G]==x and [subst(v,G) for v in inv]==x
        L=[[R.derivative(v,k) for k in range(5)] for v in G]
        determinant_G=determinant(R,L);assert determinant_G
        pushed={ij:R.mul(p(sigma),subst(R.add(*(R.mul(L[ij[0]][k],get(P,k,l),L[ij[1]][l]) for k in range(5) for l in range(5))),inv)) for ij in PAIRS}
        if target=='BN_eps':
            expected={}
            for ij,v in mod('BN_eps').items():
                value={}
                for ex,c in v.items():
                    n=ex[4];f=list(ex);f[4]=0;t={tuple(f):c}
                    for _ in range(n):t=R.mul(t,p(residual))
                    value=R.add(value,t)
                expected[ij]=value
        else:
            hand=model(target)
            expected={}
            for ij in PAIRS:
                out={}
                for ex,c in hand.get(ij,{}).items():
                    assert c.q==(F(1),)
                    value={}
                    for power,coefficient in enumerate(c.p):
                        term=R.constant(coefficient)
                        for _ in range(power):term=R.mul(term,p(residual))
                        value=R.add(value,term)
                    term={}
                    for f,d in value.items():term[f[:R.n]+tuple(v+w for v,w in zip(f[R.n:],ex[:5]))]=d
                    out=R.add(out,term)
                expected[ij]=out
        assert pushed==expected,('normalization',source,target)
        records.append(dict(source=source,target=target,coordinate_map=[R.encode(v) for v in G],inverse=[R.encode(v) for v in inv],
                            determinant=R.encode(determinant_G),tensor_scalar=sigma,residual_parameter=residual,
                            all_ten_bivector_coefficients_identical=True,both_inverse_compositions=True))
    return records


def main():
    started=time.monotonic()
    signal.signal(signal.SIGALRM,lambda sig,frame:(_ for _ in ()).throw(TimeoutError('finite independent rank-two review180s')))
    signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    destination=ROOT/OUTPUT
    if destination.exists():raise FileExistsError(destination)
    job=json.loads((ROOT/JOB).read_text())
    for item in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'],item['path']
    packet=json.loads((ROOT/INPUT).read_text())
    assert packet['all_passed']
    assert hashlib.sha256(packet['producer_bytes'].encode()).hexdigest()==packet['producer_sha256']
    source1,source2,source_review=source_bases(packet)
    maps=normalizations()
    results=[]
    for record in packet['cases']:
        name=record['type'];P=model(name)
        assert all(P.get(ij,{})==parse(record['normal_form_P0'][str(ij)]) for ij in PAIRS)
        assert not any(q.jacobi(P).values())
        assert divisor(P)==parse(record['actual_D'])
        D1=operator(P,source1);D2=operator(P,source2);S=shear(P,source1)
        assert D1==matrix(record['complete_D1_350_by30'])
        assert D2==matrix(record['complete_D2_350_by9'])
        assert S==matrix(record['complete_actual_shear_30_by6'])
        c1=certificate(D1,record['generic_D1']);c2=certificate(D2,record['generic_D2']);cs=certificate(S,record['generic_actual_shear'])
        assert all(not sum((a*b for a,b in zip(row,col) if a and b),RF()) for row in D1 for col in zip(*S))
        assert all(RF(v).q==(F(1),) for M in [D1,D2,S] for row in M for v in row)
        affine={'ND':11,'BSS_rho':11,'BN1':10,'BN0':9}[name]
        d=c1['nullity']+c2['nullity'];bound=affine+d+6-1
        assert d==record['generic_higher_fiber_upper_bound'] and bound==record['generic_PGL_image_upper_bound'] and bound<=22
        walls=[]
        if name=='BSS_rho':
            value=RF(1)
            for M,C in [(D1,record['generic_D1']),(D2,record['generic_D2']),(S,record['generic_actual_shear'])]:
                value=value*det([[M[i][j] for j in C['minor_columns']] for i in C['minor_rows']])
            assert value==scalar(record['exception_polynomial']) and value.q==(F(1),)
            remainder=value.p;roots={}
            for rho in map(F,(-3,-1,1,3)):
                multiplicity=0
                while len(remainder)>1:
                    quotient,rest=ud(remainder,(-rho,F(1)))
                    if rest:break
                    remainder=quotient;multiplicity+=1
                assert multiplicity
                roots[str(rho)]=multiplicity
            assert len(remainder)==1 and remainder[0]
            supplied={}
            for factor,multiplicity in record['all_irreducible_factors']:
                polynomial=scalar(factor);assert polynomial.q==(F(1),) and len(polynomial.p)==2
                rho=-polynomial.p[0]/polynomial.p[1];supplied[str(rho)]=multiplicity
            assert roots==supplied
            assert set(roots)=={str(F(w['rho'])) for w in record['all_exceptional_roots']}
            for wall in record['all_exceptional_roots']:
                rho=F(wall['rho']);m1=specialization(D1,rho);m2=specialization(D2,rho);ms=specialization(S,rho)
                w1=certificate(m1,wall['D1']);w2=certificate(m2,wall['D2']);ws=certificate(ms,wall['actual_shear'])
                wd=w1['nullity']+w2['nullity'];wb=10+wd+6-1
                assert wd==wall['higher_fiber_upper_bound'] and wb==wall['PGL_image_upper_bound'] and wb<=22
                walls.append(dict(rho=str(rho),D1=w1,D2=w2,actual_shear=ws,affine_symbol_dimension=10,
                                  higher_fiber_upper_bound=wd,projective_image_upper_bound=wb))
            assert evaluate(value,F(0))
        else:assert not record['all_exceptional_roots']
        results.append(dict(type=name,full_350_rows_match=True,D1=c1,D2=c2,actual_shear=cs,
                            affine_symbol_dimension=affine,higher_fiber_upper_bound=d,projective_image_upper_bound=bound,
                            all_exceptional_roots=walls,
                            exceptional_factor_multiplicities=roots if name=='BSS_rho' else {},
                            symbol_divisor_and_projective_Poisson_checked=True))
        print(name,'independent kernels',c1['nullity'],c2['nullity'],'all rho walls',len(walls),'image upper',bound,flush=True)
    # Compatibility is audited as AST text only; the old rank-three mathematics is not rerun.
    old=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_20261107.json').read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(text).body if isinstance(n,ast.FunctionDef)}
    before=functions(old['producer_bytes']);after=functions(packet['producer_bytes'])
    assert all(after.get(k)==v for k,v in before.items())
    result=dict(task_id='G',all_passed=True,
        scope='Entire m_D=2 rank-zero plane O(1) first-jet stratum rank M=2 and M(p)=0 has projective image dimension at most22. No nonlinear/E/full-G classification.',
        complements=['M(p) != 0','rank M <= 1','O(0) common-image type','m_D >= 3','non-Frobenius first jets outside m_D=2'],
        complete_sources=source_review,universal_actual_normalization_identities=maps,cases=results,
        full_source_fiber_argument='D1 kernel bounds P1. For each fixed P1, P2 solves an affine linear D2 system; its solution fiber is empty or a translate of kerD2. Higher nonlinear equations only shrink this set.',
        no_cas_or_producer_import=True,old_rank_three_AST_compatibility=True,
        producer_snapshot_sha256=packet['producer_sha256'],producer_snapshot=packet['producer_bytes'],
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUTPUT,round(result['seconds'],3),'sec',result['maximum_RSS_KiB'],'KiB',flush=True)


if __name__=='__main__':main()
