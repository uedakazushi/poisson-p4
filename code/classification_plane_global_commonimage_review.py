#!/usr/bin/env python3
"""Independent Fraction review of all four whole common-image plane sources.

Only earlier independent sparse/Fraction libraries are imported. The immutable
producer snapshot is data, never imported. Reconstructs all350 rows and all30/9
directions, and checks the entire18-parameter K0 tensor, not point samples.
"""
import ast
from fractions import Fraction as F
from itertools import combinations
import hashlib,json,resource,signal,time
from pathlib import Path
import classification_plane_rank_zero_rank_two_review as r

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_global_commonimage_review.json'
INPUT='reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_20261107.json'
OUTPUT='reports/CLASSIFICATION_PLANE_GLOBAL_COMMONIMAGE_REVIEW_20261107.json'
KINDS=('SS_COMMON_IMAGE_K1','SS_COMMON_IMAGE_K0','NIL_COMMON_IMAGE_K1','NIL_COMMON_IMAGE_K0')


def sources():
    first=[]
    for normal in (3,4):
        for i,j in combinations(range(3),2):
            for k in range(3):
                first.append(r.q.antisymmetric({(i,j):r.q.mul(r.atom(1,normal),r.atom(1,k))}))
    for normal,j in [((3,3),4),((3,4),3),((3,4),4),((4,4),3)]:
        for i in range(3):
            first.append(r.q.antisymmetric({(i,j):r.q.mul(r.atom(1,normal[0]),r.atom(1,normal[1]))}))
    second=[]
    for i,j in combinations(range(3),2):
        for normal in ((3,3),(3,4),(4,4)):
            second.append(r.q.antisymmetric({(i,j):r.q.mul(r.atom(1,normal[0]),r.atom(1,normal[1]))}))
    assert len(first)==30 and len(second)==9
    x=[r.atom(1,i)for i in range(5)]
    gauges=[]
    for normal in (3,4):
        for i in range(3):
            Y=[{}for _ in range(5)];Y[i]=x[normal]
            gauges.append(r.wedge(x,Y))
    weight=lambda ij,e:e[3]+e[4]-sum(i>=3 for i in ij)
    raw1=[(ij,e)for ij,e in r.q.COLS if weight(ij,e)==1]
    raw2=[(ij,e)for ij,e in r.q.COLS if weight(ij,e)==2]
    columns=[r.q.vec(P)for P in first+gauges]
    assert len(raw1)==36 and r.rank(list(map(list,zip(*columns))))==36
    assert all(not v or weight(ij,e)==1 for P in first+gauges for (ij,e),v in zip(r.q.COLS,r.q.vec(P)))
    assert len(raw2)==9 and r.rank(list(map(list,zip(*(r.q.vec(P)for P in second)))))==9
    return first,second,dict(full_weight_one_raw=36,Euler_gauge=6,canonical_weight_one=30,
                             full_weight_two=9,no_weight_two_gauge=True)


def model(kind):
    a,b,c,z,w=[r.atom(1,i)for i in range(5)]
    T=[b,{},{},{},{}];V=[c,{},{},{},{}];En=[{},{},{},z,w]
    ss=kind.startswith('SS');normal=[{},{},{},z,r.q.scale(w,-1)]if ss else[{},{},{},{},z]
    K=r.q.mul(z,w if ss else z)if kind.endswith('K1')else{}
    return r.padd(r.wedge(T,En),r.wedge(V,normal),r.q.antisymmetric({(3,4):K}))


def as_constant(value):
    value=r.RF(value);assert len(value.p)<=1 and value.q==(F(1),)
    return value.p[0]if value.p else F()


def coordinates(P,basis):
    coefficients=[]
    for B in basis:
        ij=next(ij for ij in r.PAIRS if B.get(ij))
        e=next(iter(B[ij]));assert B[ij][e]==1
        coefficients.append(P.get(ij,{}).get(e,r.RF()))
    rebuilt=r.padd(*[{ij:r.q.scale(B.get(ij,{}),coefficient)for ij in r.PAIRS}
                     for B,coefficient in zip(basis,coefficients)])
    assert all(rebuilt.get(ij,{})==P.get(ij,{})for ij in r.PAIRS)
    return coefficients


def explicit_zero_K_kernels(first,second):
    a,b,c,z,w=[r.atom(1,i)for i in range(5)]
    A=[]
    for normal in (z,w):
        for component in (1,2):
            for base in (b,c):A.append(r.q.antisymmetric({(0,component):r.q.mul(normal,base)}))
    for normal,component in [((z,w),3),((w,w),3),((z,z),4),((z,w),4)]:
        A.append(r.q.antisymmetric({(0,component):r.q.mul(*normal)}))
    B=[r.q.antisymmetric({(0,j):r.q.mul(*normal)})
       for j in (1,2)for normal in ((z,z),(z,w),(w,w))]
    C1=list(map(list,zip(*(coordinates(P,first)for P in A))))
    C2=list(map(list,zip(*(coordinates(P,second)for P in B))))
    assert r.rank(C1)==12 and r.rank(C2)==6
    assert all(common_direction(P)for P in A+B)
    return A,B,C1,C2


def common_direction(P):
    return all(not P.get(ij)for ij in r.PAIRS if 0 not in ij)and all(not r.q.der(P.get(ij,{}),0)for ij in r.PAIRS)


def entire_zero_K_tensor(P0,first,second):
    directions=first+second;n=5+len(directions)
    P={ij:{}for ij in r.PAIRS}
    for k,B in enumerate([P0]+directions):
        for ij in r.PAIRS:
            for exponent,value in B.get(ij,{}).items():
                e=list(exponent[:5])+[0]*len(directions)
                if k:e[4+k]=1
                P[ij]=r.q.add(P[ij],{tuple(e):as_constant(value)})
    P=r.q.antisymmetric(P)
    raw={}
    for i,j,k in combinations(range(5),3):
        raw[i,j,k]=r.q.add(*(r.q.mul(P.get((a,l),{}),r.q.der(P.get((b,c),{}),l))
                            for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(5)))
    assert not any(raw.values())
    assert common_direction(P)
    x=[r.atom(1,i,n)for i in range(5)]
    full350={quad:r.q.add(*(r.q.scale(r.q.mul(x[quad[k]],raw[quad[:k]+quad[k+1:]]),(-1)**k)
                          for k in range(4)))for quad in r.q.QUADS}
    assert not any(full350.values())
    pf=[r.pf4(P,[j for j in range(5)if j!=i])for i in range(5)]
    assert not any(pf)
    return dict(entire_higher_parameters=len(directions),strict_ten_Jacobi_polynomials_zero=True,
                projective_all350_zero=True,all_five_Pfaffians_zero=True,
                common_constant_direction='partial_a',all_coefficients_independent_of_a=True,
                every_parameter_value_rank_at_most_two=True)


def normalizations():
    """Actual Laurent normalizers of both nonzero K coefficients and inverses."""
    R=r.LaurentRing(('k',),('a','b','c','z','w'),invert=('k',))
    p=R.parse;x=list(map(R.variable,('a','b','c','z','w')))
    def substitute(poly,images):
        result={}
        for e,c in poly.items():
            term={e[:1]+(0,)*5:c}
            for i,power in enumerate(e[1:]):
                for _ in range(power):term=R.mul(term,images[i])
            result=R.add(result,term)
        return result
    records=[]
    for ss,G,Gi,sigma in [
        (True,['a','b/k','c/k','z','w'],['a','k*b','k*c','z','w'],'1/k'),
        (False,['a','b','c/k','k*z','w'],['a','b','k*c','z/k','w'],'1')]:
        G=list(map(p,G));Gi=list(map(p,Gi));assert [substitute(y,Gi)for y in G]==x
        assert [substitute(y,G)for y in Gi]==x
        P={(0,3):p('(b+c)*z'if ss else'b*z'),(0,4):p('(b-c)*w'if ss else'b*w+c*z'),(3,4):p('k*z*w'if ss else'k*z**2')}
        def get(i,j):return P.get((i,j),{})if i<j else R.scale(P.get((j,i),{}),-1)if i>j else{}
        A=[[R.derivative(y,j)for j in range(5)]for y in G]
        pushed={ij:R.mul(p(sigma),substitute(R.add(*(R.mul(A[ij[0]][k],get(k,l),A[ij[1]][l])
                               for k in range(5)for l in range(5))),Gi))for ij in r.PAIRS}
        expected={(0,3):p('(b+c)*z'if ss else'b*z'),(0,4):p('(b-c)*w'if ss else'b*w+c*z'),(3,4):p('z*w'if ss else'z**2')}
        assert all(pushed[ij]==expected.get(ij,{})for ij in r.PAIRS)
        records.append(dict(normal_type='semisimple'if ss else'nilpotent',forward=[str(v)for v in (['a','b/k','c/k','z','w']if ss else['a','b','c/k','k*z','w'])],
                            scalar=sigma,inverse_verified=True,full_ten_tensor_entries_verified=True,inverted='k',zero_complement_retained=True))
    # Canonical determinant(E,T,V) for the complete five abelian-subalgebra types.
    a,b,c=[r.atom(1,i)for i in range(3)]
    cases={'Cartan':([a,r.q.scale(b,-1),{}],[{},b,r.q.scale(c,-1)]),
           'mixed':([a,b,r.q.scale(c,-2)],[b,{},{}]),
           'regular_nilpotent':([b,c,{}],[c,{},{}]),
           'common_kernel':([c,{},{}],[{},c,{}]),
           'common_image':([b,{},{}],[c,{},{}])}
    determinant_values={}
    for name,(T,V)in cases.items():
        E=[a,b,c];value=r.q.add(*(r.q.scale(r.q.mul(E[i],r.q.mul(T[j],V[k])),sign)
                                  for i,j,k,sign in [(0,1,2,1),(1,2,0,1),(2,0,1,1),(0,2,1,-1),(2,1,0,-1),(1,0,2,-1)]))
        assert bool(value)==(name!='common_image')
        determinant_values[name]=[[list(e[:5]),str(as_constant(v))]for e,v in value.items()]
    return dict(whole_nonzero_K_normalizers=records,
                canonical_abelian_subalgebra_determinants=determinant_values,
                common_image_choice_dimension=2,ordered_independent_basis_dimension=4,
                normal_semisimple_line_dimension=2,normal_nilpotent_line_dimension=1)


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2));signal.alarm(90)
    started=time.monotonic();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    packet=json.loads((ROOT/INPUT).read_text())
    assert hashlib.sha256(packet['producer_bytes'].encode()).hexdigest()==packet['producer_sha256']=='ff542fd020159e8409b687fc0423f39f2648691b7fb40c1c77372699e3b1fd68'
    assert packet['all_passed']is False
    first,second,source=sources();A0,B0,C1,C2=explicit_zero_K_kernels(first,second)
    normalization=normalizations();records=[]
    for name in KINDS:
        saved=next(c for c in packet['cases']if c['type']==name);P=model(name)
        original=r.q.antisymmetric({ast.literal_eval(ij):r.parse(v)for ij,v in saved['normal_form_P0'].items()})
        assert all(P.get(ij,{})==original.get(ij,{})for ij in r.PAIRS)
        assert not r.divisor(P)
        M1,M2,S=r.operator(P,first),r.operator(P,second),r.shear(P,first)
        assert M1==r.matrix(saved['full_D1_350_by30'])and M2==r.matrix(saved['full_D2_350_by9'])and S==r.matrix(saved['full_actual_shear_30_by6'])
        a=saved['analysis'];certs={key:r.certificate(M,a[key])for key,M in [('D1',M1),('D2',M2),('actual_shear',S)]}
        B={'SS_COMMON_IMAGE_K1':9,'SS_COMMON_IMAGE_K0':8,'NIL_COMMON_IMAGE_K1':8,'NIL_COMMON_IMAGE_K0':7}[name]
        assert saved['affine_first_symbol_dimension']==B
        row=dict(type=name,affine_whole_symbol_dimension=B,certificates=certs)
        if name.endswith('K0'):
            for M,C in [(M1,C1),(M2,C2)]:
                assert all(not sum((a*b for a,b in zip(row,col)),r.RF())for row in M for col in zip(*C))
            assert certs['D1']['nullity']==12 and certs['D2']['nullity']==6
            assert common_direction(P)
            row['whole_source']=entire_zero_K_tensor(P,A0,B0)
            row['rank_four_G_source_empty']=True
            row['old_linear_only_image_upper']=a['PGL_image_upper_bound']
        else:
            expected=6 if name.startswith('SS')else 7
            assert certs['D1']['nullity']==expected and certs['D2']['nullity']==0
            assert B+expected+6-1==20
            row['full_projective_image_upper']=20
            row['full_D2_injective_each_P1_has_at_most_one_P2']=True
        records.append(row)
    out=dict(task_id='G',all_passed=True,independent_Fraction_only=True,producer_imports=False,CAS_imports=False,
             exact_input_bindings=job['exact_inputs_and_hashes'],producer_snapshot_sha256=packet['producer_sha256'],
             immutable_producer_linear_only_all_passed=False,full_source_quotient=source,
             whole_actual_normalizations=normalization,cases=records,
             conclusion='All four entire independent pointwise-parallel common-image globalV sources: K!=0 image<=20; K=0 every tensor rank<=2.',
             retained_complements=['dependent T=tV sources','scalar traceless-symbol-zero sources','Pi0=0 sources',
                                   'all other G branches and full E component membership'],
             elapsed_seconds=time.monotonic()-started,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUTPUT).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ('all_passed','conclusion','elapsed_seconds','peak_rss_kib')}))


if __name__=='__main__':main()
