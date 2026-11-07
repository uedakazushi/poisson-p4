#!/usr/bin/env python3
"""Independent whole normal-nilpotent b!=0 plane review over Fraction(tau).

Only previously independent sparse/Fraction arithmetic is reused. The mutable
producer is never imported; its immutable embedded snapshot is data only.
All350 rows, all30/9 source directions, all8 minor walls are reconstructed.
"""
import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import classification_plane_rank_zero_rank_two_review as r

ROOT = Path(__file__).resolve().parents[1]
INPUT = 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json'
OUTPUT = 'reports/CLASSIFICATION_PLANE_GLOBAL_NILPOTENT_NONZERO_B_REVIEW_20261107.json'
JOB = 'research/classification/jobs/plane_global_nilpotent_nonzero_b_review.json'
KINDS = ('NIL_B1_REGULAR','NIL_B1_DIAGONAL','NIL_B1_RESONANT_B','NIL_B1_RESONANT_C')


def models(kind):
    x = [r.atom(1,i) for i in range(5)]
    a,b,c,z,w = x
    V = [b,c,{}, {},{}] if kind == KINDS[0] else [b,{}, {},{},{}]
    if kind == KINDS[0]: T = [a,{},r.q.scale(c,-1),{},{}]
    else:
        rho = r.RHO if kind == KINDS[1] else r.RF(F(-2,3) if kind == KINDS[2] else F(-1,3))
        T = [r.q.scale(a,rho+1),r.q.scale(b,rho),r.q.scale(c,-1-2*rho),{},{}]
        if kind == KINDS[2]: T[0] = r.q.add(T[0],c)
        if kind == KINDS[3]: T[2] = r.q.add(T[2],b)
    En = [{},{},{},z,w]
    Nj = [{},{},{},{},z]
    P = r.padd(r.wedge(T,En),r.wedge(V,Nj),r.q.antisymmetric({(3,4):r.q.mul(z,w)}))
    e,T,V = x[:3],T[:3],V[:3]
    determinant = r.q.add(*(r.q.scale(r.q.mul(e[i],r.q.mul(T[j],V[k])),s)
                           for i,j,k,s in [(0,1,2,1),(1,2,0,1),(2,0,1,1),
                                           (0,2,1,-1),(2,1,0,-1),(1,0,2,-1)]))
    return P,determinant


def normalizer_identities():
    """Universal polynomial Ad identities, including all centralizer parameters."""
    R = r.LaurentRing(('rho','x','y','v','s','t','u','d','k'),(),invert=('d','k'))
    p = R.parse
    def mat(rows): return [[p(v) for v in row] for row in rows]
    I = mat([['1','0','0'],['0','1','0'],['0','0','1']])
    def plus(A,B): return [[R.add(a,b) for a,b in zip(ra,rb)] for ra,rb in zip(A,B)]
    def times(A,B): return [[R.add(*(R.mul(a,b) for a,b in zip(row,col))) for col in zip(*B)] for row in A]
    def scale(A,c): return [[R.mul(a,p(c)) for a in row] for row in A]
    def comm(A,B): return plus(times(A,B),scale(times(B,A),'-1'))
    def elementary(i,j):
        A=[[{} for _ in range(3)] for _ in range(3)]; A[i][j]=p('1'); return A
    J=plus(elementary(0,1),elementary(1,2)); J2=times(J,J)
    G=mat([['1','0','0'],['0','0','0'],['0','0','-1']])
    T=plus(G,plus(scale(J,'x'),scale(J2,'y')))
    g=plus(I,plus(scale(J,'x'),scale(J2,'(y+x**2)/2')))
    gi=plus(I,plus(scale(J,'-x'),scale(J2,'(-y+x**2)/2')))
    assert times(g,gi)==I and times(times(g,T),gi)==G and times(times(g,J),gi)==J
    N,B,C=elementary(0,1),elementary(0,2),elementary(2,1)
    def family(x='x',y='y',v='v'):
        return mat([['rho+1',x,y],['0','rho','0'],['0',v,'-1-2*rho']])
    T=family()
    assert comm(T,N)==N
    records=[]
    for name,E,parameter,predicted in [
        ('B',B,'s',family('x+s*v','y-s*(2+3*rho)','v')),
        ('C',C,'t',family('x-t*y','y','v+t*(1+3*rho)')),
        ('N',N,'u',family('x-u','y','v'))]:
        g,gi=plus(I,scale(E,parameter)),plus(I,scale(E,'-'+parameter))
        assert times(g,gi)==I and times(times(g,N),gi)==N
        assert times(times(g,T),gi)==predicted
        records.append(name)
    diagonal=mat([['d','0','0'],['0','d','0'],['0','0','k']])
    inverse=mat([['1/d','0','0'],['0','1/d','0'],['0','0','1/k']])
    assert times(diagonal,inverse)==I
    assert times(times(diagonal,N),inverse)==N
    assert times(times(diagonal,B),inverse)==scale(B,'d/k')
    assert times(times(diagonal,C),inverse)==scale(C,'k/d')
    # The full tracefree solution [T,N]=N has exactly rho,x,y,v free.
    cols=[]
    for i in range(3):
        for j in range(3):
            E=elementary(i,j); A=comm(E,N)
            cols.append([a.get(R.zero,F(0)) for row in A for a in row]+[F(i==j)])
    assert r.rank(list(map(list,zip(*cols))))==5
    # [T,J]=J has tracefree affine dimension two.
    cols=[]
    for i in range(3):
        for j in range(3):
            E=elementary(i,j); A=comm(E,J)
            cols.append([a.get(R.zero,F(0)) for row in A for a in row]+[F(i==j)])
    assert r.rank(list(map(list,zip(*cols))))==7
    return dict(universal_Jordan3_exponential_and_inverse=True,
                universal_Jordan21_unipotent_Ad_and_inverse=records,
                universal_diagonal_and_inverse=True,
                full_affine_solution_dimensions=dict(Jordan3=2,Jordan21=4),
                geometric_resonances=['-2/3','-1/3'])


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(175,176))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    signal.alarm(180)
    start=time.monotonic()
    job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    packet=json.loads((ROOT/INPUT).read_text())
    assert hashlib.sha256(packet['producer_bytes'].encode()).hexdigest()==packet['producer_sha256']
    assert packet['producer_sha256']=='947b14874c9ebccc9f94164816958a57829ce3136800c7c5e6dfcd411ac5bf44'
    assert packet['all_passed'] is False  # the larger21-stratum obligation remains open
    b1,b2,source=r.source_bases(packet)
    normals=normalizer_identities()
    records=[]
    for kind in KINDS:
        saved=next(c for c in packet['cases'] if c['type']==kind)
        P,base_det=models(kind)
        parsed=r.q.antisymmetric({ast.literal_eval(ij):r.parse(v) for ij,v in saved['normal_form_P0'].items()})
        assert P==parsed and base_det==r.parse(saved['base_determinant'])
        assert r.divisor(P)==r.parse(saved['actual_D']) and base_det
        M1,M2,S=r.operator(P,b1),r.operator(P,b2),r.shear(P,b1)
        assert M1==r.matrix(saved['complete_D1_350_by30'])
        assert M2==r.matrix(saved['complete_D2_350_by9'])
        assert S==r.matrix(saved['complete_actual_shear_30_by6'])
        a=saved['analysis']
        certs={key:r.certificate(M,a[key]) for key,M in [('D1',M1),('D2',M2),('actual_shear',S)]}
        d=certs['D1']['nullity']+certs['D2']['nullity']
        B=11 if kind in KINDS[:2] else 10
        assert saved['affine_first_symbol_dimension_upper']==B
        assert d==a['higher_fiber_upper_bound'] and B+d+5==a['PGL_image_upper_bound']
        roots=[]
        if kind==KINDS[1]:
            # Product of the three independently recomputed nonzero maximal minors.
            product=r.RF(1)
            for key in ('D1','D2','actual_shear'): product=product*r.scalar(a[key]['minor_determinant'])
            assert product==r.scalar(a['exception_polynomial'])
            factored=r.RF(1); wallset=set()
            for factor,multiplicity in a['all_irreducible_factors']:
                f=r.scalar(factor)
                assert len(f.p)==2 and len(f.q)==1
                rho=-f.p[0]/f.p[1]; wallset.add(rho)
                for _ in range(multiplicity):factored=factored*f
            quotient=product/factored
            assert quotient and len(quotient.p)==len(quotient.q)==1
            expected={F(2,3),F(-5,6),F(0),F(-1),F(-1,2),F(1,3),F(-1,3),F(-2,3)}
            assert wallset==expected
            assert {F(c['parameter']) for c in a['all_exceptional_roots']}==expected
            for wall in a['all_exceptional_roots']:
                rho=F(wall['parameter'])
                checked={key:r.certificate(r.specialization(M,rho),wall[key])
                         for key,M in [('D1',M1),('D2',M2),('actual_shear',S)]}
                d=checked['D1']['nullity']+checked['D2']['nullity']
                excluded=(rho==F(-1,3))
                assert wall['excluded_multiplicity_higher']==excluded
                determinant_at={ex:r.evaluate(c,rho) for ex,c in base_det.items() if r.evaluate(c,rho)}
                assert bool(determinant_at)==(not excluded)
                assert checked['D2']['nullity']==0
                assert d==(7 if rho in (F(-1,3),F(-2,3)) else 6)
                assert d==wall['higher_fiber_upper_bound'] and 10+d+5==wall['PGL_image_upper_bound']
                roots.append(dict(parameter=str(rho),certificates=checked,
                                  m2_excluded=excluded,full_image_upper=10+d+5))
        records.append(dict(type=kind,whole_symbol_dimension_upper=B,
                            certificates=certs,all_walls=roots,
                            full_source_image_upper=22 if kind in KINDS[:2] else 21))
    out=dict(task_id='G',target='entire normal-nilpotent b!=0 mD2 global-vector source image<=22',
             all_passed=True,broader_twentyone_source_packet_all_passed=False,
             independent_imports_only=True,producer_imports=False,CAS_imports=False,
             exact_input_bindings=job['exact_inputs_and_hashes'],
             immutable_producer_sha256=packet['producer_sha256'],
             source_quotient=source,actual_normalization_identities=normals,cases=records,
             full_m2_nonzero_b_image_upper=22,
             retained_complements=['b=0 sources','mD>=3 excluded leading-Pf-zero source',
                                  'other normal-type sources','full E membership and all G exhaustion'],
             elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUTPUT).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ('all_passed','full_m2_nonzero_b_image_upper','elapsed_seconds','peak_rss_kib')}))


if __name__=='__main__':main()
