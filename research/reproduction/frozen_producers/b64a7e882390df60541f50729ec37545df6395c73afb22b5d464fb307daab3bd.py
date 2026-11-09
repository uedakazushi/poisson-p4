#!/usr/bin/env python3
"""Independent Fraction all-source SNC confluence review; no producer/CAS.

Only previously independent sparse arithmetic is reused. The old producer is
bound through its immutable embedded/source-contract snapshot, never executed.
"""
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

from classification_two_lines_global_parent_review import (
    LaurentRing, determinant, matrix, eye, mm, skew, decode, sub, push,
    jacobi, actual_D,
)

ROOT=Path(__file__).resolve().parents[1]


def parent(R,p,q,k,u,v):
    H,A,B,C,D=[R.variable(n) for n in R.names[R.n:]]
    return skew(R,{
        (1,2):R.mul(p,R.add(R.mul(A,B),R.mul(u,H,H))),
        (3,4):R.mul(q,R.add(R.mul(C,D),R.mul(v,H,H))),
        (1,3):R.mul(k,A,C),(1,4):R.scale(R.mul(k,A,D),-1),
        (2,3):R.scale(R.mul(k,B,C),-1),(2,4):R.mul(k,B,D)})


def main():
    signal.alarm(110)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    start=time.monotonic()
    job=json.loads((ROOT/'research/classification/jobs/two_quadrics_chord_snc_separate_review.json').read_text())
    for a in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
    load=lambda p:json.loads((ROOT/p).read_text())
    saved=load('reports/CLASSIFICATION_TWO_QUADRICS_CHORD_SNC_CONFLUENCE_20261107.json')
    manifest=load('research/classification/inputs/two_quadrics_chord_snc_confluence_manifest_20261107.json')
    archive=load('research/classification/inputs/two_lines_pre_tan_simple_source_contracts_20261107.json')
    snapshot=archive['producer']
    assert hashlib.sha256(snapshot['bytes_utf8'].encode()).hexdigest()==snapshot['sha256']
    assert saved['producer_source']==snapshot['bytes_utf8']
    assert snapshot['sha256']=='2f3c3b6741f6ba8a05a9460eebc36410850c486865d0f6172ae7a19ea92411bb'
    for a in manifest['files']:
        if a['path']==snapshot['path']:
            assert a['sha256']==snapshot['sha256']
        else:
            assert hashlib.sha256((ROOT/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
    for a in archive['contracts']:
        assert hashlib.sha256(a['bytes_utf8'].encode()).hexdigest()==a['sha256'],a['path']

    S=LaurentRing(('p','q','k','u','v'),('H','A','B','C','D'))
    P=parent(S,*(S.variable(n) for n in ('p','q','k','u','v')))
    assert P==decode(S,saved['complete_parent_brackets'])
    assert not any(jacobi(S,P)[0].values())
    assert not any(jacobi(S,P)[1].values())
    assert actual_D(S,P)==S.parse('H*p*q*(A*B+u*H**2)*(C*D+v*H**2)')
    snc={n:S.constant(v) for n,v in [('p',1),('q',1),('k',0),('u',0),('v',0)]}
    assert sub(S,actual_D(S,P),snc)==S.parse('H*A*B*C*D')

    R=LaurentRing(('ep','ss','tt')+tuple(f'm{i}' for i in range(6)),('h','x','y','z','w'),invert=('ep',))
    P=parent(R,R.parse('-1+ep*ss'),R.constant(1),R.parse('-1+ep*tt'),R.constant(-1),R.constant(-1))
    G=matrix(R,saved['confluence_matrix']);Gi=matrix(R,saved['confluence_inverse'])
    assert mm(R,G,Gi)==eye(R)==mm(R,Gi,G)
    assert determinant(R,G)==R.parse('ep**2')
    bare=[[R.mul(R.parse('ep'),p) for p in row] for row in push(R,P,Gi,G)]
    assert bare==decode(R,saved['complete_pole_free_curve_before_magnetic_matching'])
    assert all(e[0]>=0 for row in bare for p in row for e in p)
    assert [[sub(R,p,{'ep':{}}) for p in row] for row in bare]==decode(R,saved['bare_endpoint'])

    # Whole six-magnetic inverse is explicit and constant, not an open chart.
    parameter_values=[R.parse(v) for v in ['-7*m4/5','-7*m1/5',
        '7*m2/10-7*m5/15','7*(m0+m4)/5-7*m2/10+7*m5/15',
        '7*m5/15','-7*m3/5']]
    assert parameter_values==[R.parse(v) for v in saved['magnetic_period_translation_parameters']]
    B=matrix(R,saved['magnetic_parameter_matrix']);Bi=matrix(R,saved['magnetic_parameter_matrix_inverse'])
    I6=[[R.constant(i==j) for j in range(6)] for i in range(6)]
    assert mm(R,B,Bi)==I6==mm(R,Bi,B)
    assert determinant(R,B)==R.constant(F(117649,93750))
    ms=[R.variable(f'm{i}') for i in range(6)]
    assert [R.add(*(R.mul(a,b) for a,b in zip(row,ms))) for row in B]==parameter_values
    ell=R.add(R.mul(parameter_values[2],R.parse('h')),R.mul(parameter_values[3],R.parse('x')),R.mul(parameter_values[4],R.parse('y')))
    kap=R.mul(parameter_values[5],R.parse('x'))
    T=eye(R)
    for i,p in [(3,ell),(4,kap)]:
        for j in range(5):T[i][j]=R.add(T[i][j],R.derivative(p,j))
    assert T==matrix(R,saved['full_translation_matrix'])
    scale=eye(R);scale[4][4]=R.constant(F(-5,7))
    frame=mm(R,scale,T);inverse=matrix(R,saved['full_target_frame_inverse'])
    assert frame==matrix(R,saved['full_target_frame'])
    assert mm(R,frame,inverse)==eye(R)==mm(R,inverse,frame)
    assert determinant(R,frame)==R.constant(F(-5,7))
    bound=[[sub(R,p,{'ss':parameter_values[0],'tt':parameter_values[1]}) for p in row] for row in bare]
    curve=push(R,bound,frame,inverse)
    assert curve==decode(R,saved['complete_actual_whole_magnetic_curve'])
    assert all(e[0]>=0 and e[1]==e[2]==0 for row in curve for p in row for e in p)
    strict,J=jacobi(R,curve)
    assert not any(strict.values()) and not any(J.values())
    endpoint=[[sub(R,p,{'ep':{}}) for p in row] for row in curve]
    source=next(c for c in load('reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json')['cases'] if c['name']=='contact31_chord')
    target=skew(R,{ij:R.parse(p) for ij,p in zip(combinations(range(1,5),2),source['brackets'])})
    assert endpoint==target
    for i,mon in enumerate(['h**2','h*x','h*y','x**2','x*y','y**2']):
        pos=R.names.index(f'm{i}')
        d={}
        for e,c in endpoint[3][4].items():
            if e[pos]:
                a=list(e);a[pos]-=1;d[tuple(a)]=e[pos]*c
        assert d==R.parse(mon)
    D=actual_D(R,curve)
    assert D==R.parse(saved['actual_full_quintic'])==R.parse(saved['expected_full_quintic'])
    Q1=R.parse('x*y-h**2')
    zold=R.add(R.parse('z'),R.scale(ell,-1))
    wold=R.add(R.parse('-7*w/5'),R.scale(kap,-1))
    Q2=R.add(R.mul(R.add(R.parse('x'),R.mul(R.parse('ep'),zold)),
                   R.add(R.parse('y+h'),R.mul(R.parse('ep'),wold))),R.parse('-h**2'))
    assert D==R.mul(R.parse('h'),Q1,Q2,R.parse('5/7+ep*m4'))
    assert actual_D(R,endpoint)==R.parse('5*h*(x*y-h**2)*(-h**2+h*x+x*y)/7')
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,case='contact31_chord',
        verification_method='Independent Fraction sparse Laurent polynomial reconstruction, no producer/S61/CAS imports.',
        immutable_producer_snapshot_sha256=snapshot['sha256'],all_archived_contract_hashes_verified=True,
        full_five_parameter_parent_strict10_and_actual_Pf_verified=True,parent_SNC_point='p=q=1,k=u=v=0; D=HABCD',
        all_actual_period_products_retained=True,confluence_det='ep^2',confluence_inverse_verified=True,
        final_curve_polynomial_at_ep0=True,only_nonconstant_localization='ep, only to interpret actual PGL for ep!=0',
        all_ten_strict_and_all350_projective_Jacobi_zero=True,whole_curve_actual_Pf_factorization_verified=True,
        full150_endpoint_equals_entire_six_magnetic_source=True,all_six_magnetic_derivatives_verified=True,
        magnetic_parameter_det='117649/93750',actual_target_frame_det='-5/7',no_magnetic_parameter_inverted=True,
        whole_source_membership='MPS40 finite closed union, via independently read irreducible-parent/SNC-meeting-component bridge',
        no_component_number_or_local_uniqueness_claim=True,other_source_excluded_from_scope='contact31_tangent_at_simple_intersection',
        input_sha256=job['exact_inputs_and_hashes'],reviewer_source=Path(__file__).read_text(),
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path=ROOT/'reports/CLASSIFICATION_TWO_QUADRICS_CHORD_SNC_SEPARATE_REVIEW_20261108.json'
    with path.open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print('PASS whole six-magnetic chord SNC:',round(out['seconds'],3),'seconds',flush=True)


def coupled_parent(R,p,q,a,u,v):
    H,A,B,C,D=[R.variable(n) for n in R.names[R.n:]]
    F0=R.add(R.mul(p,A,B),R.mul(u,H,H))
    b=R.add(p,R.scale(a,-1))
    return skew(R,{(1,2):F0,(3,4):R.add(R.mul(q,C,D),R.mul(v,F0)),
        (1,3):R.mul(a,A,C),(1,4):R.mul(b,A,D),
        (2,3):R.scale(R.mul(a,B,C),-1),(2,4):R.scale(R.mul(b,B,D),-1)})


def named_parents():
    results={}
    for label,path in [('T16','reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T16.json'),('T17','reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T17.json')]:
        R=LaurentRing(('p0','p1','p2','t0','t1'),('H','A','B','C','D'),invert=('p1',))
        source=json.loads((ROOT/path).read_text())['homogeneous_bivector']
        names=['A','B','C','D','H'];indices=[1,2,3,4,0];entries={}
        for ij,s in source.items():
            for i,name in enumerate(names):s=s.replace(f'x{i}',name)
            i,j=[indices[int(i)] for i in ij.split(',')]
            if i>j:i,j,s=j,i,'-('+s+')'
            entries[i,j]=R.parse(s)
        P=skew(R,entries)
        L=[R.parse(v) for v in (['0','0','(p1-p2)*B/2','(p1+p2)*C/4','-(p1+p2)*D/4'] if label=='T16' else ['0','0','-p1*B','-(p1-p2)*C/2','-(p1+p2)*D/2'])]
        E=[R.variable(n) for n in ('H','A','B','C','D')]
        gauged=[[R.add(P[i][j],R.scale(R.mul(E[i],L[j]),-1),R.mul(E[j],L[i])) for j in range(5)] for i in range(5)]
        if label=='T16':
            p,q,k=[R.parse(v) for v in ['(p2-p1)/2','p0+(p1+p2)/2','-(p1+p2)/4']]
            expected=parent(R,p,q,k,{},{});expected[1][2]=R.add(expected[1][2],R.parse('t0*H**2'));expected[2][1]=R.scale(expected[1][2],-1)
            expected[3][4]=R.add(expected[3][4],R.parse('t1*H**2'));expected[4][3]=R.scale(expected[3][4],-1)
            assert gauged==expected
            # Recover whole chord-helper parameters via the stated polynomial binding.
            T=LaurentRing(('p','q','k','u','v'),('H','A','B','C','D'))
            vals=['q+2*k','-p-2*k','p-2*k','p*u','q*v']
            binding={n:T.parse(v) for n,v in zip(R.names[:R.n],vals)}
            def convert(poly):
                out={}
                for e,c in poly.items():
                    term=T.constant(c)
                    for n,power in zip(R.names,e):
                        z=binding.get(n,T.variable(n) if n in T.names else {})
                        for _ in range(power):term=T.mul(term,z)
                    out=T.add(out,term)
                return out
            assert [[convert(z) for z in row] for row in gauged]==parent(T,*(T.variable(n) for n in ('p','q','k','u','v')))
        else:
            expected=coupled_parent(R,R.parse('p1'),R.parse('p0+p2'),R.parse('(p1-p2)/2'),R.parse('t0'),R.parse('t1/p1'))
            assert gauged==expected
            # p=rho1 is the sole completed-source denominator; the actual curve fixes p=1.
            T=LaurentRing(('p','q','a','u','v'),('H','A','B','C','D'),invert=('p',))
            binding={n:T.parse(v) for n,v in zip(R.names[:R.n],['q-p+2*a','p','p-2*a','u','p*v'])}
            def convert(poly):
                out={}
                for e,c in poly.items():
                    term=T.constant(c)
                    for n,power in zip(R.names,e):
                        z=binding.get(n,T.variable(n) if n in T.names else {})
                        if power<0:z=T.inverse_monomial(z)
                        for _ in range(abs(power)):term=T.mul(term,z)
                    out=T.add(out,term)
                return out
            assert [[convert(z) for z in row] for row in gauged]==coupled_parent(T,*(T.variable(n) for n in ('p','q','a','u','v')))
        assert not any(jacobi(R,gauged)[0].values())
        results[label]=dict(all_ten_brackets_Euler_gauge_verified=True,all_five_source_parameters_retained=True,
            completed_source_denominator='rho1 only' if label=='T17' else 'none',
            exact_named_parent_binding=True,no_component_ID_inferred_from_Pf_or_dimension=True)
    return results


def tan_simple():
    signal.alarm(110);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);start=time.monotonic()
    job=json.loads((ROOT/'research/classification/jobs/two_quadrics_tan_simple_snc_separate_review.json').read_text())
    for a in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
    load=lambda p:json.loads((ROOT/p).read_text())
    saved=load('reports/CLASSIFICATION_TWO_QUADRICS_TAN_SIMPLE_SNC_CONFLUENCE_20261107.json')
    assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()=='3c7d5eb884d27e5a2425336ab23d3dc866deff00b4ee8da1117e5ae6a61d091f'
    old=load('research/classification/inputs/two_quadrics_snc_pre_tan_review_source_contract_20261108.json')
    import ast
    now={n.name:hashlib.sha256(ast.dump(n,include_attributes=False).encode()).hexdigest() for n in ast.parse(Path(__file__).read_text()).body if isinstance(n,ast.FunctionDef)}
    assert all(now[k]==v for k,v in old['old_function_AST_sha256'].items())
    for a in old['files']:assert hashlib.sha256(a['bytes_utf8'].encode()).hexdigest()==a['sha256']
    S=LaurentRing(('p','q','a','u','v'),('H','A','B','C','D'))
    P=coupled_parent(S,*(S.variable(n) for n in ('p','q','a','u','v')))
    assert P==decode(S,saved['complete_parent_brackets'])
    assert not any(jacobi(S,P)[0].values()) and not any(jacobi(S,P)[1].values())
    assert actual_D(S,P)==S.parse('H*(p*A*B+u*H**2)*(q*C*D+v*(p*A*B+u*H**2))')
    assert sub(S,actual_D(S,P),{n:S.constant(v) for n,v in [('p',1),('q',1),('a',0),('u',0),('v',0)]})==S.parse('H*A*B*C*D')
    R=LaurentRing(('ep','qa','aa')+tuple(f'm{i}' for i in range(6)),('h','x','y','z','w'),invert=('ep',))
    P=coupled_parent(R,R.constant(1),R.parse('-1/3+ep*qa'),R.parse('2/3+ep*aa'),R.constant(-1),R.constant(F(-1,3)))
    G,Gi=matrix(R,saved['confluence_matrix']),matrix(R,saved['confluence_inverse'])
    assert mm(R,G,Gi)==eye(R)==mm(R,Gi,G) and determinant(R,G)==R.parse('-ep**2')
    bare=[[R.mul(R.parse('ep'),p) for p in row] for row in push(R,P,Gi,G)]
    assert bare==decode(R,saved['complete_pole_free_curve_before_magnetic_matching'])
    assert all(e[0]>=0 for row in bare for p in row for e in p)
    assert [[sub(R,p,{'ep':{}}) for p in row] for row in bare]==decode(R,saved['bare_endpoint'])
    ns=[R.parse(f'-8*m{i}/25') for i in range(6)]
    av=R.scale(R.add(ns[1],ns[5]),F(-1,3))
    qv=R.add(R.scale(ns[4],-1),R.scale(ns[3],-1),R.scale(ns[5],F(-1,3)),R.scale(ns[1],F(2,3)))
    A=R.scale(ns[0],F(3,2));B=R.add(ns[5],R.scale(ns[1],-2));C=R.add(R.scale(ns[2],3),R.scale(A,-1));E=R.scale(ns[3],F(3,2))
    vals=[qv,av,A,B,C,E]
    assert vals==[R.parse(p) for p in saved['magnetic_period_translation_parameters']]
    BM,BMi=matrix(R,saved['magnetic_parameter_matrix']),matrix(R,saved['magnetic_parameter_matrix_inverse'])
    I6=[[R.constant(i==j) for j in range(6)] for i in range(6)]
    assert mm(R,BM,BMi)==I6==mm(R,BMi,BM) and determinant(R,BM)==R.constant(F(-1769472,244140625))
    ms=[R.variable(f'm{i}') for i in range(6)]
    assert [R.add(*(R.mul(a,b) for a,b in zip(row,ms))) for row in BM]==vals
    ell=R.add(R.mul(A,R.parse('h')),R.mul(B,R.parse('x')),R.mul(C,R.parse('y')));kap=R.mul(E,R.parse('x'))
    T=eye(R)
    for i,p in [(3,ell),(4,kap)]:
        for j in range(5):T[i][j]=R.add(T[i][j],R.derivative(p,j))
    assert T==matrix(R,saved['full_translation_matrix'])
    Fm=eye(R);Fm[3][3]=R.constant(F(-5,4));Fm[4][3]=R.constant(F(25,24));Fm[4][4]=R.constant(F(5,2))
    frame=mm(R,Fm,T);inverse=matrix(R,saved['full_target_frame_inverse'])
    assert frame==matrix(R,saved['full_target_frame']) and mm(R,frame,inverse)==eye(R)==mm(R,inverse,frame)
    assert determinant(R,frame)==R.constant(F(-25,8))
    bound=[[sub(R,p,{'qa':qv,'aa':av}) for p in row] for row in bare]
    curve=push(R,bound,frame,inverse)
    assert curve==decode(R,saved['complete_actual_whole_magnetic_curve'])
    assert all(e[0]>=0 and e[1]==e[2]==0 for row in curve for p in row for e in p)
    strict,J=jacobi(R,curve);assert not any(strict.values()) and not any(J.values())
    endpoint=[[sub(R,p,{'ep':{}}) for p in row] for row in curve]
    source=next(c for c in load('reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json')['cases'] if c['name']=='contact31_tangent_at_simple_intersection')
    assert endpoint==skew(R,{ij:R.parse(p) for ij,p in zip(combinations(range(1,5),2),source['brackets'])})
    for i,mon in enumerate(['h**2','h*x','h*y','x**2','x*y','y**2']):
        pos=R.names.index(f'm{i}');d={}
        for e,c in endpoint[3][4].items():
            if e[pos]:a=list(e);a[pos]-=1;d[tuple(a)]=e[pos]*c
        assert d==R.parse(mon)
    assert actual_D(R,curve)==R.parse(saved['actual_full_quintic'])==R.parse(saved['expected_full_quintic'])
    assert actual_D(R,endpoint)==R.parse('-25*h*(h*x-y**2)*(h*x+x*y-y**2)/24')
    named=named_parents()
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,case=source['name'],method='Independent Fraction whole coupled parent and confluence, no producer/S61/CAS import.',all_old_FunctionDef_AST_preserved=True,old_source_and_job_bytes_archived=True,producer_embedded_sha256='3c7d5eb884d27e5a2425336ab23d3dc866deff00b4ee8da1117e5ae6a61d091f',whole_five_parameter_coupled_parent_strict10_Pf_and_SNC=True,original_coupled_relations='b=p-a; F=pAB+uH2; G=qCD+vF',confluence_det='-ep^2',final_curve_polynomial_at_ep0=True,full150_endpoint_equals_entire_six_source=True,whole_curve_strict10_all350_actualPf_verified=True,magnetic_parameter_det='-1769472/244140625',target_frame_det='-25/8',no_magnetic_parameter_inverted=True,zero_magnetic_boundary_included=True,whole_source_MPS40_by_same_irreducible_parent_SNC_component_bridge=True,named_parents=named,named_parent_conclusion='whole chord source inT16 and whole tan-simple source inT17, using accepted completed-source component theorem',input_sha256=job['exact_inputs_and_hashes'],reviewer_source=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,new_component_count=0)
    path=ROOT/'reports/CLASSIFICATION_TWO_QUADRICS_TAN_SIMPLE_SNC_SEPARATE_REVIEW_20261108.json'
    with path.open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print('PASS whole tan-simple + namedT16/T17:',round(out['seconds'],3),'seconds',flush=True)


if __name__=='__main__':
    import sys
    tan_simple() if '--tan-simple' in sys.argv else main()
