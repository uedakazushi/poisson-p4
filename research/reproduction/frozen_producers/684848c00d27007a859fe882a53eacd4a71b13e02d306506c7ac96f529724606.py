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


def translation_nonplane_source_audit(saved):
    """Reconstruct all21 coefficients, both complete charts and actual orbit47.

    The saved equations and triangular steps are data, never executed producer
    code. Every equation is first matched with the freshly expanded Jacobi.
    Ambient D96 and its seven normals are unnecessary for positive membership.
    """
    from itertools import combinations_with_replacement
    from collections import Counter
    S=LaurentRing(tuple(f'n{i}' for i in range(21)),('h','x','u','v','w'),invert=('n17',))
    p=S.parse
    raw=skew(S,{(1,2):p('x**2'),(1,3):p('h*x'),
        (2,3):p('h**2*n0+h*n1*x+h*n2*w+(1-n4)*h*u+n3*x**2+n4*v*x+n5*w*x+n6*w**2'),
        (2,4):p('h**2*n7-h*n11*u+h*n8*x+h*n9*w+n10*x**2+n11*v*x+n12*w*x+n13*w**2'),
        (3,4):p('h**2*n14+h*n15*x+h*n16*w-h*n18*u+n17*x**2+n18*v*x+n19*w*x+n20*w**2')})
    old=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json').read_text())
    original=next(z for z in old['cases'] if z['case']=='j111_rank2')
    assert raw==decode(S,original['brackets'])
    strict,projective=jacobi(S,raw)
    assert all(not v for k,v in strict.items() if k!=(2,3,4))
    coefficients={}
    for e,c in strict[2,3,4].items():
        a=e[S.n:];t=e[:S.n]+(0,)*5
        coefficients.setdefault(a,{})[t]=c
    equations=[p(v) for v in original['remaining_coefficient_equations']]
    token=lambda a:tuple(sorted(a.items()))
    assert Counter(map(token,equations))==Counter(map(token,coefficients.values()))
    records=[]
    def closed_bindings(values):
        for _ in range(24):
            new={k:sub(S,v,values) for k,v in values.items()}
            if new==values:return values
            values=new
        raise AssertionError('nontriangular binding')
    def poly_divide_monomial(a,b):return S.mul(a,S.inverse_monomial(b))
    for branch,row in enumerate(saved['cases']):
        values={'n11':S.constant(1),'n18':S.constant(branch)}
        for step in row['necessity_triangular_unit_steps']:
            expr=sub(S,equations[step['original_equation_index']],values)
            assert expr==p(step['previously_substituted_equation'])
            pivot=step['pivot'];pos=S.names.index(pivot);linear={};constant={}
            for e,c in expr.items():
                assert e[pos] in (0,1)
                if e[pos]:a=list(e);a[pos]-=1;linear[tuple(a)]=c
                else:constant[e]=c
            unit=p(step['allowed_unit'])
            ratio=poly_divide_monomial(linear,unit)
            assert set(ratio)=={S.zero} and ratio[S.zero]
            answer=poly_divide_monomial(S.scale(constant,-1),linear)
            assert answer==p(step['solution'])
            values[pivot]=answer;values=closed_bindings(values)
        assert all(not sub(S,e,values) for e in equations)
        assert all(values.get(f'n{i}',p(f'n{i}'))==p(row['complete_original_n0_to_n20_binding'][f'n{i}']) for i in range(21))
        P=[[sub(S,e,values) for e in rr] for rr in raw]
        assert P==decode(S,row['whole_universal_brackets'])
        assert not any(jacobi(S,P)[0].values()) and not any(jacobi(S,P)[1].values())
        assert actual_D(S,P)==p(row['actual_quintic'])
        if branch==0:
            partial=[[sub(S,e,{'n17':{}}) for e in rr] for rr in P]
            # Entire cubic becomes h-divisible when the only inverted unit vanishes.
            cubic=S.add(S.scale(S.mul(p('h'),partial[2][4]),-1),S.mul(p('x'),partial[3][4]))
            assert all(e[S.n]>0 for e in cubic)
        T=matrix(S,row['whole_source_cubic_normalizer']);Ti=matrix(S,row['whole_source_cubic_normalizer_inverse'])
        assert mm(S,T,Ti)==eye(S)==mm(S,Ti,T) and determinant(S,T)==S.constant(1)
        transformed=push(S,P,T,Ti)
        assert transformed[1][2]==p('x**2') and transformed[1][3]==p('h*x') and not transformed[1][4]
        assert actual_D(S,transformed)==S.mul(p('h*x'),p(row['normal_form_cubic']))
        # Exactly the previously fixed source point; no new point or anchor.
        R=LaurentRing((),('h','x','u','v','w'))
        specialized=[[sub(S,e,{k:p(v) for k,v in row['canonical_original_source_values'].items()}) for e in rr] for rr in P]
        point=[[{e[S.n:]:c for e,c in q.items()} for q in rr] for rr in specialized]
        assert point==decode(R,row['canonical_point_brackets'])
        mons=[]
        for factors in combinations_with_replacement(range(5),2):
            e=[0]*5
            for i in factors:e[i]+=1
            mons.append(tuple(e))
        mons=sorted(mons,reverse=True)
        def vector(A):return [A[i][j].get(e,F(0)) for i,j in combinations(range(5),2) for e in mons]
        xs=[R.variable(n) for n in R.names]
        vertical=[];orbit=[]
        for i in range(5):
            for j in range(5):
                V=skew(R,{(a,b):R.add(R.mul(xs[a],xs[j]) if b==i else {},R.scale(R.mul(xs[b],xs[j]),-1) if a==i else {}) for a,b in combinations(range(5),2)})
                O=skew(R,{(a,b):R.add(point[j][b] if a==i else {},point[a][j] if b==i else {},R.scale(R.mul(xs[j],R.derivative(point[a][b],i)),-1)) for a,b in combinations(range(5),2)})
                vertical.append(vector(V));orbit.append(vector(O))
        columns=vertical+orbit+[vector(point)]
        M=[list(v) for v in zip(*columns)]
        Fsaved=[[F(0) for _ in range(51+len(row['free_source_parameters']))] for _ in range(150)]
        for i,j,c in row['full_F_sparse']:Fsaved[i][j]=F(c)
        assert M==[rr[:51] for rr in Fsaved]
        def rank(rows):
            basis={}
            for rowv in rows:
                v={i:F(c) for i,c in enumerate(rowv) if c}
                for k in sorted(basis):
                    if k in v:
                        q=v[k]
                        for j,c in basis[k].items():
                            a=v.get(j,F(0))-q*c
                            if a:v[j]=a
                            else:v.pop(j,None)
                if v:
                    k=min(v);q=v[k];basis[k]={j:c/q for j,c in v.items()}
            return len(basis)
        assert rank(list(zip(*vertical)))==24 and rank(M)==47
        direction_values={k:p(v) for k,v in row['canonical_original_source_values'].items()}
        family=[]
        for n in row['free_source_parameters']:
            pos=S.names.index(n);d=[]
            for i,j in combinations(range(5),2):
                derivative={}
                for e,c in P[i][j].items():
                    if e[pos]:a=list(e);a[pos]-=1;derivative[tuple(a)]=c*e[pos]
                evaluated=sub(S,derivative,direction_values)
                d.extend(evaluated.get((0,)*S.n+e,F(0)) for e in mons)
            family.append(d)
        FM=[list(v) for v in zip(*(columns+family))]
        assert FM==Fsaved and rank(FM)==47
        records.append(dict(branch=row['branch'],whole21_reconstruction_and_all_Jacobi=True,necessity_triangular_steps=len(row['necessity_triangular_unit_steps']),only_source_unit=row['source_units'],all_source_normalizer_and_inverse=True,actual_Euler_rank=24,actual_GL_scalar_Euler_rank=47,all_real_source_rank=47,projective_orbit_dimension=22,full_ambient_D_or_normals_not_needed=True))
    return records


def translation_nonplane_parent(R,A,L,a,kind):
    h,x,u,v,w=[R.variable(n) for n in R.names[R.n:]]
    if kind=='perazzo':
        return skew(R,{(1,2):R.scale(R.mul(x,u),2),(1,3):R.mul(x,v),(1,4):R.scale(R.mul(x,w),-1),(2,3):R.mul(a,u,v),(2,4):R.add(R.scale(R.mul(L,h,v),-1),R.scale(R.mul(R.add(a,R.constant(2)),u,w),-1)),(3,4):R.add(R.mul(A,h,x),R.mul(a,v,w))})
    if kind=='essential_four':
        return skew(R,{(1,2):R.scale(R.mul(x,u),3),(1,3):R.mul(x,v),(1,4):R.scale(R.mul(x,w),-1),(2,3):R.mul(a,u,v),(2,4):R.add(R.scale(R.mul(L,v,v),-1),R.scale(R.mul(R.add(a,R.constant(3)),u,w),-1)),(3,4):R.add(R.mul(A,h,x),R.mul(a,v,w))})
    raise AssertionError(kind)


def translation_nonplane_arc(saved,kind):
    R=LaurentRing(('ep','kap'),('h','x','u','v','w'),invert=('ep',))
    p=R.parse
    amp='ampL' if kind=='perazzo' else 'ampB'
    S=LaurentRing(('ampA',amp,'period'),('H','X','U','V','W'))
    parent=translation_nonplane_parent(S,*(S.variable(n) for n in ('ampA',amp,'period')),kind)
    assert parent==decode(S,saved['complete_parent_three_parameter_brackets'])
    assert not any(jacobi(S,parent)[0].values()) and not any(jacobi(S,parent)[1].values())
    assert actual_D(S,parent)==S.parse(saved['parent_actual_quintic'])
    assert sub(S,actual_D(S,parent),{'ampA':{},amp:{},'period':S.constant(1)})==S.parse(('4' if kind=='perazzo' else '6')+'*H*X*U*V*W')
    if kind=='perazzo':
        A,L,a=p('ep**(-2)-2*kap/ep'),p('-2/ep+4*kap'),p('-2+2*kap*ep')
    elif kind=='essential_four':
        A,L,a=p('2/ep**2-6*kap/ep'),p('-6*ep'),p('-3+6*kap*ep')
    else:raise AssertionError(kind)
    assert [A,L,a]==[p(saved['actual_parent_parameter_binding'][n]) for n in ('ampA',amp,'period')]
    P=translation_nonplane_parent(R,A,L,a,kind)
    shear=matrix(R,saved['actual_shear_matrix']);si=matrix(R,saved['actual_shear_inverse'])
    expected_shear=eye(R)
    for i,j,value in [(2,0,'1/ep'),(3,1,'1/ep'),(4,0,'1/ep')]:expected_shear[i][j]=p(value)
    if kind=='perazzo':expected_shear[4][1]=p('-1/ep')
    assert shear==expected_shear
    assert mm(R,shear,si)==eye(R)==mm(R,si,shear) and determinant(R,shear)==p('1')
    internal=[[R.mul(p('ep'),f) for f in row] for row in push(R,P,si,shear)]
    assert internal==decode(R,saved['complete_pole_free_internal_curve'])
    G=matrix(R,saved['constant_target_PGL']);Gi=matrix(R,saved['constant_target_inverse_PGL'])
    assert mm(R,G,Gi)==eye(R)==mm(R,Gi,G) and determinant(R,G)==p(saved['constant_target_determinant'])
    curve=push(R,internal,G,Gi)
    assert curve==decode(R,saved['complete_actual_regular_kappa_curve'])
    assert all(e[0]>=0 for row in curve for f in row for e in f)
    assert not any(jacobi(R,curve)[0].values()) and not any(jacobi(R,curve)[1].values())
    endpoint=[[sub(R,f,{'ep':{}}) for f in row] for row in curve]
    assert endpoint==decode(R,saved['complete_endpoint'])
    source=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_COMPLETION_20261108.json').read_text())
    row=next(a for a in source['cases'] if a['branch']==('L_h_minus_x' if kind=='perazzo' else 'L_h'))
    assert endpoint==decode(R,row['canonical_actual_kappa_model'])
    fixed=matrix(R,row['canonical_cubic_normalizer'])
    inverse_fixed=eye(R)
    if kind=='essential_four':inverse_fixed[3][4]=p('2')
    assert mm(R,fixed,inverse_fixed)==eye(R)==mm(R,inverse_fixed,fixed)
    point_recovered=push(R,[[sub(R,e,{'kap':p('1')}) for e in rr] for rr in endpoint],inverse_fixed,fixed)
    assert point_recovered==decode(R,row['canonical_point_brackets'])
    assert actual_D(R,curve)==p(saved['actual_full_quintic'])
    assert actual_D(R,endpoint)==p(saved['endpoint_quintic'])
    return dict(branch=row['branch'],complete_parent_three_parameters_and_SNC=True,whole_ep_kappa_full150_actual_curve=True,full350_and_strict10_zero=True,all_actual_Pf_quintic_coefficients=True,all_shear_and_target_inverses=True,all_negative_ep_powers_cancelled=True,kappa_not_inverted=True,kappa_zero_retained=True,same_predeclared_canonical_model=True)


def translation_nonplane():
    signal.alarm(110)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    start=time.monotonic()
    job=json.loads((ROOT/'research/classification/jobs/two_hyperplanes_translation_nonplane_separate_review.json').read_text())
    for a in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
    snapshot=json.loads((ROOT/'research/classification/inputs/two_quadrics_snc_pre_translation_review_source_contracts_20261108.json').read_text())
    import ast
    current=ast.parse(Path(__file__).read_text())
    old=ast.parse(snapshot['reviewer']['bytes_utf8'])
    olddefs={n.name:ast.dump(n,include_attributes=False) for n in old.body if isinstance(n,ast.FunctionDef)}
    newdefs={n.name:ast.dump(n,include_attributes=False) for n in current.body if isinstance(n,ast.FunctionDef)}
    assert all(newdefs[k]==v for k,v in olddefs.items())
    assert hashlib.sha256(snapshot['reviewer']['bytes_utf8'].encode()).hexdigest()==snapshot['reviewer']['sha256']
    for a in snapshot['contracts']:assert hashlib.sha256(a['bytes_utf8'].encode()).hexdigest()==a['sha256']
    saved=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_COMPLETION_20261108.json').read_text())
    producer='code/classification_two_hyperplanes_cubic_nilpotent.py'
    assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['input_sha256'][producer]
    source=translation_nonplane_source_audit(saved)
    arcs=[]
    for kind,name in [('perazzo','PERAZZO'),('essential_four','ESSENTIAL_FOUR')]:
        data=json.loads((ROOT/f'reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_{name}_SS_CONFLUENCE_20261108.json').read_text())
        assert hashlib.sha256(data['producer_source'].encode()).hexdigest()==data['input_sha256'][producer]
        arcs.append(translation_nonplane_arc(data,kind))
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,method='Independent Fraction whole-source/orbit/arcs; no producer, S61 engine or CAS import.',whole_source_checks=source,actual_arcs=arcs,old_FunctionDef_AST_and_contracts_preserved=True,geometric_bridge='Independent13flag+9fixedfibre upper and actual22orbit lower imply whole irreducible source saturation equals canonical22orbit closure; parent SNC-component theorem covers whole arcs and hence whole sources.',ambient_componenthood_or_local_uniqueness_not_claimed=True,new_component_count=0,mathematical_input_bindings=job['exact_inputs_and_hashes'],reviewer_source=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path=ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_SEPARATE_REVIEW_20261108.json'
    with path.open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print('PASS translation nonplane whole sources and arcs:',round(out['seconds'],3),'seconds',flush=True)


def j21_hx_parent(R,A,b,a):
    h,x,u,v,w=[R.variable(n) for n in R.names[R.n:]]
    return skew(R,{(1,2):R.scale(R.mul(x,u),3),(1,3):R.mul(x,v),(1,4):R.scale(R.mul(x,w),-1),
        (2,3):R.mul(a,u,v),
        (2,4):R.add(R.scale(R.mul(b,v,v),-1),R.scale(R.mul(a,u,w),-1)),
        (3,4):R.add(R.mul(A,h,h),R.mul(a,v,w))})


def j21_hx():
    """One complete six-source inverse normalization and one full(s,k) arc."""
    from collections import Counter
    import ast
    signal.alarm(110)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    start=time.monotonic()
    job=json.loads((ROOT/'research/classification/jobs/two_hyperplanes_j21_h_x_ss_separate_review.json').read_text())
    for entry in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest()==entry['sha256'],entry['path']
    old=json.loads((ROOT/'research/classification/inputs/two_quadrics_snc_pre_j21_hx_review_source_contracts_20261108.json').read_text())
    olddefs={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(old['reviewer']['bytes_utf8']).body if isinstance(n,ast.FunctionDef)}
    defs={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(__file__).read_text()).body if isinstance(n,ast.FunctionDef)}
    assert all(defs[k]==v for k,v in olddefs.items())
    assert hashlib.sha256(old['reviewer']['bytes_utf8'].encode()).hexdigest()==old['reviewer']['sha256']
    for a in old['contracts']:assert hashlib.sha256(a['bytes_utf8'].encode()).hexdigest()==a['sha256']
    saved=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_J21_H_X_SS_CONFLUENCE_20261108.json').read_text())
    assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['input_sha256']['code/classification_two_hyperplanes_cubic_nilpotent.py']
    prior=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json').read_text())
    row=next(a for a in prior['cases'] if a['case']=='j21_h_x')
    data=json.loads((ROOT/'code/data/classification_two_hyperplanes_cubic_nilpotent_sources.json').read_text())
    branch=next(a for a in data['branches'] if a['name']=='j21_h_x_irreducible')
    S=LaurentRing(tuple(f'n{i}' for i in range(14)),('h','x','u','v','w'))
    p=S.parse
    raw=skew(S,{(1,2):p('v*x'),(1,3):p('h*x'),(1,4):p('x**2'),
        (2,3):p('h**2*n0+h*n1*x-h*n3*w-2*h*n4*u+n2*x**2+n3*v*x+n4*v**2'),
        (2,4):p('h**2*n5+h*n10*v+h*n11*w+h*n6*x-h*n8*w-2*h*n9*u+n12*w*x+2*n13*u*x-n13*v*w+n7*x**2+n8*v*x+n9*v**2-v*w'),
        (3,4):p('h**2*n10+h*n11*x-h*n13*w-h*w+n12*x**2+n13*v*x')})
    assert raw==decode(S,row['brackets'])
    sj,_=jacobi(S,raw)
    assert all(not v for k,v in sj.items() if k!=(2,3,4))
    co={}
    for e,c in sj[2,3,4].items():co.setdefault(e[S.n:],{})[e[:S.n]+(0,)*5]=c
    token=lambda a:tuple(sorted(a.items()))
    assert Counter(map(token,co.values()))==Counter(token(p(v)) for v in row['remaining_coefficient_equations'])
    values={n:p(v) for n,v in saved['complete_original_binding'].items()}
    assert values=={n:p(v) for n,v in branch['substitution'].items()}
    P=[[sub(S,f,values) for f in rr] for rr in raw]
    assert P==decode(S,saved['complete_original_full_source'])
    assert not any(jacobi(S,P)[0].values()) and not any(jacobi(S,P)[1].values())
    T=eye(S)
    A=p('n12');B=p('n12*n8+n7');C=p('n0+n10*n12')
    T[2][0]=S.scale(C,-1);T[2][1]=S.scale(B,-1);T[2][4]=S.scale(A,-1)
    T[3][1]=S.scale(A,-1);T[4][0]=p('-n10');T[4][1]=p('-n8')
    Ti=matrix(S,saved['actual_source_inverse_normalizer'])
    assert T==matrix(S,saved['actual_source_normalizer'])
    assert mm(S,T,Ti)==eye(S)==mm(S,Ti,T) and determinant(S,T)==p('1')
    model=skew(S,{(1,2):p('x*v'),(1,3):p('h*x'),(1,4):p('x**2'),
        (2,3):p('v**2/2-h*u'),(2,4):p('n5*h**2+n8*x*v-v*w'),(3,4):p('n8*h*x-h*w')})
    assert push(S,P,T,Ti)==model and push(S,model,Ti,T)==P
    assert actual_D(S,model)==p('h*x*(-n5*h**3+x*(v**2/2-h*u))')
    # The normalizer preserves V0=v*d_u+h*d_v+x*d_w, in every parameter.
    cofactor=[{}, {},p('v'),p('h'),p('x')]
    moved=[S.add(*(S.mul(a,b) for a,b in zip(rr,cofactor))) for rr in T]
    xs=[S.variable(n) for n in S.names[S.n:]]
    invcoords={n:S.add(*(S.mul(a,b) for a,b in zip(rr,xs))) for n,rr in zip(S.names[S.n:],Ti)}
    assert [sub(S,v,invcoords) for v in moved]==cofactor

    Q=LaurentRing(('ampA','ampB','period'),('H','X','U','V','W'))
    parent=j21_hx_parent(Q,*(Q.variable(n) for n in ('ampA','ampB','period')))
    assert parent==decode(Q,saved['whole_parent_three_parameter_brackets'])
    assert not any(jacobi(Q,parent)[0].values()) and not any(jacobi(Q,parent)[1].values())
    assert actual_D(Q,parent)==Q.parse('H*X*(3*ampA*H**2*U+ampB*V**3+3*period*U*V*W)')
    assert sub(Q,actual_D(Q,parent),{'ampA':{},'ampB':{},'period':Q.constant(1)})==Q.parse('3*H*X*U*V*W')

    R=LaurentRing(('ep','amp','kap'),('h','x','u','v','w'),invert=('ep',))
    p=R.parse
    a=p('1-3*kap*ep');ampA=R.mul(a,p('1/(2*ep**3)'));ampB=p('-8*amp*ep**4')
    assert [a,ampA,ampB]==[p(saved['actual_parent_parameter_binding'][n]) for n in ('period','ampA','ampB')]
    P=j21_hx_parent(R,ampA,ampB,a)
    G=matrix(R,[['1','0','0','0','0'],['0','1','0','0','0'],['0','1/(3*ep)','0','0','1'],['1/(2*ep**2)','0','1','1/(2*ep)','0'],['-1/ep','0','0','1','0']])
    Gi=matrix(R,saved['actual_inverse_coordinate_map'])
    assert G==matrix(R,saved['actual_coordinate_map'])
    assert mm(R,G,Gi)==eye(R)==mm(R,Gi,G) and determinant(R,G)==p('1')
    curve=[[R.mul(p('ep'),f) for f in rr] for rr in push(R,P,Gi,G)]
    assert curve==decode(R,saved['complete_pole_free_actual_curve'])
    assert all(all(a>=0 for a in e[:3]) for rr in curve for f in rr for e in f)
    assert not any(jacobi(R,curve)[0].values()) and not any(jacobi(R,curve)[1].values())
    endpoint=[[sub(R,f,{'ep':{}}) for f in rr] for rr in curve]
    assert endpoint==decode(R,saved['complete_endpoint'])==skew(R,{(1,2):p('x*v'),(1,3):p('h*x'),(1,4):p('x**2'),(2,3):p('v**2/2-h*u'),(2,4):p('amp*h**2+kap*x*v-v*w'),(3,4):p('kap*h*x-h*w')})
    assert actual_D(R,curve)==p(saved['actual_full_curve_quintic'])
    xs=[R.variable(n) for n in R.names[R.n:]]
    oldcoords={n:R.add(*(R.mul(a,b) for a,b in zip(rr,xs))) for n,rr in zip(R.names[R.n:],G)}
    assert actual_D(R,curve)==R.mul(p('ep**2'),sub(R,actual_D(R,P),oldcoords))
    assert actual_D(R,endpoint)==p('h*x*(-amp*h**3+x*(v**2/2-h*u))')
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,whole_source='j21_h_x_irreducible all six source parameters',method='Independent Fraction one whole(s,k) actualSS arc; no producer, S61 or CAS imported.',original14source_and_all_Jacobi_reconstructed=True,original_irreducible_six_binding_matches_accepted_source=True,all_six_forward_inverse_full150=True,source_cofactor_preserved_exactly=True,source_parameter_inversions=[],source_normalizer_determinant='1',parent_full_three_parameter_strict10_all350_Pf_and_SNC=True,actual_parent_parameter_product_A_equals_a_over_2ep3_retained=True,full_ep_amp_kap_curve_all150_polynomial=True,whole_curve_strict10_all350_actualPf_and_parent_pullback=True,all_s_and_k_zero_boundaries_retained=True,all_actual_PGL_inverses_and_determinant=True,full_original_six_source_MPS_by_irreducible_parent_SNC_and_closedness=True,old_FunctionDef_AST_and_current_contract_bytes_preserved=True,new_point_anchor_or_normal_jet=False,new_component_count=0,ambient_componenthood_and_local_uniqueness_not_claimed=True,input_sha256=job['exact_inputs_and_hashes'],reviewer_source=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path=ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_J21_H_X_SS_SEPARATE_REVIEW_20261108.json'
    with path.open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print('PASS j21(h,x) whole six source actualSS:',round(out['seconds'],3),'seconds',flush=True)


def j21_zero_hx():
    """Whole normalized n6!=0 source; polynomial modular twist and SS open."""
    from collections import Counter
    import ast
    signal.alarm(110)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    start=time.monotonic()
    job=json.loads((ROOT/'research/classification/jobs/two_hyperplanes_j21_zero_hx_modular_separate_review.json').read_text())
    for entry in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest()==entry['sha256'],entry['path']
    old=json.loads((ROOT/'research/classification/inputs/two_quadrics_snc_pre_j21_zero_hx_review_source_contracts_20261108.json').read_text())
    defs=lambda src:{n.name:ast.dump(n,include_attributes=False) for n in ast.parse(src).body if isinstance(n,ast.FunctionDef)}
    now=defs(Path(__file__).read_text())
    assert all(now[k]==v for k,v in defs(old['reviewer']['bytes_utf8']).items())
    assert hashlib.sha256(old['reviewer']['bytes_utf8'].encode()).hexdigest()==old['reviewer']['sha256']
    for entry in old['contracts']:assert hashlib.sha256(entry['bytes_utf8'].encode()).hexdigest()==entry['sha256']
    saved=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_CONFLUENCE_20261108.json').read_text())
    assert hashlib.sha256(saved['producer_source'].encode()).hexdigest()==saved['input_sha256']['code/classification_two_hyperplanes_cubic_nilpotent.py']
    prior=json.loads((ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json').read_text())
    row=next(a for a in prior['cases'] if a['case']=='j21_0_hx')
    S=LaurentRing(tuple(f'n{i}' for i in range(19))+('amplitude',),('h','x','u','v','w'),invert=('amplitude',))
    p=S.parse
    raw=skew(S,{(1,2):p('v*x'),(1,4):p('h*x+x**2'),
        (2,3):p('h**2*n0+h*n1*x+h*n2*v-h*n6*u+n3*x**2+n4*v*x+n5*v**2-n6*u*x+n6*v*w'),
        (2,4):p('h**2*n7-h*n13*u+h*n14*w+h*n15*u-h*n16*w+h*n8*x+h*n9*v+n10*x**2+n11*v*x+n12*v**2-n13*u*x+n13*v*w+n16*w*x+n17*u*x+n18*u*v-u*x'),
        (3,4):p('h**2*n14-h**2*n16+h*n14*x+h*n15*v+n16*x**2+n17*v*x+n18*v**2')})
    assert raw==decode(S,row['brackets'])
    strict,_=jacobi(S,raw)
    assert all(not v for k,v in strict.items() if k!=(2,3,4))
    co={}
    for e,c in strict[2,3,4].items():co.setdefault(e[S.n:],{})[e[:S.n]+(0,)*5]=c
    token=lambda a:tuple(sorted(a.items()))
    assert Counter(map(token,co.values()))==Counter(token(p(v)) for v in row['remaining_coefficient_equations'])
    scale=eye(S);inverse=eye(S)
    for i in (2,3):scale[i][i]=p('amplitude');inverse[i][i]=p('1/amplitude')
    scaled=push(S,raw,scale,inverse)
    assert mm(S,scale,inverse)==eye(S)==mm(S,inverse,scale)
    cubic=S.add(S.mul(p('v'),scaled[3][4]),S.mul(p('h+x'),scaled[2][3]))
    assert S.derivative(cubic,2)==p('-amplitude*n6*(h+x)**2')
    assert scaled[1][2]==p('x*v') and not scaled[1][3] and scaled[1][4]==p('x*(h+x)')
    values={n:p(v) for n,v in saved['complete_original_binding'].items()}
    P=[[sub(S,f,values) for f in rr] for rr in raw]
    assert P==decode(S,saved['complete_original_source'])
    assert not any(jacobi(S,P)[0].values()) and not any(jacobi(S,P)[1].values())
    T=eye(S);T[2][0]=p('-n0');T[2][1]=p('-n3');T[4][0]=p('n16+n2');T[4][1]=p('n16+n4');T[4][3]=p('n5')
    Ti=matrix(S,saved['actual_source_inverse'])
    assert T==matrix(S,saved['actual_source_normalizer'])
    assert mm(S,T,Ti)==eye(S)==mm(S,Ti,T) and determinant(S,T)==p('1')
    model=skew(S,{(1,2):p('x*v'),(1,4):p('x*(h+x)'),
        (2,3):p('-(h+x)*u+v*w-n16*(h+x)*v'),
        (2,4):p('(h+x)*u-v*w+n16*(h+x)*w+2*n16*x*v'),
        (3,4):p('x*v+n16*(h+x)**2')})
    assert push(S,P,T,Ti)==model and push(S,model,Ti,T)==P
    assert actual_D(S,model)==p('h*x*(-(h+x)**2*u+(h+x)*v*w+x*v**2)')
    cofactor=[{}, {},p('v'),{},p('h+x')];xs=[S.variable(n) for n in S.names[S.n:]]
    invcoords={n:S.add(*(S.mul(a,b) for a,b in zip(rr,xs))) for n,rr in zip(S.names[S.n:],Ti)}
    moved=[S.add(*(S.mul(a,b) for a,b in zip(rr,cofactor))) for rr in T]
    assert [sub(S,v,invcoords) for v in moved]==cofactor

    R=LaurentRing(('ep','kap'),('h','x','u','v','w'),invert=('ep',));p=R.parse
    Q=skew(R,{(1,2):p('x*v'),(1,3):p('ep*x*v'),(1,4):p('x*(h+x+ep*w)'),
        (2,3):p('-(h+x+ep*w)*u+v*w'),(2,4):p('(h+x+ep*w)*u-v*w'),(3,4):p('x*v')})
    assert Q==decode(R,saved['actual_base_curve'])
    assert not any(jacobi(R,Q)[0].values()) and not any(jacobi(R,Q)[1].values())
    V=[{}, {},p('v'),p('ep*v'),p('h+x+ep*w')]
    Z=[R.add(*(R.derivative(Q[i][j],i) for i in range(1,5))) for j in range(5)]
    assert V==[p(v) for v in saved['Hamiltonian_field']]
    assert Z==[p(v) for v in saved['actual_affine_modular_field']]
    assert R.add(*(R.derivative(V[i],i) for i in range(5)))==p('2*ep')
    assert all(not R.add(*(R.add(R.mul(V[j],R.derivative(Z[i],j)),R.scale(R.mul(Z[j],R.derivative(V[i],j)),-1)) for j in range(5))) for i in range(5))
    def Lie(field):
        return [R.add(*(R.add(R.mul(field[l],R.derivative(Q[i][j],l)),R.scale(R.mul(Q[l][j],R.derivative(field[i],l)),-1),R.scale(R.mul(Q[i][l],R.derivative(field[j],l)),-1)) for l in range(5))) for i,j in combinations(range(5),2)]
    assert not any(Lie(V)) and not any(Lie(Z))
    curve=skew(R,{(i,j):R.add(Q[i][j],R.mul(p('kap'),R.add(R.mul(V[i],Z[j]),R.scale(R.mul(V[j],Z[i]),-1)))) for i,j in combinations(range(5),2)})
    assert curve==decode(R,saved['complete_actual_polynomial_kappa_epsilon_curve'])
    assert all(e[0]>=0 and e[1]>=0 for rr in curve for f in rr for e in f)
    assert not any(jacobi(R,curve)[0].values()) and not any(jacobi(R,curve)[1].values())
    Dbase=p('h*x*(x*v**2+(h+x+ep*w+ep*v)*(-(h+x+ep*w)*u+v*w))')
    assert actual_D(R,Q)==Dbase and actual_D(R,curve)==R.mul(p('1+2*kap*ep'),Dbase)==p(saved['actual_full_curve_quintic'])
    endpoint=[[sub(R,f,{'ep':{}}) for f in rr] for rr in curve]
    assert endpoint==decode(R,saved['complete_canonical_kappa_endpoint'])
    assert actual_D(R,endpoint)==p(saved['endpoint_actual_quintic'])
    G=matrix(R,[['1','0','0','0','0'],['0','1','0','0','0'],['0','0','1','-1/ep','0'],['0','0','0','1','0'],['1/ep','1/ep','0','0','1']])
    Gi=matrix(R,saved['complete_actual_semisimple_inverse'])
    assert G==matrix(R,saved['complete_actual_semisimple_coordinate_map'])
    assert determinant(R,G)==p('1') and mm(R,G,Gi)==eye(R)==mm(R,Gi,G)
    semi=push(R,curve,G,Gi)
    assert push(R,semi,Gi,G)==curve
    scalar=p('ep*(1+2*ep*kap)')
    assert not semi[1][2] and semi[1][3]==R.mul(scalar,p('x*v')) and semi[1][4]==R.mul(scalar,p('x*w'))
    # Complete normalized150 entries have only the declared e*c denominators.
    denominators=[]
    def rational_node(node):
        one=R.constant(1)
        if isinstance(node,ast.Constant) and type(node.value)is int:return R.constant(node.value),one
        if isinstance(node,ast.Name):return R.variable(node.id),one
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
            a,b=rational_node(node.operand);return R.scale(a,-1),b
        assert isinstance(node,ast.BinOp)
        a,b=rational_node(node.left);c,d=rational_node(node.right)
        if isinstance(node.op,ast.Add):return R.add(R.mul(a,d),R.mul(b,c)),R.mul(b,d)
        if isinstance(node.op,ast.Sub):return R.add(R.mul(a,d),R.scale(R.mul(b,c),-1)),R.mul(b,d)
        if isinstance(node.op,ast.Mult):return R.mul(a,c),R.mul(b,d)
        if isinstance(node.op,ast.Div):denominators.append(c);return R.mul(a,d),R.mul(b,c)
        assert isinstance(node.op,ast.Pow) and set(c)<= {R.zero} and d==one
        k=int(c.get(R.zero,0));assert F(k)==c.get(R.zero,F(0))
        if k<0:a,b=b,a;denominators.append(b)
        return R.mul(*([a]*abs(k))),R.mul(*([b]*abs(k)))
    for key,text in saved['complete_semisimple_tensor'].items():
        i,j=map(int,key.split(','));num,den=rational_node(ast.parse(text,mode='eval').body)
        assert R.mul(semi[i][j],den)==R.mul(scalar,num)
    present={tuple(map(int,k.split(','))) for k in saved['complete_semisimple_tensor']}
    assert all(not semi[i][j] for i,j in combinations(range(5),2) if(i,j)not in present)
    allowed=[p(f'ep**{m}*(1+2*ep*kap)**{n}') for m in range(7) for n in range(4)]
    def allowed_denominator(den):
        for a in allowed:
            if set(den)==set(a):
                e0=next(iter(a))
                if all(den[e]*a[e0]==a[e]*den[e0] for e in a):return True
        return False
    for den in denominators:
        assert den and allowed_denominator(den)
    # Direct eigenfield, including its zero-eigenspace translation complement.
    xs=[R.variable(n) for n in R.names[R.n:]]
    invcoords={n:R.add(*(R.mul(a,b) for a,b in zip(rr,xs))) for n,rr in zip(R.names[R.n:],Gi)}
    moved=[sub(R,R.add(*(R.mul(a,b) for a,b in zip(rr,V))),invcoords) for rr in G]
    assert moved==[{}, {},{},p('ep*v'),p('ep*w')]
    out=dict(task_id='E_reduced_nonnormal',all_passed=True,method='Independent Fraction whole six-source modular arc; no producer, S61 or CAS import.',source_branch='j21_0_hx n6!=0',original19source_and_all17_Jacobi_coefficients_reconstructed=True,n6_amplitude_normalizer_whole_chart_verified=True,n6_zero_complement_retained=True,accepted_six_binding_universal_Jacobi_and_full150_inverse=True,source_cofactor_exactly_preserved=True,source_normalizer_parameter_inversions=[],full_epsilon_Q_and_modular_twist_kappa_polynomial=True,Hamiltonian_modular_Lie_Poisson_commuting_identities=True,whole_curve_strict10_original350_actualPf=True,complete_semisimple150_rational_forward_inverse_verified=True,only_semisimple_open='ep*(1+2kap*ep)!=0',zero_eigen_translation_is_zero=True,all_kappa_endpoint_and_all_six_source_closure_MPS=True,semisimple_application_requires_full_existing_diag011_support=True,no_new_anchor_point_normal_jet=True,new_component_count=0,old_FunctionDef_AST_and_exact_current_contracts_preserved=True,input_sha256=job['exact_inputs_and_hashes'],reviewer_source=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    path=ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_SEPARATE_REVIEW_20261108.json'
    with path.open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print('PASS j21(0,h+x) whole six modularSS:',round(out['seconds'],3),'seconds',flush=True)


if __name__=='__main__':
    import sys
    if '--j21-zero-hx' in sys.argv:j21_zero_hx()
    elif '--j21-hx' in sys.argv:j21_hx()
    elif '--translation-nonplane' in sys.argv:translation_nonplane()
    elif '--tan-simple' in sys.argv:tan_simple()
    else:main()
