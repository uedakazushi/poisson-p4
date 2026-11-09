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


if __name__=='__main__':
    import sys
    if '--translation-nonplane' in sys.argv:translation_nonplane()
    elif '--tan-simple' in sys.argv:tan_simple()
    else:main()
