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


if __name__=='__main__':main()
