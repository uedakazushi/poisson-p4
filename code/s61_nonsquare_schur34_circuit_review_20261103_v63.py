#!/usr/bin/env python3
"""Independent semantic audit of all entries of the exact Schur34 DAG."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PARENT='S61_BINARY_NONSQUARE_SCHUR34_CIRCUIT_20261103_v62'
UNIT='S61_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v63_REVIEW'
STEM='S61_NONSQUARE_SCHUR34_CIRCUIT_20261103_v63_REVIEW'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def review(output):
    start=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    pd=ROOT/'reports'/PARENT;parent=json.loads((pd/'report.json').read_text())
    unitpath=ROOT/'reports'/UNIT/'report.json';unit=json.loads(unitpath.read_text())
    assert parent['all_passed'] and unit['all_passed']
    for packet in (parent,unit):
        for name,digest in packet['input_sha256'].items():assert sha(ROOT/name)==digest,name
    graphpath=pd/'circuit.json.gz';graph=json.loads(gzip.decompress(graphpath.read_bytes()))
    source=graph['source_matrix'];assert source['path']==unit['actual_source']['path']
    assert source['sha256']==unit['actual_source']['sha256']==sha(ROOT/source['path'])
    with (ROOT/source['path']).open() as f:assert list(map(int,f.readline().split()))==[24,58,1]
    cp=graph['cofactor_input']
    assert cp['path']=='reports/S61_BINARY_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v62/cofactor.txt.gz'
    assert sha(ROOT/cp['path'])==cp['sha256']==unit['input_sha256'][cp['path']]
    with gzip.open(ROOT/cp['path'],'rt') as f:
        assert int(f.readline())==58;coflines=list(f);assert len(coflines)==58
    active=[j for j,line in enumerate(coflines) if int(line.split()[0])]
    assert len(active)==33 and sum(int(line.split()[0]) for line in coflines)==975
    assert graph['zero_cofactor_components']==[j for j in range(58) if j not in active]
    assert graph['source_shape']==[24,58] and graph['result_shape']==[23,58] and graph['pivot_row']==7
    rows=list(range(7))+list(range(8,24))
    assert graph['remaining_old_rows']==rows and graph['remaining_old_columns']==list(range(58))
    assert graph['product_rows_in_result']==list(range(7)) and graph['genuine_cubic_rows_in_result']==list(range(7,23))
    assert graph['field']=='QQ[a,e,P]' and graph['serialized_source_variable_meaning']=='d denotes P=ad'
    assert graph['localization']=='a*P*F*B0*C*J!=0; J=5P-6e-4; B0=aB'
    assert graph['source_entries_complete'] and graph['no_degree_or_jet_truncation']
    assert graph['pivot_unit_sparse']==[[[1,0,3],'5'],[[1,1,2],'-6'],[[1,0,2],'-4']]
    # Decode expressions bottom-up into independent formal term signatures.
    # No producer build/evaluate import, polynomial expansion, or numerical test.
    nodes=graph['nodes'];signatures=[];leaves={};qleaves={};gindices=[];edges=0
    for i,node in enumerate(nodes):
        op=node['op']
        if op=='source_entry':
            assert set(node)=={'op','row','column'}
            r,c=node['row'],node['column'];assert 0<=r<24 and 0<=c<58 and (r,c) not in leaves
            leaves[r,c]=i;signature=('A',r,c)
        elif op=='cofactor_component':
            assert set(node)=={'op','component'};j=node['component']
            assert j in active and j not in qleaves;qleaves[j]=i;signature=('q',j)
        elif op=='pivot_unit':
            assert set(node)=={'op'};gindices.append(i);signature=('g',)
        else:
            assert op in ('mul','sub','sum') and set(node)=={'op','args'}
            args=node['args'];assert all(isinstance(j,int) and 0<=j<i for j in args)
            assert len(args)==(33 if op=='sum' else 2)
            signature=(op,tuple(signatures[j] for j in args));edges+=len(args)
        signatures.append(signature)
    assert len(nodes)==6244 and len(leaves)==1392 and sorted(qleaves)==active and len(gindices)==1
    assert graph['pivot_unit_node']==gindices[0]
    image=[('sum',tuple(('mul',(('A',i,j),('q',j))) for j in active)) for i in range(24)]
    assert len(graph['image_roots'])==24
    for i,index in enumerate(graph['image_roots']):assert signatures[index]==image[i]
    assert len(graph['roots'])==23 and all(len(r)==58 for r in graph['roots'])
    for nr,i in enumerate(rows):
        for j,index in enumerate(graph['roots'][nr]):
            expected=('sub',(('mul',(('g',),('A',i,j))),('mul',(image[i],('A',7,j)))))
            assert signatures[index]==expected,(i,j)
    reachable=set();stack=[x for r in graph['roots'] for x in r]+graph['image_roots']
    while stack:
        i=stack.pop()
        if i in reachable:continue
        reachable.add(i);stack.extend(nodes[i].get('args',[]))
    assert len(reachable)==6244
    # Different full symbolic fixture; only s is inverted. Neither q component
    # e nor te-1 is a named unit. Pivot row is2, instead of author's row1.
    x,t,e,h=s.symbols('s t e h')
    A=s.Matrix([[1,0,0],[0,0,h],[x*t,-x,0]])
    q=s.Matrix([e,t*e-1,0]);v=A*q;g=x;pivot=2;keep=[0,1]
    assert s.expand(v[pivot]-g)==0
    B=s.Matrix([[s.expand(g*A[i,j]-v[i]*A[pivot,j]) for j in range(3)] for i in keep])
    assert B==s.Matrix([[x*(1-t*e),x*e,0],[0,0,x*h]]).applyfunc(s.expand)
    left=s.diag(g,g,1)
    for i in keep:left[i,pivot]=-v[i]
    right=s.eye(4)
    for j in range(3):right[3,j]=-A[pivot,j]/g
    assert s.expand(left.det()-g**2)==0 and right.det()==1
    block=(left*A.row_join(v)*right).extract([pivot]+keep,[3,0,1,2]).applyfunc(s.expand)
    expected=s.zeros(3,4);expected[0,0]=g
    for i in range(2):
        for j in range(3):expected[i+1,j+1]=B[i,j]
    assert block==expected
    minors=[s.expand(B[:,list(c)].det()) for c in itertools.combinations(range(3),2)]
    assert minors[0]==0 and s.expand(minors[1]+t*minors[2]-x**2*h)==0
    assert s.expand(A.det()-x*h)==0
    # Every nonzero minor is a multiple of s²h; the displayed combination gives
    # s²h, so this is equality of ideals, preserving the unknown h factor.
    assert all(s.cancel(m/(x**2*h)).is_polynomial(t,e) for m in minors)
    fixture=dict(source_matrix=[[str(v) for v in r] for r in A.tolist()],cofactor=[str(v) for v in q],pivot_row=pivot,
        result_matrix=[[str(v) for v in r] for r in B.tolist()],only_inverted_variable='s',
        left_determinant='s^2',right_determinant='1',source_determinant='s*h',
        result_maximal_minors=[str(v) for v in minors],minor_ideal='(s^2*h)',
        generating_combination='minor(0,2)+t*minor(1,2)=s^2*h',unknown_factor_h_preserved=True)
    files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pd/'report.json',graphpath,unitpath]
    files += [ROOT/n for packet in (parent,unit) for n in packet['input_sha256']]
    report=dict(id=STEM,all_passed=True,method='independent formal expression signatures, complete DAG references, distinct full polynomial fixture',
        actual_source_sha256=source['sha256'],source_shape=[24,58],result_shape=[23,58],all_source_leaves=1392,
        all_cofactor_leaves=33,all_topological_nodes=6244,all_edges=edges,all_output_entries=1334,all_full_image_entries=24,
        every_node_reachable=True,all_nodes_and_formulas_verified_without_producer_import=True,
        unit_gate=UNIT,unit_target='a*P^2*J',remaining_old_rows=rows,product_rows=7,genuine_cubic_rows=16,
        cokernel_equivalent_over_any_commutative_base_where_g_is_unit=True,no_unknown_core_factor_inverted=True,
        constructive_exact_circuit_pivots=34,expanded_native_cache_pivots=33,no_large_output_expansion=True,
        independent_fixture=fixture,
        not_claimed=['product41 or full cubic57 or quartic coverage','a new expanded-native34 cache','a replacement23x57 matrix'],
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',STEM,report['seconds'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();review(a.output)
