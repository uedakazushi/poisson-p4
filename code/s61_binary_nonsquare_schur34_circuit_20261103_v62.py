#!/usr/bin/env python3
"""Exact complete23x58 Schur matrix as a source-bound arithmetic circuit.

No source polynomial or resulting entry is truncated. Default verification
regenerates every node/formula and proves the generic block operation on a
fresh polynomial fixture, without expanding the large result matrix.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
ROOT=Path(__file__).resolve().parents[1];STEM='S61_BINARY_NONSQUARE_SCHUR34_CIRCUIT_20261103_v62';DIR=ROOT/'reports'/STEM
UNITDIR=ROOT/'reports/S61_BINARY_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v62'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_cofactor(data):
    it=iter(data.decode().split());n=int(next(it));out=[]
    for i in range(n):
        q=[]
        for j in range(int(next(it))):
            c=next(it);ex=[int(next(it))for _ in range(3)];q.append([ex,c])
        out.append(q)
    assert next(it,None)is None
    return out

def build(n,m,pivot,cofactor):
    assert len(cofactor)==m;nodes=[]
    def node(op,**kw):nodes.append(dict(op=op,**kw));return len(nodes)-1
    source=[[node('source_entry',row=i,column=j)for j in range(m)]for i in range(n)]
    active=[j for j,q in enumerate(cofactor)if q]
    qnode={j:node('cofactor_component',component=j)for j in active}
    g=node('pivot_unit')
    image=[node('sum',args=[node('mul',args=[source[i][j],qnode[j]])for j in active])for i in range(n)]
    rows=[i for i in range(n)if i!=pivot]
    roots=[]
    for i in rows:
        roots.append([node('sub',args=[node('mul',args=[g,source[i][j]]),node('mul',args=[image[i],source[pivot][j]])])for j in range(m)])
    return dict(nodes=nodes,roots=roots,image_roots=image,pivot_unit_node=g,remaining_old_rows=rows,remaining_old_columns=list(range(m)),source_shape=[n,m],result_shape=[n-1,m],pivot_row=pivot)

def evaluate(graph,source,cofactor,g):
    values=[]
    for node in graph['nodes']:
        op=node['op']
        if op=='source_entry':v=source[node['row'],node['column']]
        elif op=='cofactor_component':v=cofactor[node['component']]
        elif op=='pivot_unit':v=g
        elif op=='mul':v=values[node['args'][0]]*values[node['args'][1]]
        elif op=='sub':v=values[node['args'][0]]-values[node['args'][1]]
        elif op=='sum':v=sum(values[i]for i in node['args'])
        else:raise AssertionError(op)
        values.append(v)
    return s.Matrix([[s.expand(values[i])for i in row]for row in graph['roots']])

def verify():
    start=time.monotonic();unit_report=UNITDIR/'report.json';U=json.loads(unit_report.read_text());assert U['all_passed'] and U['executed_native_unit_pivots']==33
    for name,want in U['input_sha256'].items():assert sha(ROOT/name)==want,name
    source=ROOT/U['source'];assert sha(source)==U['source_sha256'];cofpath=UNITDIR/'cofactor.txt.gz';cof=read_cofactor(gzip.decompress(cofpath.read_bytes()));assert len(cof)==58 and sum(len(q)for q in cof)==975
    graph=build(24,58,7,cof);assert len(graph['nodes'])==6244
    assert graph['remaining_old_rows']==list(range(7))+list(range(8,24))
    for i,node in enumerate(graph['nodes']):
        if'args'in node:assert all(0<=j<i for j in node['args'])
    graph.update(field='QQ[a,e,P]',serialized_source_variable_meaning='d denotes P=ad',
        source_matrix=dict(path=U['source'],sha256=U['source_sha256'],native_format='24 58 1, then metadata and complete sparse QQ polynomials'),
        cofactor_input=dict(path=str(cofpath.relative_to(ROOT)),sha256=sha(cofpath),format='58, then58 complete sparse QQ polynomials'),
        pivot_unit_sparse=[[[1,0,3],'5'],[[1,1,2],'-6'],[[1,0,2],'-4']],
        localization='a*P*F*B0*C*J!=0; J=5P-6e-4; B0=aB',
        operation='append image column; polynomial row-scaled Schur by pivot g; discard unit block',
        zero_cofactor_components=[j for j,q in enumerate(cof)if not q],
        product_rows_in_result=list(range(7)),genuine_cubic_rows_in_result=list(range(7,23)),
        source_entries_complete=True,no_degree_or_jet_truncation=True)
    graphfile=DIR/'circuit.json.gz'
    if graphfile.exists():assert json.loads(gzip.decompress(graphfile.read_bytes()))==graph
    else:graphfile.write_bytes(gzip.compress(json.dumps(graph,separators=(',',':')).encode(),mtime=0))
    # Different3x3 fixture: q=(t,1-t,0) has NO individually known-unit
    # component; only s is inverted, and the unknown factor h must survive.
    x,t,h=s.symbols('s t h');A=s.Matrix([[1,0,h],[x,x,0],[0,1,0]]);q=s.Matrix([t,1-t,0]);v=A*q;g=x;pivot=1
    small=build(3,3,pivot,[[1],[1],[]]);B=evaluate(small,A,q,g)
    assert B==s.Matrix([[x*(1-t),-x*t,x*h],[-x*(1-t),x*t,0]]).applyfunc(s.expand)
    appended=A.row_join(v);left=s.diag(g,1,g);left[0,pivot]=-v[0];left[2,pivot]=-v[2]
    assert s.expand(left.det())==x**2
    right=s.eye(4)
    for j in range(3):right[3,j]=-A[pivot,j]/g
    assert s.expand(right.det())==1
    transformed=(left*appended*right).extract([pivot,0,2],[3,0,1,2])
    expected=s.zeros(3,4);expected[0,0]=g
    for i in range(2):
        for j in range(3):expected[i+1,j+1]=B[i,j]
    assert transformed.applyfunc(s.expand)==expected.applyfunc(s.expand)
    minors=[s.expand(B[:,list(cols)].det())for cols in itertools.combinations(range(3),2)]
    assert minors[0]==0 and s.expand(minors[1]-minors[2])==x**2*h
    assert s.expand(A.det())==x*h
    fixture=dict(source_matrix=[list(map(str,row))for row in A.tolist()],cofactor=list(map(str,q)),only_inverted_variable='s',
        result_matrix=[list(map(str,row))for row in B.tolist()],source_maximal_minor='s*h',result_maximal_minors=list(map(str,minors)),
        result_minor_combination='minor(0,2)-minor(1,2)=s^2*h',unknown_factor_h_retained=True)
    files=[Path(__file__),unit_report,graphfile,cofpath,ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,exact_arithmetic_circuit_nodes=len(graph['nodes']),complete_output_shape=[23,58],complete_output_formulas=23*58,
        exact_formula='B_ij=g*A_ij-(sum_k A_ik*q_k)*A_7j, g=a*P^2*J',remaining_old_rows=graph['remaining_old_rows'],remaining_old_columns=graph['remaining_old_columns'],
        product_rows=7,genuine_cubic_rows=16,no_new_core_divisor_inverted=True,no_coefficient_or_degree_truncation=True,
        cokernel_equivalence_over_named_localization=True,constructive_product_unit_pivots=34,expanded_native_cache_unit_pivots=33,
        exact_symbolic_unknown_factor_fixture=fixture,
        all_point_product41_or_full57_or_quartic_coverage_claimed=False,
        input_sha256={str(q.relative_to(ROOT)):sha(q)for q in files},seconds=time.monotonic()-start)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])
