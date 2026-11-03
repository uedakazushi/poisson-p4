#!/usr/bin/env python3
"""Exact T40=A9 graph quotient at a=-8/9 and an actual smooth full-scheme point.

The finite parametrization algorithm below copies the frozen Q51 algorithm,
with its center supplied as an explicit argument. No frozen file is changed.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from q_local_q51_quotient_20261103_v63 import evaluate,ZERO

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_parametrization(center,source,parent,raw,residuals,relations,weights,jac,rows,columns):
    free=[i for i in range(35)if i not in columns]
    names=[parent['input_names'][i]for i in free]
    dag=DAG(names);assigned=[dag.zero]*35;assigned[0]=dag.q(center)
    for index,name in zip(free,names):assigned[index]=dag.inputs[name]
    def clone():
        mapped=[]
        for op,*data in source.nodes:
            if op=='q':value=dag.q(data[0])
            elif op=='v':value=assigned[data[0]]
            elif op=='+':value=dag.add(*(mapped[i]for i in data[0]))
            elif op=='*':value=dag.mul(*(mapped[i]for i in data[0]))
            elif op=='inv':value=dag.inv(mapped[data[0]])
            mapped.append(value)
        return mapped
    blocks=[]
    for degree in range(1,7):
        selected_rows=[i for i in rows if weights[relations[i]]==degree]
        selected_columns=[i for i in columns if parent['input_weights'][i]==degree]
        assert len(selected_rows)==len(selected_columns)
        matrix=jac[selected_rows,selected_columns];inverse=matrix.inv();mapped=clone()
        current=[mapped[relations[i]]for i in selected_rows]
        for j,index in enumerate(selected_columns):
            assigned[index]=dag.neg(dag.add(*(dag.mul(dag.q(inverse[j,k]),current[k])for k in range(len(current)))))
        blocks.append(dict(weight=degree,rows=selected_rows,columns=selected_columns,
                           matrix=[[str(v)for v in row]for row in matrix.tolist()],
                           inverse=[[str(v)for v in row]for row in inverse.tolist()],
                           determinant=str(matrix.det())))
    assert sum(len(block['rows'])for block in blocks)==25
    final=clone()
    outputs=[*assigned,*[final[i]for i in raw],*[final[i]for i in residuals],*[final[i]for i in relations]]
    # Remove discarded intermediate-stage circuits before exact expansion/save.
    used=set(outputs);used.update(dag.inputs.values());used.update([dag.zero,dag.one]);pending=list(used)
    while pending:
        index=pending.pop();op,*data=dag.nodes[index]
        dependencies=data[0]if op in ('+','*')else[data[0]]if op=='inv'else[]
        for child in dependencies:
            if child not in used:used.add(child);pending.append(child)
    renumber={old:new for new,old in enumerate(sorted(used))};nodes=[]
    for old in sorted(used):
        op,*data=dag.nodes[old]
        nodes.append((op,tuple(renumber[i]for i in data[0]))if op in ('+','*')else
                     (op,renumber[data[0]])if op=='inv'else(op,data[0]))
    outputs=[renumber[i]for i in outputs]
    R,*variables=ring(','.join(names),s.QQ);values=[]
    for op,*data in nodes:
        if op=='q':value=R(s.Rational(data[0]))
        elif op=='v':value=variables[data[0]]
        elif op=='+':value=sum((values[i]for i in data[0]),R.zero)
        elif op=='*':
            value=R.one
            for i in data[0]:value*=values[i]
        elif op=='inv':
            assert values[data[0]].is_ground and values[data[0]]
            value=R.one/values[data[0]]
        values.append(value)
    inputs=outputs[:35];raw_nodes=outputs[35:185];jacobi_nodes=outputs[185:535];graph_nodes=outputs[535:]
    assert len(jacobi_nodes)==350 and len(graph_nodes)==45
    assert all(not values[i]for i in [*jacobi_nodes,*graph_nodes])
    for j,index in enumerate(free):assert values[inputs[index]]==variables[j]
    for j,input_index in enumerate(free):
        raw_index=parent['positive_raw_indices'][parent['free_positive_columns'][input_index-1]]
        assert values[raw_nodes[raw_index]]==variables[j]
    return dict(input_names=names,input_weights=[parent['input_weights'][i]for i in free],
                parent_input_indices=free,nodes=nodes,parent35_input_nodes=inputs,
                raw150_nodes=raw_nodes,parent350_residual_nodes=jacobi_nodes,
                Q52_graph_relation_nodes=graph_nodes,blocks=blocks,
                exact_QQ9_all350_and45_zero=True,
                physical_free_coordinate_left_inverse=True,
                total_raw150_nonzero_polynomial_terms=sum(len(values[i])for i in raw_nodes))


def build(output):
    started=time.monotonic();center=Fraction(-8,9)
    parent_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    source_path=ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/graph-model.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    source=json.loads(gzip.decompress(source_path.read_bytes()))
    assert source['parent_archive_sha256']==sha(parent_path)
    assert source['nodes'][:len(parent['nodes'])]==parent['nodes']
    weights=[];dependencies=[]
    for op,*data in source['nodes']:
        if op=='q':weight=None if not Fraction(data[0])else 0;dep=0
        elif op=='v':weight=parent['input_weights'][data[0]];dep=1<<data[0]
        elif op in ('+','*'):
            dep=0
            for i in data[0]:dep|=dependencies[i]
            if op=='+':
                options={weights[i]for i in data[0]if weights[i]is not None};assert len(options)<=1
                weight=next(iter(options),None)
            else:weight=None if any(weights[i]is None for i in data[0])else sum(weights[i]for i in data[0])
        elif op=='inv':dep=dependencies[data[0]];assert dep&~1==0;weight=0
        weights.append(weight);dependencies.append(dep)
    relations=source['Q52_graph_relation_nodes'];assert weights[relations[25]]==0 and dependencies[relations[25]]==1
    a=s.Symbol('a');K=s.QQ.frac_field(a)
    scalar=evaluate(source['nodes'],K,[K.gens[0]]+[K.zero]*34)
    assert s.cancel(K.to_sympy(scalar[relations[25]])-(9*a+8)/(2*a**2))==0
    assert s.sympify(parent['reduced_a_open']).subs(a,s.Rational(-8,9))!=0
    point=[]
    for i in range(35):
        v=list(ZERO);v[i+1]=Fraction(1)
        if i==0:v[0]=center
        point.append(tuple(v))
    jets=evaluate(source['nodes'],None,point,True)
    assert all(not jets[i][0]for i in relations)
    jac=s.Matrix([jets[i][1:]for i in relations])
    rows=list(DomainMatrix.from_Matrix(jac.T).rref()[1]);assert len(rows)==26
    columns=list(DomainMatrix.from_Matrix(jac[rows,:]).rref()[1]);assert len(columns)==26
    graph_minor=DomainMatrix.from_Matrix(jac[rows,columns]).det()
    assert str(graph_minor)=='83385908498332845/3378528649216'
    full_at_center=s.Matrix([jets[i][1:]for i in parent['all350_residual_nodes']])
    assert DomainMatrix.from_Matrix(full_at_center).rank()==1
    dag=DAG(parent['input_names']);dag.nodes=source['nodes']
    parameter=finite_parametrization(center,dag,parent,parent['raw150_reconstruction_nodes'],
                                     parent['all350_residual_nodes'],relations,weights,jac,rows,columns)
    parameter['id']='Q_T40_FINITE_QUOTIENT_20261103_V63'
    parameter['parent_archive_sha256']=sha(parent_path)
    parameter['source_graph_archive_sha256']=sha(source_path)
    nine=[Fraction(1)]*9
    point_dag=DAG(parameter['input_names']);point_dag.nodes=parameter['nodes']
    parameter_values=point_dag.evaluate(dict(zip(parameter['input_names'],nine)))
    actual=[Fraction(parameter_values[i])for i in parameter['parent35_input_nodes']]
    assert actual[0]==center
    dense_point=[]
    for i,c in enumerate(actual):
        v=list(ZERO);v[0]=c;v[i+1]=Fraction(1);dense_point.append(tuple(v))
    dense_jets=evaluate(parent['nodes'],None,dense_point,True)
    assert all(not dense_jets[i][0]for i in parent['all350_residual_nodes'])
    full_dependent=s.Matrix([[dense_jets[i][j+1]for j in columns]for i in parent['all350_residual_nodes']])
    dense_rows=list(DomainMatrix.from_Matrix(full_dependent.T).rref()[1]);assert len(dense_rows)==26
    dense_minor=DomainMatrix.from_Matrix(full_dependent[dense_rows,:]).det();assert dense_minor
    output.mkdir(parents=True,exist_ok=False)
    archive=output/'T40-nine-parameter-model.json.gz';encoded=(json.dumps(parameter,separators=(',',':'))+'\n').encode()
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    certificate=dict(center=str(center),Q52_graph_Jacobian45x35=[[str(v)for v in row]for row in jac.tolist()],
                     Q52_graph_selected_rows=rows,Q52_graph_selected_columns=columns,Q52_graph_minor=str(graph_minor),
                     nine_dense_point=list(map(str,nine)),dense_parent35_point=list(map(str,actual)),
                     dense_full350x26_dependent_Jacobian=[[str(v)for v in row]for row in full_dependent.tolist()],
                     dense_selected_equation_rows=dense_rows,dependent_parent_input_columns=columns,
                     dense_minor=str(dense_minor),triangular_blocks=parameter['blocks'])
    (output/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    note=ROOT/'research/notes/Q_T40_FINITE_QUOTIENT_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,source_path,ROOT/'code/q_local_q51_quotient_20261103_v63.py',
            ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/report.json',
            ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json',
            ROOT/'reports/Q52_UNIVERSAL_SECTION_REVIEW_20261101.json',
            ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json']
    report=dict(id=parameter['id'],all_passed=True,center=str(center),
                full_Q52_relation25='(9*a+8)/(2*a**2)',only_a_dependency_checked=True,
                graph_Q52_Jacobian_rank_at_Qa_center=26,exact_graph26_minor=str(graph_minor),
                full_parent_Jacobian_rank_at_Qa_center=1,whole_scheme_regular_at_Qa_center_claimed=False,
                actual_global_quotient='R/(G52) = QQ[nine parameters] on parent D(a)-open',
                nine_parameters=parameter['input_names'],nine_weights=parameter['input_weights'],
                triangular_block_sizes=[len(b['rows'])for b in parameter['blocks']],
                all350_and45_exact_QQ9_zero=True,physical_free_coordinate_left_inverse=True,
                raw150_polynomial_term_count=parameter['total_raw150_nonzero_polynomial_terms'],
                full_parent_Jacobian_dependent_rank_at_dense_A9_point=26,
                exact_dense26_minor=str(dense_minor),
                full_scheme_regular_at_dense_A9_point_dimension=9,
                new_branch_through_Qstar_claimed=False,full_parent_branch_exhaustion_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in [archive,output/'certificate.json']},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact T40 A9 quotient and full-scheme smooth dense point;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
