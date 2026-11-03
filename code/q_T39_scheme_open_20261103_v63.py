#!/usr/bin/env python3
"""Exact full-Jacobian minor giving a controlled reduced T39 scheme open."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from q_local_q51_quotient_20261103_v63 import evaluate,ZERO
from q_normality_open_certificate_20261103_v63 import epsilon_evaluate

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output):
    started=time.monotonic()
    parent_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    base=[Fraction(2,9)]+[Fraction(0)]*34
    for name,c in zip(['z11','z12','z23','z33'],[4,6,-4,1]):base[parent['input_names'].index(name)]=Fraction(c)
    point=[]
    for i,c in enumerate(base):
        value=list(ZERO);value[0]=c;value[i+1]=Fraction(1);point.append(tuple(value))
    jets=evaluate(parent['nodes'],None,point,True)
    assert all(not jets[i][0]for i in parent['all350_residual_nodes'])
    matrix=s.Matrix([jets[i][1:]for i in parent['all350_residual_nodes']])
    rows=list(DomainMatrix.from_Matrix(matrix.T).rref()[1])
    columns=list(DomainMatrix.from_Matrix(matrix[rows,:]).rref()[1])
    assert len(rows)==len(columns)==26
    minor=DomainMatrix.from_Matrix(matrix[rows,columns]).det()
    assert str(minor)=='-36115228490770579741605258797776896000000000000/13'
    dag=DAG(parent['input_names'])
    dag.nodes=[(op,tuple(data[0]))if op in ('+','*')else(op,data[0])for op,*data in parent['nodes']]
    dag.lookup={node:i for i,node in enumerate(dag.nodes)}
    dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
    dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
    derivatives=[]
    for index,(op,*data)in enumerate(parent['nodes']):
        if op=='q':derivative={}
        elif op=='v':derivative={data[0]:dag.one}if data[0]in columns else{}
        elif op=='+':
            active=set().union(*(derivatives[i]for i in data[0]))
            derivative={j:dag.add(*(derivatives[i][j]for i in data[0]if j in derivatives[i]))for j in active}
        elif op=='*':
            terms={}
            for position,child in enumerate(data[0]):
                others=data[0][:position]+data[0][position+1:]
                for variable,value in derivatives[child].items():
                    terms.setdefault(variable,[]).append(dag.mul(value,*others))
            derivative={j:dag.add(*terms[j])for j in terms}
        elif op=='inv':
            derivative={j:dag.neg(dag.mul(value,index,index))for j,value in derivatives[data[0]].items()}
        derivative={j:value for j,value in derivative.items()if value!=dag.zero}
        derivatives.append(derivative)
    derivative_nodes=[[derivatives[parent['all350_residual_nodes'][i]].get(j,dag.zero)for j in columns]for i in rows]
    assert json.loads(json.dumps(dag.nodes[:len(parent['nodes'])]))==parent['nodes']
    assert sum(node[0]=='inv'for node in dag.nodes)==sum(node[0]=='inv'for node in parent['nodes'])
    values=dag.evaluate(dict(zip(parent['input_names'],base)))
    assert s.Matrix([[values[i]for i in row]for row in derivative_nodes])==matrix[rows,columns]
    weights=[]
    for op,*data in dag.nodes:
        if op=='q':weight=None if not Fraction(data[0])else 0
        elif op=='v':weight=parent['input_weights'][data[0]]
        elif op=='+':
            options={weights[i]for i in data[0]if weights[i]is not None};assert len(options)<=1
            weight=next(iter(options),None)
        elif op=='*':weight=None if any(weights[i]is None for i in data[0])else sum(weights[i]for i in data[0])
        elif op=='inv':assert weights[data[0]]==0;weight=0
        weights.append(weight)
    for r,i in enumerate(rows):
        for c,j in enumerate(columns):
            assert weights[derivative_nodes[r][c]]in(None,parent['residual_weight_bounds'][i]-parent['input_weights'][j])
    degree=sum(parent['residual_weight_bounds'][i]for i in rows)-sum(parent['input_weights'][j]for j in columns)
    assert degree==264 and degree//6==44
    epsilon_inputs=[(Fraction(2,9),)]+[()for _ in range(34)]
    for name,c in zip(['z11','z12','z23','z33'],[4,6,-4,1]):
        epsilon_inputs[parent['input_names'].index(name)]=(Fraction(0),Fraction(c))
    epsilon_values=epsilon_evaluate(dag.nodes,epsilon_inputs)
    assert all(not epsilon_values[i]for i in parent['all350_residual_nodes'])
    for row in derivative_nodes:
        for index in row:
            poly=epsilon_values[index];weight=weights[index]
            assert not poly or(weight%6==0 and all(not c or 6*i==weight for i,c in enumerate(poly)))
    output.mkdir(parents=True,exist_ok=False)
    payload=dict(id='Q_T39_SCHEME_OPEN_20261103_V63',input_names=parent['input_names'],
                 input_weights=parent['input_weights'],parent_archive_sha256=sha(parent_path),
                 nodes=dag.nodes,unchanged_parent_node_prefix=len(parent['nodes']),
                 parent_all350_residual_nodes=parent['all350_residual_nodes'],
                 selected_equation_rows=rows,selected_input_columns=columns,
                 selected_derivative_matrix_nodes=derivative_nodes,
                 J_definition='det(selected_derivative_matrix_nodes)',J_relative_weight=degree)
    archive=output/'Jacobian-minor-model.json.gz';encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode()
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    certificate=dict(full350x35_Jacobian_at_Rstar_eps1=[[str(v)for v in row]for row in matrix.tolist()],
                     selected_rows=rows,selected_columns=columns,minor=str(minor))
    (output/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    note=ROOT/'research/notes/Q_T39_SCHEME_OPEN_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,
            ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',
            ROOT/'code/q_local_q51_quotient_20261103_v63.py',
            ROOT/'code/q_normality_open_certificate_20261103_v63.py',
            ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/report.json',
            ROOT/'reports/Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63/report.json',
            ROOT/'reports/Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63_REVIEW.json']
    report=dict(id=payload['id'],all_passed=True,
                full350x35_Jacobian_exact_QQ_rank_at_Rstar_eps1=26,
                exact26_minor=str(minor),parent_node_prefix_unchanged=True,
                added_derivative_nodes=len(dag.nodes)-len(parent['nodes']),no_new_inverses=True,
                selected676_symbolic_derivatives_match_independent_Fraction_jets_at_eps1=True,
                every_derivative_entry_weight_checked=True,minor_relative_weight=degree,
                exact_QQ_eps_minor_identity='J(Rstar_eps)=c_J*eps^44',
                all676_entries_checked_over_untruncated_QQ_eps=True,
                entire_nonzero_resonance_ray_full_scheme_regular_dimension=9,
                precise_scheme_open='R[1/(Delta*J)] = (R/(G51))[1/(Delta*J)] = QQ[nine parameters][1/(Delta*J)|_T39]',
                product_Delta_J_exact_ray_power=103,
                equality_on_all_Delta_open_claimed=False,full35_branch_exhaustion_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in [archive,output/'certificate.json']},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS full scheme T39 open; exact rank26, J=c*eps44;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
