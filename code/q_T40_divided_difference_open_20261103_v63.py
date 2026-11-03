#!/usr/bin/env python3
"""A finite F=M*N certificate for an explicit full T40 scheme open."""
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
from q_normality_open_certificate_20261103_v63 import epsilon_evaluate

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def ancestors(nodes,outputs):
    used=set(outputs);pending=list(used)
    while pending:
        op,*data=nodes[pending.pop()]
        previous=data[0]if op in ('+','*')else[data[0]]if op=='inv'else[]
        for i in previous:
            if i not in used:used.add(i);pending.append(i)
    return used


def polynomial_values(nodes,needed,names,input_map):
    R,*variables=ring(','.join(names),s.QQ);values={}
    for index in sorted(needed):
        op,*data=nodes[index]
        if op=='q':value=R(s.Rational(data[0]))
        elif op=='v':
            assert data[0]in input_map
            value=variables[input_map[data[0]]]
        elif op=='+':value=sum((values[i]for i in data[0]),R.zero)
        elif op=='*':
            value=R.one
            for i in data[0]:value*=values[i]
        elif op=='inv':
            assert values[data[0]].is_ground and values[data[0]]
            value=R.one/values[data[0]]
        values[index]=value
    return values


def build(output):
    started=time.monotonic()
    parent_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    parameter_path=ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/T40-nine-parameter-model.json.gz'
    jacobian_path=ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/certificate.json'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    parameter=json.loads(gzip.decompress(parameter_path.read_bytes()))
    jacobian=json.loads(jacobian_path.read_text());rows=jacobian['dense_selected_equation_rows'];columns=jacobian['dependent_parent_input_columns']
    assert len(rows)==len(columns)==26
    assert columns==[i for i in range(35)if i not in parameter['parent_input_indices']]
    # Replay every saved polynomial identity defining the actual A9 quotient.
    saved=polynomial_values(parameter['nodes'],set(range(len(parameter['nodes']))),
                            parameter['input_names'],{i:i for i in range(9)})
    assert all(not saved[i]for i in [*parameter['parent350_residual_nodes'],*parameter['Q52_graph_relation_nodes']])
    dag=DAG(parent['input_names'])
    dag.nodes=[(op,tuple(data[0]))if op in ('+','*')else(op,data[0])for op,*data in parent['nodes']]
    dag.lookup={node:i for i,node in enumerate(dag.nodes)}
    dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
    dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
    parameter_needed=ancestors(parameter['nodes'],parameter['parent35_input_nodes'])
    parameter_images={}
    for index in sorted(parameter_needed):
        op,*data=parameter['nodes'][index]
        if op=='q':value=dag.q(data[0])
        elif op=='v':value=dag.inputs[parameter['input_names'][data[0]]]
        elif op=='+':value=dag.add(*(parameter_images[i]for i in data[0]))
        elif op=='*':value=dag.mul(*(parameter_images[i]for i in data[0]))
        elif op=='inv':value=dag.inv(parameter_images[data[0]])
        parameter_images[index]=value
    graph_inputs=[parameter_images[i]for i in parameter['parent35_input_nodes']]
    normalized_generators=[dag.sub(dag.inputs[parent['input_names'][i]],graph_inputs[i])for i in columns]
    equation_nodes=[parent['all350_residual_nodes'][i]for i in rows]
    needed=ancestors(parent['nodes'],equation_nodes)
    graph_values={};cofactors={};dependent_positions={j:k for k,j in enumerate(columns)}
    for index in sorted(needed):
        op,*data=parent['nodes'][index]
        if op=='q':graph=dag.q(data[0]);cofactor={}
        elif op=='v':
            graph=graph_inputs[data[0]]
            cofactor={dependent_positions[data[0]]:dag.one}if data[0]in dependent_positions else{}
            if data[0]not in dependent_positions:assert graph==index
        elif op=='+':
            graph=dag.add(*(graph_values[i]for i in data[0]))
            active=set().union(*(cofactors[i]for i in data[0]))
            cofactor={j:dag.add(*(cofactors[i][j]for i in data[0]if j in cofactors[i]))for j in active}
        elif op=='*':
            graph=dag.mul(*(graph_values[i]for i in data[0]));terms={}
            # Product difference: original prefix, changed factor, graph suffix.
            for position,child in enumerate(data[0]):
                prefix=data[0][:position];suffix=[graph_values[i]for i in data[0][position+1:]]
                for variable,value in cofactors[child].items():
                    terms.setdefault(variable,[]).append(dag.mul(value,*prefix,*suffix))
            cofactor={j:dag.add(*terms[j])for j in terms}
        elif op=='inv':
            graph=dag.inv(graph_values[data[0]])
            assert dag.nodes[graph][0]=='q' and Fraction(dag.nodes[graph][1])
            cofactor={j:dag.neg(dag.mul(value,index,graph))for j,value in cofactors[data[0]].items()}
        graph_values[index]=graph
        cofactors[index]={j:value for j,value in cofactor.items()if value!=dag.zero}
    matrix_nodes=[[cofactors[node].get(j,dag.zero)for j in range(26)]for node in equation_nodes]
    graph_equations=[graph_values[node]for node in equation_nodes]
    graph_needed=ancestors(dag.nodes,graph_equations)
    graph_polynomials=polynomial_values(dag.nodes,graph_needed,parameter['input_names'],
                                        {index:j for j,index in enumerate(parameter['parent_input_indices'])})
    assert all(not graph_polynomials[i]for i in graph_equations)
    assert json.loads(json.dumps(dag.nodes[:len(parent['nodes'])]))==parent['nodes']
    assert sum(node[0]=='inv'for node in dag.nodes)==sum(node[0]=='inv'for node in parent['nodes'])
    base={name:Fraction(value)for name,value in zip(parent['input_names'],jacobian['dense_parent35_point'])}
    assert base['a']==Fraction(-8,9)
    values=dag.evaluate(base)
    assert all(not values[i]for i in normalized_generators)
    matrix=s.Matrix([[values[i]for i in row]for row in matrix_nodes])
    jac=s.Matrix([[s.Rational(v)for v in row]for row in jacobian['dense_full350x26_dependent_Jacobian']])
    assert matrix==jac[rows,:]
    minor=DomainMatrix.from_Matrix(matrix).det();assert str(minor)==jacobian['dense_minor']
    # Dense exact QQ tests supplement the universal divided-difference induction.
    dense_tests=[]
    for a_value in [Fraction(1),Fraction(2),Fraction(3,7)]:
        dense={name:Fraction(((5*i+3)%13)-6,i%4+2)for i,name in enumerate(parent['input_names'])}
        dense['a']=a_value;current=dag.evaluate(dense)
        assert any(current[i]for i in equation_nodes)
        for r,node in enumerate(equation_nodes):
            assert current[node]-current[graph_equations[r]]==sum((current[c]*current[n]for c,n in zip(matrix_nodes[r],normalized_generators)),Fraction())
        dense_tests.append(str(a_value))
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
        for c,j in enumerate(columns):assert weights[matrix_nodes[r][c]]in(None,parent['residual_weight_bounds'][i]-parent['input_weights'][j])
    degree=sum(parent['residual_weight_bounds'][i]for i in rows)-sum(parent['input_weights'][j]for j in columns)
    assert degree==275
    epsilon_inputs=[]
    for name,weight in zip(parent['input_names'],parent['input_weights']):
        c=base[name]
        epsilon_inputs.append(tuple([Fraction(0)]*weight+[c])if c else())
    epsilon_values=epsilon_evaluate(dag.nodes,epsilon_inputs)
    assert all(not epsilon_values[i]for i in normalized_generators)
    for row in matrix_nodes:
        for index in row:
            poly=epsilon_values[index];weight=weights[index]
            assert not poly or all(not c or i==weight for i,c in enumerate(poly))
    output.mkdir(parents=True,exist_ok=False)
    payload=dict(id='Q_T40_DIVIDED_DIFFERENCE_OPEN_20261103_V63',input_names=parent['input_names'],
                 input_weights=parent['input_weights'],parent_archive_sha256=sha(parent_path),
                 nine_parameter_archive_sha256=sha(parameter_path),nodes=dag.nodes,
                 unchanged_parent_node_prefix=len(parent['nodes']),
                 dependent_parent_input_indices=columns,graph_parent35_input_nodes=graph_inputs,
                 normalized_graph_generator_nodes=normalized_generators,
                 selected_parent350_equation_indices=rows,selected_parent_equation_nodes=equation_nodes,
                 graph_selected_equation_nodes=graph_equations,cofactor26_matrix_nodes=matrix_nodes,
                 exact_identity='selected_parent_equations = cofactor26_matrix * normalized_graph_generators',
                 determinant_definition='C=det(cofactor26_matrix_nodes)',determinant_relative_weight=degree,
                 source_ancestor_count=len(needed))
    archive=output/'cofactor-model.json.gz';encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode()
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    note=ROOT/'research/notes/Q_T40_DIVIDED_DIFFERENCE_OPEN_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,parameter_path,jacobian_path,
            ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',
            ROOT/'code/q_normality_open_certificate_20261103_v63.py',
            ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/report.json']
    report=dict(id=payload['id'],all_passed=True,
                exact_divided_difference_induction_all_needed_nodes=True,
                parent_source_ancestor_count=len(needed),cofactor_matrix_shape=[26,26],
                parent_prefix_unchanged=True,added_nodes=len(dag.nodes)-len(parent['nodes']),
                new_inverses_added=False,all350_graph_identities_replayed_exact_QQ9=True,
                all26_graph_equation_remainders_independently_exact_QQ9_zero=True,
                normalized_generators_are_polynomial_graph_kernel=True,
                dense_QQ_supplementary_a_values=dense_tests,
                cofactor_matrix_at_dense_T40_point_matches_full_Jacobian_dependent_columns=True,
                exact_cofactor_determinant_at_dense_T40_point=str(minor),
                determinant_relative_weight=degree,exact_weighted_ray_power=degree,
                exact_scheme_identity='I[1/C]=(N)[1/C], hence R[1/C]=QQ[nine parameters][1/C|_graph]',
                normality_classification_used_for_this_scheme_identity=False,
                every_other_component_or_nonreduced_part_lies_in_C_zero=True,
                full35_branch_or_primary_exhaustion_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={'cofactor-model.json.gz':sha(archive)},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact divided-difference T40 scheme open, no normality assumption;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
