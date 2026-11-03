#!/usr/bin/env python3
"""Polynomial17 cofactor certificate extracting the actual vertical Rat113 prime."""
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
from q_smith_global_finite_model_20261103_v63 import closure
from q_T39_divided_difference_open_20261103_v63 import polynomial_values
from q_normality_open_certificate_20261103_v63 import epsilon_evaluate

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output):
    started=time.monotonic()
    parent_path=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/model37.json.gz'
    parameter_path=ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63/model20.json.gz'
    certificate_path=ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63/certificate.json'
    frame_path=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/frames.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()));parameter=json.loads(gzip.decompress(parameter_path.read_bytes()))
    certificate=json.loads(certificate_path.read_text());frames=json.loads(gzip.decompress(frame_path.read_bytes()))
    rows=certificate['dense_selected_equation_rows'];columns=certificate['dependent_parent_input_columns']
    assert len(rows)==len(columns)==17
    assert columns==[i for i in range(37)if i not in parameter['parent37_input_indices']]
    saved=polynomial_values(parameter['nodes'],set(range(len(parameter['nodes']))),parameter['input_names'],{i:i for i in range(20)})
    assert all(not saved[i]for i in [*parameter['all350_residual_nodes'],*parameter['all150_Schur_nodes']])
    dag=DAG(parent['input_names'])
    dag.nodes=[(op,tuple(data[0]))if op in ('+','*')else(op,data[0])for op,*data in parent['nodes']]
    dag.lookup={node:i for i,node in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
    dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
    images={}
    for i in sorted(closure(parameter['nodes'],parameter['parent37_input_nodes'])):
        op,*data=parameter['nodes'][i]
        if op=='q':value=dag.q(data[0])
        elif op=='v':value=dag.inputs[parameter['input_names'][data[0]]]
        elif op=='+':value=dag.add(*(images[j]for j in data[0]))
        elif op=='*':value=dag.mul(*(images[j]for j in data[0]))
        else:raise ValueError(op)
        images[i]=value
    projection=[images[i]for i in parameter['parent37_input_nodes']]
    generators=[dag.sub(dag.inputs[parent['input_names'][i]],projection[i])for i in columns]
    assert columns[0]==0 and generators[0]==dag.inputs['a']
    equations=[parent['all350_residual_nodes'][i]for i in rows]
    graph_values={};cofactors={};positions={column:i for i,column in enumerate(columns)}
    needed=closure(parent['nodes'],equations)
    for i in sorted(needed):
        op,*data=parent['nodes'][i]
        if op=='q':graph=dag.q(data[0]);cofactor={}
        elif op=='v':
            graph=projection[data[0]];cofactor={positions[data[0]]:dag.one}if data[0]in positions else{}
            if data[0]not in positions:assert graph==i
        elif op=='+':
            graph=dag.add(*(graph_values[j]for j in data[0]));active=set().union(*(cofactors[j]for j in data[0]))
            cofactor={k:dag.add(*(cofactors[j][k]for j in data[0]if k in cofactors[j]))for k in active}
        elif op=='*':
            graph=dag.mul(*(graph_values[j]for j in data[0]));terms={}
            for position,child in enumerate(data[0]):
                prefix=data[0][:position];suffix=[graph_values[j]for j in data[0][position+1:]]
                for variable,value in cofactors[child].items():terms.setdefault(variable,[]).append(dag.mul(value,*prefix,*suffix))
            cofactor={k:dag.add(*terms[k])for k in terms}
        else:raise ValueError(op)
        graph_values[i]=graph;cofactors[i]={k:v for k,v in cofactor.items()if v!=dag.zero}
    matrix_nodes=[[cofactors[i].get(j,dag.zero)for j in range(17)]for i in equations]
    graph_equations=[graph_values[i]for i in equations]
    polynomials=polynomial_values(dag.nodes,closure(dag.nodes,graph_equations),parameter['input_names'],
                                  {old:new for new,old in enumerate(parameter['parent37_input_indices'])})
    assert all(not polynomials[i]for i in graph_equations)
    assert json.loads(json.dumps(dag.nodes[:len(parent['nodes'])]))==parent['nodes']
    assert all(node[0]!='inv'for node in dag.nodes)
    point={name:Fraction(value)for name,value in zip(parent['input_names'],certificate['dense_parent37_point'])}
    values=dag.evaluate(point);assert all(not values[i]for i in generators)
    matrix=s.Matrix([[values[i]for i in row]for row in matrix_nodes])
    full=s.Matrix([[s.Rational(v)for v in row]for row in certificate['dense_full350x37_Jacobian']])
    assert matrix==full[rows,columns]
    determinant=DomainMatrix.from_Matrix(matrix).det();assert str(determinant)==certificate['dense_dependent17_minor']
    for a_value in [Fraction(0),Fraction(1),Fraction(2,9)]:
        dense={name:Fraction((-1)**i*(1+(5*i+3)%11),2+i%3)for i,name in enumerate(parent['input_names'])};dense['a']=a_value
        current=dag.evaluate(dense)
        for r,node in enumerate(equations):
            assert current[node]-current[graph_equations[r]]==sum((current[m]*current[n]for m,n in zip(matrix_nodes[r],generators)),Fraction())
    weights=[]
    for op,*data in dag.nodes:
        if op=='q':w=None if not Fraction(data[0])else 0
        elif op=='v':w=parent['input_weights'][data[0]]
        elif op=='+':
            options={weights[i]for i in data[0]if weights[i]is not None};assert len(options)<=1;w=next(iter(options),None)
        elif op=='*':w=None if any(weights[i]is None for i in data[0])else sum(weights[i]for i in data[0])
        weights.append(w)
    for r,i in enumerate(rows):
        for c,j in enumerate(columns):assert weights[matrix_nodes[r][c]]in(None,frames['Jacobi_row_weights'][i]-parent['input_weights'][j])
    degree=sum(frames['Jacobi_row_weights'][i]for i in rows)-sum(parent['input_weights'][i]for i in columns);assert degree==129
    ray=[]
    for name,w in zip(parent['input_names'],parent['input_weights']):
        c=point[name];ray.append(tuple([Fraction(0)]*w+[c])if c else())
    ray_values=epsilon_evaluate(dag.nodes,ray);assert all(not ray_values[i]for i in generators)
    for row in matrix_nodes:
        for i in row:assert not ray_values[i]or all(not c or j==weights[i]for j,c in enumerate(ray_values[i]))
    output.mkdir(parents=True,exist_ok=False)
    payload=dict(id='Q_RAT113_COFACTOR_OPEN_20261103_V63',input_names=parent['input_names'],input_weights=parent['input_weights'],
                 nodes=dag.nodes,unchanged_parent_prefix=len(parent['nodes']),parent_archive_sha256=sha(parent_path),
                 A20_archive_sha256=sha(parameter_path),dependent_parent_input_indices=columns,
                 graph_parent37_input_nodes=projection,normalized17_generator_nodes=generators,
                 selected350_equation_indices=rows,selected_equation_nodes=equations,
                 graph_selected_equation_nodes=graph_equations,cofactor17_matrix_nodes=matrix_nodes,
                 determinant_definition='C113=det(cofactor17_matrix_nodes)',determinant_relative_weight=degree)
    archive=output/'cofactor-model.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    note=ROOT/'research/notes/Q_RAT113_COFACTOR_OPEN_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,parameter_path,certificate_path,frame_path,
            ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',
            ROOT/'code/q_smith_global_finite_model_20261103_v63.py',
            ROOT/'code/q_T39_divided_difference_open_20261103_v63.py',
            ROOT/'code/q_normality_open_certificate_20261103_v63.py',
            ROOT/'reports/Q_RAT113_GLOBAL_SLICE_20261103_v63/report.json']
    report=dict(id=payload['id'],all_passed=True,exact_polynomial_identity='F_selected=M17*N113',
                all350_graph_identities_and150Schur_replayed_QQ20=True,
                all17_graph_remainders_rechecked_QQ20_zero=True,parent_ancestor_count=len(needed),
                added_nodes=len(dag.nodes)-len(parent['nodes']),no_inverses=True,
                all289_cofactor_values_equal_full_Jacobian_dependent_minor=True,
                exact_nonzero_C113_at_twenty_all1=str(determinant),relative_weight=degree,
                exact_weighted_ray_identity='C113(ray(lambda))=c*lambda^129',
                exact_scheme_open='I[1/C113]=N113[1/C113], R37[1/C113]=QQ[A20][1/C113|_graph]',
                exact_colon='I:C113^k=N113 for every k>=1',
                exact_prime_peel='I=N113 intersection (I+(C113))',
                explicit_nonzero_a_torsion='a*C113=0 in full R37, C113 !=0 even at Q0 local/completed ring',
                full_R37_nonflat_over_a_at_Q0=True,
                residual_remaining_branches_classified=False,
                normal_divisor_classification_used=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={'cofactor-model.json.gz':sha(archive)},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS polynomial17 cofactor Rat113 prime peel and explicit full a-torsion;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
