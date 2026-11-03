#!/usr/bin/env python3
"""Exact Q51/Q52 pullbacks in the finite Q model; no branch exhaustion."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring

from mps_s61_arithmetic_circuits_20261103_v62 import DAG,transform,identity
from s61_model import S61,X,J,vertical_rref,bivector_from_coefficients,coefficients,pushforward

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
N=36
ZERO=(Fraction(0),)*N
ONE=(Fraction(1),)+(Fraction(0),)*(N-1)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def jconstant(value):return(Fraction(value),)+(Fraction(0),)*(N-1)
def jadd(*values):return tuple(sum((v[i]for v in values),Fraction())for i in range(N))
def jmul(left,right):
    return(left[0]*right[0],)+tuple(left[0]*right[i]+left[i]*right[0]for i in range(1,N))
def jinverse(value):
    assert value[0]
    return(1/value[0],)+tuple(-value[i]/value[0]**2 for i in range(1,N))


def evaluate(nodes,ground,point,first_order=False):
    values=[]
    for op,*data in nodes:
        if op=='q':v=jconstant(data[0])if first_order else ground.from_sympy(s.Rational(data[0]))
        elif op=='v':v=point[data[0]]
        elif op=='+':v=jadd(*(values[i]for i in data[0]))if first_order else sum((values[i]for i in data[0]),ground.zero)
        elif op=='*':
            v=ONE if first_order else ground.one
            for i in data[0]:v=jmul(v,values[i])if first_order else v*values[i]
        elif op=='inv':v=jinverse(values[data[0]])if first_order else 1/values[data[0]]
        else:raise ValueError(op)
        values.append(v)
    return values


def graph_polynomials(path):
    data=json.loads(path.read_text());assert data['all_passed']
    u=s.symbols('u1:62');beta=[u[i]for i in data['weight1_u_indices_zero_based']]
    polynomials=[s.Poly(s.sympify(text),*beta,domain=s.QQ)for text in data['quadratic_section_q']]
    return data,polynomials


def compose_graph(dag,parameters,data,polynomials):
    out=[]
    for target,poly in zip(data['weight2_u_indices_zero_based'],polynomials):
        terms=[]
        for mon,c in poly.terms():
            factors=[dag.q(c)]
            for index,power in zip(data['weight1_u_indices_zero_based'],mon):
                factors.extend([parameters[index]]*power)
            terms.append(dag.mul(*factors))
        out.append(dag.sub(parameters[target],dag.add(*terms)))
    return out


def finite_parametrization(source,parent,raw,residuals,relations,weights,jac,rows,columns):
    free=[i for i in range(35)if i not in columns]
    names=[parent['input_names'][i]for i in free]
    dag=DAG(names);assigned=[dag.zero]*35;assigned[0]=dag.q(Fraction(2,9))
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
                Q51_graph_relation_nodes=graph_nodes,blocks=blocks,
                exact_QQ9_all350_and45_zero=True,
                physical_free_coordinate_left_inverse=True,
                total_raw150_nonzero_polynomial_terms=sum(len(values[i])for i in raw_nodes))


def build(output):
    started=time.monotonic()
    parent_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    dag=DAG(parent['input_names'])
    dag.nodes=[(op,tuple(data[0]))if op in ('+','*')else(op,data[0])for op,*data in parent['nodes']]
    dag.lookup={node:i for i,node in enumerate(dag.nodes)}
    dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
    dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
    raw=parent['raw150_reconstruction_nodes']
    a,c,d=raw[90],raw[60],raw[75]
    assert a==dag.inputs['a'] and raw[105]==dag.q(Fraction(-2,9))
    assert raw[120]==raw[135]==dag.zero
    M=identity(dag,5);Mi=identity(dag,5)
    for i in range(1,5):
        for j in range(1,5):M[i][j]=Mi[i][j]=dag.zero
    M[1][1]=dag.one;M[2][4]=dag.inv(a);M[3][2]=dag.one
    M[3][4]=dag.neg(dag.mul(c,dag.inv(a)))
    M[4][3]=dag.q(Fraction(-9,2));M[4][4]=dag.mul(dag.q(Fraction(9,2)),d,dag.inv(a))
    Mi[1][1]=dag.one;Mi[2][2]=c;Mi[2][3]=dag.one
    Mi[3][2]=d;Mi[3][4]=dag.q(Fraction(-2,9));Mi[4][2]=a
    # Independently verify the displayed universal frame and its inverse.
    aa,cc,dd=s.symbols('a c d')
    C=s.Matrix([[0,cc,dd,aa],[-cc,0,-s.Rational(2,9),0],[-dd,s.Rational(2,9),0,0],[-aa,0,0,0]])
    frame=s.Matrix([[1,0,0,0],[0,0,0,1/aa],[0,1,0,-cc/aa],[0,0,-s.Rational(9,2),9*dd/(2*aa)]])
    inverse=s.Matrix([[1,0,0,0],[0,cc,1,0],[0,dd,0,-s.Rational(2,9)],[0,aa,0,0]])
    assert (frame*inverse-s.eye(4)).applyfunc(s.cancel)==s.zeros(4)
    assert (frame*C*frame.T-J[1:5,1:5]).applyfunc(s.cancel)==s.zeros(4)
    assert s.factor(frame.det())==-s.Rational(9,2)/aa

    based=transform(dag,raw,M,Mi)
    model=S61();rref,pivots=vertical_rref()
    def quotient(vector):
        return[dag.sub(value,dag.add(*(dag.mul(dag.q(rref[i,k]),vector[p])
                      for i,p in enumerate(pivots)if rref[i,k])))for k,value in enumerate(vector)]
    initial=quotient(based)
    selected=[initial[i]for i in model.certificate['perspective_orbit_pivots_raw']]
    perspective_inverse=model.perspective_matrix().inv()
    alpha=[dag.neg(dag.add(*(dag.mul(dag.q(perspective_inverse[i,j]),selected[j])for j in range(4))))for i in range(4)]
    PM=identity(dag,5);PI=identity(dag,5)
    for i in range(4):PM[0][i+1]=alpha[i];PI[0][i+1]=dag.neg(alpha[i])
    normal=quotient(transform(dag,based,PM,PI))
    parameters=[normal[i]for i in model.certificate['raw_free_indices']]
    paths={label:ROOT/f'reports/{label}_UNIVERSAL_SECTION_20261101.json'for label in ['Q51','Q52']}
    graphs={label:graph_polynomials(path)for label,path in paths.items()}
    relations={label:compose_graph(dag,parameters,*data)for label,data in graphs.items()}
    assert json.loads(json.dumps(dag.nodes[:len(parent['nodes'])]))==parent['nodes']
    dependencies=[];weights=[]
    for index,(op,*data)in enumerate(dag.nodes):
        if op=='q':dep=0;weight=None if not Fraction(data[0])else 0
        elif op=='v':dep=1<<data[0];weight=parent['input_weights'][data[0]]
        elif op in ('+','*'):
            assert all(j<index for j in data[0]);dep=0
            for j in data[0]:dep|=dependencies[j]
            if op=='+':
                options={weights[i]for i in data[0]if weights[i]is not None};assert len(options)<=1
                weight=next(iter(options),None)
            else:weight=None if any(weights[i]is None for i in data[0])else sum(weights[i]for i in data[0])
        elif op=='inv':
            assert data[0]<index;dep=dependencies[data[0]];assert dep&~1==0
            assert weights[data[0]]==0;weight=0
        dependencies.append(dep);weights.append(weight)
    new_inverses=[node for node in dag.nodes[len(parent['nodes']):]if node[0]=='inv']
    assert all(node==('inv',a)for node in new_inverses)
    assert all(dependencies[indices[25]]==1 and weights[indices[25]]==0 for indices in relations.values())

    point=[]
    for i in range(35):
        v=list(ZERO);v[i+1]=Fraction(1)
        if i==0:v[0]=Fraction(2,9)
        point.append(tuple(v))
    jets=evaluate(dag.nodes,None,point,True)
    assert all(jets[i]==ZERO for i in parent['all350_residual_nodes'])
    normalized_jets=[jets[i]for i in normal];parameter_jets=[jets[i]for i in parameters]
    assert all(normalized_jets[i]==ZERO for i in model.certificate['perspective_orbit_pivots_raw'])
    def encoded_jet(poly):
        terms=[]
        for mon,nu,de in poly:
            v=jconstant(Fraction(nu,de))
            for i in mon:v=jmul(v,parameter_jets[i])
            terms.append(v)
        return jadd(*terms)
    reconstructed=[jconstant(v)for v in coefficients(X[0]**2*J)]
    for i,v in zip(model.certificate['raw_free_indices'],parameter_jets):reconstructed[i]=v
    for item in model.data['eliminated_raw_variables']:reconstructed[item['raw_index']]=encoded_jet(item['polynomial'])
    assert reconstructed==normalized_jets
    assert all(encoded_jet(poly)==ZERO for poly in model.reduced['equations'])
    assert all(not jets[i][0]for i in relations['Q51'])
    jac=s.Matrix([jets[i][1:]for i in relations['Q51']])
    rows=list(DomainMatrix.from_Matrix(jac.T).rref()[1]);assert len(rows)==26
    columns=list(DomainMatrix.from_Matrix(jac[rows,:]).rref()[1]);assert len(columns)==26
    minor=DomainMatrix.from_Matrix(jac[rows,columns]).det()
    assert str(minor)=='-50031545098999707/95420416'

    # Exact zero-section arithmetic in QQ(a), independently normalized physically.
    a_symbol=s.Symbol('a');K=s.QQ.frac_field(a_symbol)
    field_values=evaluate(dag.nodes,K,[K.gens[0]]+[K.zero]*34)
    scalar_raw=[K.to_sympy(field_values[i])for i in raw]
    assert all(not field_values[i]for i in parent['all350_residual_nodes'])
    physical=bivector_from_coefficients(scalar_raw)
    scalar_frame=s.zeros(5);scalar_frame[0,0]=scalar_frame[1,1]=1
    scalar_frame[2,4]=1/a_symbol;scalar_frame[3,2]=1;scalar_frame[4,3]=-s.Rational(9,2)
    physical_normal=model.normalize(pushforward(physical,scalar_frame))
    physical_checks=model.verify_point(physical_normal)
    assert all(s.cancel(K.to_sympy(field_values[node])-value)==0
               for node,value in zip(normal,physical_normal['coefficients']))
    assert physical_normal['perspective']==[0]*4
    scalar_relations={label:[s.factor(K.to_sympy(field_values[i]))for i in indices]
                      for label,indices in relations.items()}
    assert scalar_relations['Q51']==[s.factor((9*a_symbol-2)/(2*a_symbol**2))if i==25 else 0 for i in range(45)]
    assert scalar_relations['Q52']==[s.factor((9*a_symbol+8)/(2*a_symbol**2))if i==25 else 0 for i in range(45)]
    localized=s.sympify(parent['reduced_a_open'])
    assert all(localized.subs(a_symbol,value)!=0 for value in [s.Rational(2,9),-s.Rational(8,9)])
    parametrization=finite_parametrization(dag,parent,raw,parent['all350_residual_nodes'],
                                         relations['Q51'],weights,jac,rows,columns)

    payload=dict(id='Q_LOCAL_Q51_QUOTIENT_20261103_V63',input_names=parent['input_names'],
                 parent_archive_sha256=sha(parent_path),nodes=dag.nodes,
                 unchanged_parent_node_prefix=len(parent['nodes']),
                 parent_all350_residual_nodes=parent['all350_residual_nodes'],
                 raw150_parent_reconstruction_nodes=raw,frame_nodes=M,inverse_frame_nodes=Mi,
                 perspective_parameter_nodes=alpha,normal150_nodes=normal,normal61_nodes=parameters,
                 Q51_graph_relation_nodes=relations['Q51'],Q52_graph_relation_nodes=relations['Q52'])
    output.mkdir(parents=True,exist_ok=False)
    archive=output/'graph-model.json.gz';encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode()
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    parameter_archive=output/'T39-nine-parameter-model.json.gz'
    encoded=(json.dumps(parametrization,separators=(',',':'))+'\n').encode()
    with parameter_archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    certificate=dict(Q51_jacobian_45x35=[[str(v)for v in row]for row in jac.tolist()],
                     selected_rows=rows,selected_columns=columns,selected_minor=str(minor),
                     scalar_section_raw150=list(map(str,scalar_raw)),
                     scalar_section_normal61=list(map(str,physical_normal['parameters'])),
                     scalar_section_relations={label:list(map(str,values))for label,values in scalar_relations.items()},
                     frame=[[str(v)for v in row]for row in frame.tolist()],
                     inverse_frame=[[str(v)for v in row]for row in inverse.tolist()])
    certificate['T39_triangular_blocks']=parametrization['blocks']
    (output/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    note=ROOT/'research/notes/Q_LOCAL_Q51_QUOTIENT_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,ROOT/'code/s61_model.py',
            ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',*paths.values(),
            ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json',
            ROOT/'reports/Q_RESONANCE_Q51_INCIDENCE_20261103_v62.json',
            ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',
            ROOT/'reports/Q51_UNIVERSAL_SECTION_REVIEW_20261101.json',
            ROOT/'reports/Q52_UNIVERSAL_SECTION_REVIEW_20261101.json']
    report=dict(id=payload['id'],all_passed=True,full_Q51_relations=45,full_Q52_relations=45,
                parent_prefix_unchanged=True,added_node_count=len(dag.nodes)-len(parent['nodes']),
                only_new_inverse_is_a=True,new_inverse_count=len(new_inverses),fixed_e0_Pfaffian='-2*a/9',
                frame_determinant='-9/(2*a)',graph_Q51_jacobian_rank_at_Qstar=26,
                exact26_minor=str(minor),all150_first_order_reconstruction=True,
                all134_first_order_normalized_residuals_zero=True,
                all350_parent_first_order_residuals_zero=True,
                scalar_section_exact_QQ_a_physical_checks=physical_checks,
                scalar_section_Q51_only_nonzero_relation='(9*a-2)/(2*a**2)',
                scalar_section_Q52_only_nonzero_relation='(9*a+8)/(2*a**2)',
                scalar_section_nonzero_relation_ordinal_zero_based=25,
                same_graph_relation25_depends_only_on_a_in_full35_model=True,
                full_Q51_relation25_identity='(9*a-2)/(2*a**2)',
                full_Q52_relation25_identity='(9*a+8)/(2*a**2)',
                every_node_weight_homogeneous=True,
                exact_finite_global_T39_quotient_is_QQ_polynomial_ring_in9=True,
                T39_nine_free_parameters=parametrization['input_names'],
                T39_nine_parameter_weights=parametrization['input_weights'],
                T39_triangular_block_sizes=[len(b['rows'])for b in parametrization['blocks']],
                T39_all350_and45_exact_QQ9_zero=True,
                T39_explicit_physical_coordinate_left_inverse=True,
                T39_raw150_total_nonzero_polynomial_terms=parametrization['total_raw150_nonzero_polynomial_terms'],
                T39_local_quotient_regular_dimension=9,
                T39_local_smoothness_argument='Dimension >=26-17=9 from the known invariant component and the etale group slice; tangent dimension35-26=9; equality implies a regular local ring.',
                full35_other_components_excluded=False,full35_primary_decomposition_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in [archive,parameter_archive,output/'certificate.json']},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact Q51 local quotient rank26; symbolic Qa graph tests Q51/Q52;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
