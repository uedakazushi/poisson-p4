#!/usr/bin/env python3
"""Global34 exact quotient over QQ[a,1/(a*(9a-2))] of frozen Smith35."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from q_smith_global_finite_model_20261103_v63 import closure,physical,PAIRS,EXPS2,EXPS3,TRIPLES

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate_field(nodes,point,K):
    values=[]
    for op,*data in nodes:
        if op=='q':value=K.from_sympy(s.Rational(data[0]))
        elif op=='v':value=point[data[0]]
        elif op=='+':value=sum((values[i]for i in data[0]),K.zero)
        elif op=='*':
            value=K.one
            for i in data[0]:value*=values[i]
        elif op=='inv':value=K.one/values[data[0]]
        values.append(value)
    return values


def physical_field(raw,K):
    P={pair:{mon:raw[15*j+k]for k,mon in enumerate(EXPS2)if raw[15*j+k]}for j,pair in enumerate(PAIRS)}
    def entry(i,j):return P.get((i,j),{})if i<j else{m:-v for m,v in P.get((j,i),{}).items()}if i>j else{}
    derivatives={}
    for i in range(5):
        for j in range(5):
            for k in range(5):
                value={}
                for mon,c in entry(i,j).items():
                    if mon[k]:changed=list(mon);changed[k]-=1;value[tuple(changed)]=c*K.convert(mon[k])
                derivatives[i,j,k]=value
    for j in range(5):
        divergence=defaultdict(lambda:K.zero)
        for i in range(5):
            for mon,c in derivatives[i,j,i].items():divergence[mon]+=c
        assert all(not c for c in divergence.values())
    result=[]
    for i,j,k in TRIPLES:
        value=defaultdict(lambda:K.zero)
        for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
            for ell in range(5):
                for mon,co in entry(b,ell).items():
                    for other,second in derivatives[c,d,ell].items():value[tuple(x+y for x,y in zip(mon,other))]+=co*second
        result.extend(value.get(mon,K.zero)for mon in EXPS3)
    return result


def build(output):
    started=time.monotonic();folder=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63'
    parent_path=folder/'model35.json.gz';frame_path=folder/'frames.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()));frames=json.loads(gzip.decompress(frame_path.read_bytes()))
    records=frames['source_Smith_pivots'];index=next(i for i,p in enumerate(records)if p['kind']=='u')
    pivot=records[index];assert pivot['weight']==6
    removed=parent['free_weighted_columns'].index(pivot['column'])+1
    G=parent['Smith_transformed75_residual_nodes'][index]
    block=next(b for b in frames['Jacobi_complement_Smith_blocks']if b['weight']==6)
    smith_row=block['columns'].index(pivot['column'])
    combination=[(row,s.sympify(coefficient))for row,coefficient in zip(block['rows'],block['left'][smith_row])if s.sympify(coefficient)]
    a=s.Symbol('a');K=s.QQ.frac_field(a);diagonal=s.sympify(pivot['diagonal']);scale=s.cancel(diagonal/(9*a-2))
    assert scale.is_Rational and scale
    R,*variables=ring(','.join(parent['input_names'][1:]),K);values={}
    ancestors=closure(parent['nodes'],[G,*[parent['all350_residual_nodes'][row]for row,c in combination]])
    for node in sorted(ancestors):
        op,*data=parent['nodes'][node]
        if op=='q':value=R(K.from_sympy(s.Rational(data[0])))
        elif op=='v':value=R(K.gens[0])if data[0]==0 else variables[data[0]-1]
        elif op=='+':value=sum((values[i]for i in data[0]),R.zero)
        elif op=='*':
            value=R.one
            for i in data[0]:value*=values[i]
        elif op=='inv':assert values[data[0]].is_ground;value=R.one/values[data[0]]
        values[node]=value
    polynomial=values[G];unit_mon=tuple(int(i==removed-1)for i in range(34))
    assert polynomial==sum((R(K.from_sympy(c))*values[parent['all350_residual_nodes'][row]]for row,c in combination),R.zero)
    assert polynomial[unit_mon]==K.from_sympy(diagonal)
    assert all(not mon[removed-1]or mon==unit_mon for mon in polynomial)
    assert all(not s.denom(s.cancel(K.to_sympy(c))).has(a)for c in polynomial.values())
    rhs=polynomial.copy();del rhs[unit_mon]
    retained=[i for i in range(35)if i!=removed]
    dag=DAG(['a']+[f'z{i}'for i in range(33)])
    mapping={old:dag.inputs[dag.names[new]]for new,old in enumerate(retained)};mapping[removed]=dag.zero
    def clone(required):
        translated={}
        for i in sorted(closure(parent['nodes'],required)):
            op,*data=parent['nodes'][i]
            if op=='q':value=dag.q(data[0])
            elif op=='v':value=mapping[data[0]]
            elif op=='+':value=dag.add(*(translated[j]for j in data[0]))
            elif op=='*':value=dag.mul(*(translated[j]for j in data[0]))
            elif op=='inv':value=dag.inv(translated[data[0]])
            translated[i]=value
        return translated
    r0=clone([G])[G];u=dag.add(dag.mul(dag.q(9),dag.inputs['a']),dag.q(-2))
    solution=dag.neg(dag.mul(dag.q(1/scale),dag.inv(u),r0));mapping[removed]=solution
    translated=clone(list(range(len(parent['nodes']))))
    model=dict(id='Q_SMITH_GLOBAL_NONRESONANT_MODEL_20261103_V63',
               input_names=dag.names,input_weights=[parent['input_weights'][i]for i in retained],
               localized_base='QQ[a,1/(a*(9*a-2))]',nodes=dag.nodes,
               parent_archive_sha256=sha(parent_path),parent_removed_input_index=removed,
               parent_retained_input_indices=retained,parent_input_reconstruction_nodes=[mapping[i]for i in range(35)],
               complete_parent_node_composition_map=[translated[i]for i in range(len(parent['nodes']))],
               raw150_reconstruction_nodes=[translated[i]for i in parent['raw150_reconstruction_nodes']],
               all350_residual_nodes=[translated[i]for i in parent['all350_residual_nodes']],
               positive108_reconstruction_nodes=[translated[i]for i in parent['positive108_reconstruction_nodes']],
               Smith_transformed75_residual_nodes=[translated[i]for i in parent['Smith_transformed75_residual_nodes']],
               eliminated_relation_node=translated[G],eliminated_input_expression_node=solution)
    assert len(model['input_names'])==34
    weights=[];dependencies=[];inverses=[]
    for i,(op,*data)in enumerate(dag.nodes):
        if op=='q':weight=None if not Fraction(data[0])else 0;dep=0
        elif op=='v':weight=model['input_weights'][data[0]];dep=1<<data[0]
        elif op in ('+','*'):
            dep=0
            for j in data[0]:dep|=dependencies[j]
            if op=='+':
                options={weights[j]for j in data[0]if weights[j]is not None};assert len(options)<=1
                weight=next(iter(options),None)
            else:weight=None if any(weights[j]is None for j in data[0])else sum(weights[j]for j in data[0])
        elif op=='inv':assert dependencies[data[0]]==1 and weights[data[0]]==0;dep=1;weight=0;inverses.append(i)
        weights.append(weight);dependencies.append(dep)
    assert len(inverses)==2 and weights[solution]==6
    scalar_values=evaluate_field(dag.nodes,[K.gens[0]]+[K.zero]*33,K)
    inverse_arguments={str(s.factor(K.to_sympy(scalar_values[dag.nodes[i][1]])))for i in inverses}
    assert inverse_arguments=={'a','9*a - 2'}
    tests=[];parent_dag=DAG(parent['input_names']);parent_dag.nodes=parent['nodes']
    for value in [Fraction(1),Fraction(-4,9),Fraction(1,3),Fraction(44,81),Fraction(-8,9)]:
        point={name:Fraction((-1)**i*(1+(7*i+3)%11),2+i%3)for i,name in enumerate(dag.names)};point['a']=value
        current=dag.evaluate(point);old_point={name:current[node]for name,node in zip(parent['input_names'],model['parent_input_reconstruction_nodes'])}
        old_values=parent_dag.evaluate(old_point)
        assert all(current[translated[i]]==old_values[i]for i in range(len(parent['nodes'])))
        assert current[model['eliminated_relation_node']]==0
        raw=[current[i]for i in model['raw150_reconstruction_nodes']]
        assert physical(raw)==[current[i]for i in model['all350_residual_nodes']]
        tests.append(dict(a=str(value),all_parent_nodes_and150_and350_match=True,actual_physical_Jacobi_and_divergence_match=True,
                          eliminated_full_relation_zero=True,nonzero_residuals=sum(bool(current[i])for i in model['all350_residual_nodes'])))
    number_field=s.QQ.algebraic_field(s.sqrt(19));h=number_field.unit
    a_value=(number_field.convert(10)+number_field.convert(2)*h)/number_field.convert(9)
    assert number_field.convert(27)*a_value*a_value-number_field.convert(60)*a_value+number_field.convert(8)==number_field.zero
    point=[a_value]+[number_field.from_sympy(s.Rational((-1)**i*(1+(7*i+3)%11),2+i%3))for i in range(1,34)]
    current=evaluate_field(dag.nodes,point,number_field)
    old_point=[current[i]for i in model['parent_input_reconstruction_nodes']]
    old_values=evaluate_field(parent['nodes'],old_point,number_field)
    assert all(current[translated[i]]==old_values[i]for i in range(len(parent['nodes'])))
    assert current[model['eliminated_relation_node']]==number_field.zero
    raw=[current[i]for i in model['raw150_reconstruction_nodes']]
    assert physical_field(raw,number_field)==[current[i]for i in model['all350_residual_nodes']]
    output.mkdir(parents=True,exist_ok=False)
    def save(path,data):
        with path.open('wb')as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(data,separators=(',',':'))+'\n').encode())
    save(output/'model34.json.gz',model)
    elimination=dict(parent_input=parent['input_names'][removed],parent_input_index=removed,weight=6,
                     complete_relation_parent_node=G,complete_diagonal=str(diagonal),diagonal_over_u=str(scale),
                     relation_parent350_linear_combination=[[row,str(c)]for row,c in combination],
                     complete_polynomial_terms=[[list(mon),str(K.to_sympy(c))]for mon,c in sorted(polynomial.items())],
                     rhs_terms=[[list(mon),str(K.to_sympy(c))]for mon,c in sorted(rhs.items())],
                     coefficient_ring='QQ[a]',solution='-rhs/(a-2/9)',finite_jet_truncation_used=False)
    save(output/'elimination.json.gz',elimination)
    note=ROOT/'research/notes/Q_SMITH_GLOBAL_NONRESONANT_MODEL_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,frame_path,folder/'report.json',
            ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.json',
            ROOT/'code/q_smith_global_finite_model_20261103_v63.py',ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py']
    report=dict(id=model['id'],all_passed=True,variables=34,positive_variables=33,
                complete_high_order_relation_terms=len(polynomial),rhs_terms=len(rhs),
                exact_relation='(a-2/9)*z30 + r0(a,other33)=0',
                relation_weight=6,relation_max_ordinary_degree=max(map(sum,polynomial)),
                all_coefficients_polynomial_in_a=True,unknown_has_only_one_pure_linear_term=True,
                only_inverse_arguments=sorted(inverse_arguments),all_nodes_weight_homogeneous=True,
                source_parent_nodes=len(parent['nodes']),model_nodes=len(dag.nodes),
                all150_and350_parent_outputs_retained=True,complete_parent_node_composition_map_saved=True,
                exact_quotient_isomorphism='P/I ≅ B/substituted_I since eliminated full Smith relation is a unit times z30-solution',
                dense_QQ_tests=tests,
                quadratic_old_artifact_field='QQ(sqrt(19)), a=(10+2sqrt(19))/9; algebraic identity covers both embeddings',
                both_roots_27a2_minus60a_plus8_exact_field_validation=True,
                number_field_all_parent_nodes_and150_and350_and_physical_Jacobi_match=True,
                frozen_parent37_and35_modified=False,full_branch_or_flatness_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in output.glob('*.gz')},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS global34 only a*(9a-2) inverses, exact full relation and QQ/sqrt19 validation;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
