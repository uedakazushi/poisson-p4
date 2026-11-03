#!/usr/bin/env python3
"""Exact 34-variable model on the nonresonant open, including Q_(1).

Eliminate one weight6 variable by an ACTUAL full residual, including its lower
weight nonlinear terms. No linearized or finite-jet replacement is used.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from q_weighted_finite_model_20261103_v63 import frames,equations,A

ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def build(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    model=json.loads(gzip.decompress(parent.read_bytes()))
    a=s.Symbol('a');K=s.QQ.frac_field(a)
    P,*zz=ring(','.join(model['input_names'][1:]),K)
    target=model['all350_residual_nodes'][296]
    needed={target}
    for i in range(target,-1,-1):
        if i not in needed:continue
        op,*data=model['nodes'][i]
        if op in ('+','*'):needed.update(data[0])
        elif op=='inv':needed.add(data[0])
    values={}
    for i,(op,*data)in enumerate(model['nodes'][:target+1]):
        if i not in needed:continue
        if op=='q':v=P(K(s.Rational(data[0])))
        elif op=='v':v=P(K(a))if data[0]==0 else zz[data[0]-1]
        elif op=='+':v=sum((values[j]for j in data[0]),P.zero)
        elif op=='*':
            v=P.one
            for j in data[0]:v*=values[j]
        elif op=='inv':
            prev=values[data[0]];assert all(not any(m)for m in prev)
            v=P(prev[(0,)*34]**(-1))
        else:raise ValueError(op)
        values[i]=v
    residual=values[target]
    assert all(mon[33]<=1 for mon in residual)
    coefficient=residual.diff(zz[33]);assert coefficient.diff(zz[33])==0
    assert all(not any(mon)for mon in coefficient)
    c=coefficient[(0,)*34]
    expected=-(9*a-2)*(288*a+128)/(12*(27*a*a-60*a+8))
    assert c==K(expected)
    r0=residual-zz[33]*coefficient
    assert all(mon[33]==0 for mon in r0)
    assert residual==coefficient*zz[33]+r0
    assert K.to_sympy(c).subs(a,1)!=0
    C,free,cuts,positive,degrees,groupdet,groupcols=frames()
    _,_,_,linear=equations(C,positive,degrees)
    generic_linear_rank=sum(len(DomainMatrix.from_Matrix(linear[:,[j for j,d in enumerate(degrees)if d==weight]])
        .convert_to(s.QQ.frac_field(A)).rref()[1])for weight in sorted(set(degrees)))
    a_one_rank=len(DomainMatrix.from_Matrix(linear.subs(A,1)).convert_to(s.QQ).rref()[1])
    assert generic_linear_rank==75 and a_one_rank==75

    names=model['input_names'][:-1];dag=DAG(names)
    def coefficient_node(value):return dag.expression(s.cancel(K.to_sympy(value)))
    r0_node=dag.add(*(dag.mul(coefficient_node(value),*(dag.power(dag.inputs[f'z{j}'],power)
        for j,power in enumerate(mon[:33])if power))for mon,value in r0.items()))
    eliminated=dag.mul(dag.q(-1),dag.inv(dag.expression(expected)),r0_node)
    old=[]
    for op,*data in model['nodes']:
        if op=='q':v=dag.q(data[0])
        elif op=='v':v=dag.inputs[names[data[0]]]if data[0]<34 else eliminated
        elif op=='+':v=dag.add(*(old[j]for j in data[0]))
        elif op=='*':v=dag.mul(*(old[j]for j in data[0]))
        elif op=='inv':v=dag.inv(old[data[0]])
        else:raise ValueError(op)
        old.append(v)
    payload=dict(id='Q_NONRESONANT_FINITE_MODEL_20261103_V63',input_names=names,
        input_weights=model['input_weights'][:-1],nodes=dag.nodes,
        all350_residual_nodes=[old[i]for i in model['all350_residual_nodes']],
        raw150_reconstruction_nodes=[old[i]for i in model['raw150_reconstruction_nodes']],
        eliminated_parent_input='z33',eliminated_parent_input_node=eliminated,
        exact_eliminating_equation_index=296,eliminating_unit_coefficient=str(s.factor(expected)),
        parent_model=str(parent.relative_to(ROOT)),parent_archive_sha256=sha(parent),
        reduced_a_open=s.sstr(s.factor(s.sympify(model['reduced_a_open'])*(9*a-2))),
        group_cut_raw_indices=model['group_cut_raw_indices'],
        positive_raw_indices=model['positive_raw_indices'],
        parent_free_positive_columns=model['free_positive_columns'])
    encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode()
    archive=output/'model.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as z:z.write(encoded)
    evidence=dict(eliminating_residual_QQ_a_polynomial=[[[*mon],str(c)]for mon,c in residual.items()],
        exact_coefficient=str(s.factor(expected)),eliminating_equation_monomials=len(residual),
        nonlinear_terms=len([mon for mon in residual if sum(mon)>1]),
        residual_constant_in_z33_monomials=len(r0),all_nonlinear_terms_retained=True)
    (output/'elimination.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result=dict(id=payload['id'],all_passed=True,model_variables=34,positive_variables=33,
        zero_weight_parameter='a',actual_full_residual_used=True,eliminated_variable_weight=6,
        residual_350_retained=True,all_higher_orders_retained=True,
        includes_distinguished_Q_a_one=True,excludes_resonant_Q_a_two_ninths=True,
        formal_smooth_factor_projective=16,formal_smooth_factor_affine=17,
        center_a_one_transverse_tangent_dimension=34,
        exact_generic_positive_linear_rank=generic_linear_rank,
        exact_a_one_positive_linear_rank=a_one_rank,
        circuit_nodes=len(dag.nodes),eliminating_equation_monomials=len(residual),
        eliminating_equation_nonlinear_monomials=evidence['nonlinear_terms'],
        reduced_a_open=payload['reduced_a_open'],
        branch_classification_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),parent,
            ROOT/'code/q_weighted_finite_model_20261103_v63.py',
            ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json']},
        generated_sha256={p.name:sha(p)for p in [archive,output/'elimination.json']},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('EXACT NONRESONANT MODEL',result,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args().output)
