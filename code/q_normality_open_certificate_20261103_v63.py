#!/usr/bin/env python3
"""A finite normality-open certificate in the exact35-variable Q model.

The frozen parent is copied and only new output circuits are appended. The
minor is specified by a66x66 matrix of circuit outputs, without expanding its
determinant in35 variables. This does not classify the full local branches.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MON2=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MON2))
MON4=[tuple(m.count(i)for i in range(3))for m in itertools.combinations_with_replacement(range(3),4)]
MON6=[tuple(m.count(i)for i in range(3))for m in itertools.combinations_with_replacement(range(3),6)]
MON10=[tuple(m.count(i)for i in range(3))for m in itertools.combinations_with_replacement(range(3),10)]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def trim(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return tuple(p)


def add(p,q):
    out=[Fraction(0)]*max(len(p),len(q))
    for i,c in enumerate(p):out[i]+=c
    for i,c in enumerate(q):out[i]+=c
    return trim(out)


def multiply(p,q):
    if not p or not q:return()
    out=[Fraction(0)]*(len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,d in enumerate(q):out[i+j]+=c*d
    return trim(out)


def epsilon_evaluate(nodes,inputs):
    values=[]
    for op,*data in nodes:
        if op=='q':value=trim([Fraction(data[0])])
        elif op=='v':value=inputs[data[0]]
        elif op=='+':
            value=()
            for j in data[0]:value=add(value,values[j])
        elif op=='*':
            value=(Fraction(1),)
            for j in data[0]:value=multiply(value,values[j])
        elif op=='inv':
            previous=values[data[0]];assert len(previous)==1 and previous[0]
            value=(1/previous[0],)
        else:raise ValueError(op)
        values.append(value)
    return values


def node_tuple(node):
    op,*data=node
    return(op,tuple(data[0]))if op in ('+','*')else(op,data[0])


def build(output,model_path):
    started=time.monotonic();parent=json.loads(gzip.decompress(model_path.read_bytes()))
    resonance_path=ROOT/'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.json'
    resonance=json.loads(resonance_path.read_text());assert resonance['all_passed']
    index={value:i for i,value in enumerate(BASIS)}
    formula={(h,m):defaultdict(int)for h in range(5)for m in MON4}
    # F=sum_i(-1)^i*x_i*Pf(P with row/column i omitted). Differentiate
    # in ambient x_h, THEN restrict x3=x4=0.
    for omitted in range(5):
        j,k,l,m=[i for i in range(5)if i!=omitted]
        for left,right,sign in [((j,k),(l,m),1),((j,l),(k,m),-1),((j,m),(k,l),1)]:
            for mon in MON2:
                for other in MON2:
                    total=sorted([omitted,*mon,*other])
                    for h in set(total):
                        remaining=total.copy();remaining.remove(h)
                        if any(v>2 for v in remaining):continue
                        exponents=tuple(remaining.count(i)for i in range(3))
                        pair=tuple(sorted([index[left,mon],index[right,other]]))
                        formula[h,exponents][pair]+=(-1)**omitted*sign*total.count(h)
    formula={key:{mon:c for mon,c in terms.items()if c}for key,terms in formula.items()}

    dag=DAG(parent['input_names'])
    dag.nodes=[node_tuple(node)for node in parent['nodes']]
    dag.lookup={node:i for i,node in enumerate(dag.nodes)};assert len(dag.lookup)==len(dag.nodes)
    dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
    dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)}
    raw=parent['raw150_reconstruction_nodes'];partial_nodes=[]
    for h in range(5):
        partial_nodes.append([dag.add(*(dag.mul(dag.q(c),raw[i],raw[j])
                              for(i,j),c in sorted(formula[h,m].items())))for m in MON4])
    scale=162
    scaled_nodes=[[dag.mul(dag.q(scale),node)for node in row]for row in partial_nodes]
    assert json.loads(json.dumps(dag.nodes[:len(parent['nodes'])]))==parent['nodes']
    assert sum(node[0]=='inv'for node in dag.nodes)==sum(node[0]=='inv'for node in parent['nodes'])

    inputs=[(Fraction(2,9),)]+[()for _ in range(34)]
    for name,c in zip(['z11','z12','z23','z33'],[4,6,-4,1]):
        inputs[parent['input_names'].index(name)]=(Fraction(0),Fraction(c))
    values=epsilon_evaluate(dag.nodes,inputs)
    assert all(values[i]==()for i in parent['all350_residual_nodes'])
    X=s.symbols('x0:5');eps=s.Symbol('eps')
    actual_F=s.sympify(resonance['degeneracy_polynomial'])
    for h in range(5):
        physical=s.Poly(s.diff(actual_F,X[h]).subs({X[3]:0,X[4]:0}),*X[:3])
        for j,mon in enumerate(MON4):
            polynomial=s.Poly(physical.coeff_monomial(mon),eps)
            expected=trim([Fraction(polynomial.nth(k))for k in range(max(polynomial.degree(),0)+1)])
            assert values[partial_nodes[h][j]]==expected

    weights=[]
    for op,*data in dag.nodes:
        if op=='q':weight=None if not Fraction(data[0])else 0
        elif op=='v':weight=parent['input_weights'][data[0]]
        elif op=='+':
            options={weights[j]for j in data[0]if weights[j]is not None};assert len(options)<=1
            weight=next(iter(options),None)
        elif op=='*':weight=None if any(weights[j]is None for j in data[0])else sum(weights[j]for j in data[0])
        elif op=='inv':assert weights[data[0]]==0;weight=0
        weights.append(weight)
    for h in range(5):
        for j,mon in enumerate(MON4):assert weights[scaled_nodes[h][j]]in(None,h+mon[1]+2*mon[2])
    positions4={mon:i for i,mon in enumerate(MON4)}
    matrix_nodes=[];row_weights=[]
    for h in range(5):
        for shift in MON6:
            row=[]
            for mon in MON10:
                difference=tuple(x-y for x,y in zip(mon,shift))
                row.append(scaled_nodes[h][positions4[difference]]if min(difference)>=0 else dag.zero)
            matrix_nodes.append(row);row_weights.append(h-shift[1]-2*shift[2])
    assert len(matrix_nodes)==140 and all(len(row)==66 for row in matrix_nodes)
    matrix=s.Matrix([[sum(values[i])for i in row]for row in matrix_nodes])
    assert all(c.q==1 for c in matrix)
    pivot_rows=list(DomainMatrix.from_Matrix(matrix.T).convert_to(s.GF(32003)).rref()[1])
    assert len(pivot_rows)==66
    determinant=DomainMatrix.from_Matrix(matrix[pivot_rows,:]).convert_to(s.ZZ).det()
    assert determinant and int(determinant)%32003==9347
    degree=sum(row_weights[i]for i in pivot_rows)+sum(mon[1]+2*mon[2]for mon in MON10)
    assert degree==354 and degree%6==0
    for h in range(5):
        for j,mon in enumerate(MON4):
            weight=h+mon[1]+2*mon[2];polynomial=values[scaled_nodes[h][j]]
            assert not polynomial or(weight%6==0 and all(not c or k==weight//6 for k,c in enumerate(polynomial)))
    payload=dict(id='Q_NORMALITY_OPEN_CERTIFICATE_20261103_V63',
                 parent_model=str(model_path.relative_to(ROOT)),parent_archive_sha256=sha(model_path),
                 input_names=parent['input_names'],input_weights=parent['input_weights'],nodes=dag.nodes,
                 unchanged_parent_node_prefix=len(parent['nodes']),
                 raw150_reconstruction_nodes=raw,all350_parent_residual_nodes=parent['all350_residual_nodes'],
                 plane='x3=x4=0',quartic_monomial_exponents=MON4,
                 raw150_quadratic_formulas=[[[[i,j],c]for(i,j),c in sorted(formula[h,m].items())]
                                           for h in range(5)for m in MON4],
                 partial_coefficient_nodes=partial_nodes,scaled_partial_coefficient_nodes=scaled_nodes,
                 constant_scale=scale,matrix_rows=[[h,list(mon)]for h in range(5)for mon in MON6],
                 matrix_columns=MON10,matrix_nodes=matrix_nodes,selected_minor_rows=pivot_rows,
                 minor_relative_weight=degree,minor_definition='det(matrix_nodes[selected_minor_rows,:])')
    output.mkdir(parents=True,exist_ok=False)
    encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode();archive=output/'matrix-model.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    note=ROOT/'research/notes/Q_NORMALITY_OPEN_CERTIFICATE_20261103_v63.md'
    paths=[Path(__file__),note,model_path,resonance_path,
           ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',
           ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json',
           ROOT/'reports/Q_RESONANCE_Q51_INCIDENCE_20261103_v62.json',
           ROOT/'reports/NORMAL_QUINTIC_RECOGNITION_AUDIT_20261103_v62.json',
           ROOT/'reports/NORMAL_ACTUAL_EIGHT_SUPPORT_AUDIT_20261103_v62.json']
    report=dict(id=payload['id'],all_passed=True,
                fixed_plane='x3=x4=0',five_ambient_partials_restricted_after_differentiation=True,
                Macaulay_degree=10,matrix_shape=[140,66],constant_QQ_scale=scale,
                selected_minor_rows=pivot_rows,exact_integer_minor_at_Rstar_eps_one=str(determinant),
                nonzero_minor_mod32003=int(determinant)%32003,
                universal_minor_relative_weight=degree,eta_axis_minor_exact_power=degree//6,
                exact_eta_axis_identity='Delta(Qstar+eps*Rstar-75eps^2*T)=integer_minor*eps^59',
                all75_quartic_coefficients_verified_in_untruncated_QQ_eps=True,
                parent_node_prefix_unchanged=True,additional_circuit_nodes=len(dag.nodes)-len(parent['nodes']),
                new_circuit_inverses_added=False,
                Delta_nonzero_implies_actual_divisor_geometrically_normal=True,
                normal_open_reduced_support_in_T39_uses_audited_classification=True,
                any_reduced_component_not_contained_in_T39_avoids_this_minor_open=True,
                Delta_zero_locus_exhaustion_claimed=False,full35_branch_or_primary_classification_claimed=False,
                scalar_Qa_axis_has_Delta_zero=True,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
                evidence_sha256={'matrix-model.json.gz':sha(archive)},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS finite normality minor140x66, exact eps^59 witness;',report['seconds'],'seconds',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve(),args.model.resolve())
