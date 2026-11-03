#!/usr/bin/env python3
"""Independent exact review of the nonresonant34-variable finite model.

Imports neither model producer nor its DAG implementation. Expands the actual
eliminating row by sparse dictionaries, compares all500 circuit outputs in an
independent arithmetic normal form, audits the localization and linear jets,
and reconstructs physical divergence/Jacobi over dense rational inputs.
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
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
R,a=ring('a',s.QQ);K=R.to_field()


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def needed_nodes(nodes,target):
    needed={target};stack=[target]
    while stack:
        op,*data=nodes[stack.pop()]
        refs=data[0]if op in ('+','*')else[data[0]]if op=='inv'else[]
        for j in refs:
            if j not in needed:needed.add(j);stack.append(j)
    return needed


def sparse_expansion(nodes,target,positive_count):
    """Ordinary full polynomials, with no degree bound or truncation."""
    zero=(0,)*positive_count
    def constant(c):return {zero:c}if c else{}
    def add(values):
        out={}
        for value in values:
            for mon,c in value.items():
                coefficient=out.get(mon,K.zero)+c
                if coefficient:out[mon]=coefficient
                else:out.pop(mon,None)
        return out
    def multiply(values):
        out=constant(K.one)
        for value in values:
            if not value:return{}
            product={}
            for mon,c in out.items():
                for other,d in value.items():
                    key=tuple(x+y for x,y in zip(mon,other))
                    coefficient=product.get(key,K.zero)+c*d
                    if coefficient:product[key]=coefficient
                    else:product.pop(key,None)
            out=product
        return out
    needed=needed_nodes(nodes,target);values={}
    for index,(op,*data)in enumerate(nodes):
        if index not in needed:continue
        if op=='q':value=constant(K(s.Rational(data[0])))
        elif op=='v':
            if data[0]==0:value=constant(K(a))
            else:
                mon=list(zero);mon[data[0]-1]=1;value={tuple(mon):K.one}
        elif op=='+':value=add(values[j]for j in data[0])
        elif op=='*':value=multiply(values[j]for j in data[0])
        elif op=='inv':
            previous=values[data[0]];assert len(previous)==1 and zero in previous
            value=constant(1/previous[zero])
        else:raise ValueError(op)
        values[index]=value
    return values[target],len(needed)


class ArithmeticNormalForm:
    """Exact linear combinations of products of opaque additive atoms.

    Rational functions of a are evaluated in QQ(a). Products are associative
    and commutative. Sums combine identical products. Multiplication of two
    genuine sums is retained as a product of additive atoms, avoiding a full
    expansion in33 positive variables. All rewrites are ring identities;
    dictionary equality is the test, not a numerical or hash identity.
    """
    def __init__(self,count):
        self.atoms=[('v',i)for i in range(count)]
        self.lookup={atom:i for i,atom in enumerate(self.atoms)}
        self.definitions={};self.normal_cache={}

    @staticmethod
    def scalar(c):return {():c}if c else{}

    @staticmethod
    def variable(i):return {(i,):K.one}

    @staticmethod
    def add(values):
        out={}
        for value in values:
            for mon,c in value.items():
                coefficient=out.get(mon,K.zero)+c
                if coefficient:out[mon]=coefficient
                else:out.pop(mon,None)
        return out

    def multiply(self,values):
        coefficient=K.one;factors=[]
        for value in values:
            if not value:return{}
            if len(value)==1:
                mon,c=next(iter(value.items()));coefficient*=c;factors.extend(mon)
            else:
                identity=id(value)
                if identity in self.normal_cache:leading,atom=self.normal_cache[identity]
                else:
                    leading=value[min(value)]
                    normalized={mon:c/leading for mon,c in value.items()}
                    key=('sum',tuple((mon,str(c))for mon,c in sorted(normalized.items())))
                    if key not in self.lookup:
                        atom=len(self.atoms);self.lookup[key]=atom;self.atoms.append(key)
                        self.definitions[atom]=normalized
                    atom=self.lookup[key];self.normal_cache[identity]=leading,atom
                coefficient*=leading;factors.append(atom)
        if not coefficient:return{}
        if len(factors)==1 and factors[0]in self.definitions:
            return {mon:coefficient*c for mon,c in self.definitions[factors[0]].items()}
        return {tuple(sorted(factors)):coefficient}

    def evaluate(self,nodes,inputs):
        # Retain all dictionaries while using the identity-keyed normalization
        # cache, so object ids cannot be recycled during either evaluation.
        values=[]
        for op,*data in nodes:
            if op=='q':value=self.scalar(K(s.Rational(data[0])))
            elif op=='v':value=inputs[data[0]]
            elif op=='+':value=self.add(values[j]for j in data[0])
            elif op=='*':value=self.multiply(values[j]for j in data[0])
            elif op=='inv':
                previous=values[data[0]];assert len(previous)==1 and()in previous
                value=self.scalar(1/previous[()])
            else:raise ValueError(op)
            values.append(value)
        return values


def denominator_supported(poly,localization):
    assert poly
    while not poly.is_ground:
        divisor=poly.gcd(localization);assert not divisor.is_ground,str(poly)
        poly=poly.exquo(divisor)


def audit_graph(model,localization):
    ground=[];weights=[];dependencies=[];inverse_count=0;inverse_pairs=set()
    for index,(op,*data)in enumerate(model['nodes']):
        if op=='q':
            value=K(s.Rational(data[0]));dependency=0;weight=0 if value else None
        elif op=='v':
            dependency=1<<data[0];weight=model['input_weights'][data[0]]
            value=K(a)if data[0]==0 else None
        elif op in ('+','*'):
            assert all(0<=j<index for j in data[0])
            dependency=0
            for j in data[0]:dependency|=dependencies[j]
            if op=='+':
                options={weights[j]for j in data[0]if weights[j]is not None}
                assert len(options)<=1;weight=next(iter(options),None)
            else:weight=None if any(weights[j]is None for j in data[0])else sum(weights[j]for j in data[0])
            if dependency&~1:value=None
            elif op=='+':value=sum((ground[j]for j in data[0]),K.zero)
            else:
                value=K.one
                for j in data[0]:value*=ground[j]
        elif op=='inv':
            j=data[0];assert 0<=j<index and dependencies[j]&~1==0 and weights[j]==0
            previous=ground[j];assert previous
            denominator_supported(previous.numer,localization)
            denominator_supported(previous.denom,localization)
            value=1/previous;dependency=dependencies[j];weight=0;inverse_count+=1
            inverse_pairs.add((str(previous.numer),str(previous.denom)))
        else:raise ValueError(op)
        ground.append(value);weights.append(weight);dependencies.append(dependency)
    return dict(weights=weights,inverse_count=inverse_count,distinct_inverse_functions=len(inverse_pairs))


def linear_jets(nodes,positive_count):
    """Constants and first positive-variable derivatives over QQ(a)."""
    values=[]
    for op,*data in nodes:
        if op=='q':value=K(s.Rational(data[0])),{}
        elif op=='v':value=(K(a),{})if data[0]==0 else(K.zero,{data[0]-1:K.one})
        elif op=='+':
            constant=K.zero;linear={}
            for j in data[0]:
                c,v=values[j];constant+=c
                for k,x in v.items():linear[k]=linear.get(k,K.zero)+x
            value=constant,{k:x for k,x in linear.items()if x}
        elif op=='*':
            constant=K.one;linear={}
            for j in data[0]:
                c,v=values[j];updated={k:x*c for k,x in linear.items()if x*c}
                for k,x in v.items():updated[k]=updated.get(k,K.zero)+constant*x
                constant*=c;linear={k:x for k,x in updated.items()if x}
            value=constant,linear
        elif op=='inv':
            c,v=values[data[0]];assert c and not v;value=1/c,{}
        else:raise ValueError(op)
        values.append(value)
    return values


def evaluate_fraction(nodes,inputs):
    values=[]
    for op,*data in nodes:
        if op=='q':value=Fraction(data[0])
        elif op=='v':value=inputs[data[0]]
        elif op=='+':value=sum((values[j]for j in data[0]),Fraction(0))
        elif op=='*':
            value=Fraction(1)
            for j in data[0]:value*=values[j]
        elif op=='inv':assert values[data[0]];value=1/values[data[0]]
        else:raise ValueError(op)
        values.append(value)
    return values


PAIRS=list(itertools.combinations(range(5),2))
MON2=list(itertools.combinations_with_replacement(range(5),2))
TRIPLES=list(itertools.combinations(range(5),3))
MON3=list(itertools.combinations_with_replacement(range(5),3))


def physical_equations(raw):
    """Direct coefficient multiplication/differentiation of the raw tensor."""
    tensor={(i,j):[]for i in range(5)for j in range(5)}
    divergence=defaultdict(Fraction)
    for ((i,j),mon),coefficient in zip(itertools.product(PAIRS,MON2),raw):
        if not coefficient:continue
        tensor[i,j].append((mon,coefficient));tensor[j,i].append((mon,-coefficient))
        for direction,target,sign in [(i,j,1),(j,i,-1)]:
            if direction in mon:
                remaining=list(mon);remaining.remove(direction)
                divergence[target,remaining[0]]+=sign*mon.count(direction)*coefficient
    assert all(not value for value in divergence.values())
    jacobi=[]
    for i,j,k in TRIPLES:
        coefficients=defaultdict(Fraction)
        for left,right,last in [(i,j,k),(j,k,i),(k,i,j)]:
            for direction in range(5):
                for mon,c in tensor[left,direction]:
                    for other,d in tensor[right,last]:
                        if direction not in other:continue
                        remaining=list(other);remaining.remove(direction)
                        coefficients[tuple(sorted((*mon,*remaining)))]+=c*d*other.count(direction)
        jacobi.extend(coefficients[mon]for mon in MON3)
    assert len(jacobi)==350
    return jacobi


def digest(values):return hashlib.sha256(json.dumps(list(map(str,values)),separators=(',',':')).encode()).hexdigest()


def review(model_path,output):
    started=time.monotonic();model_path=model_path.resolve()
    child=json.loads(gzip.decompress(model_path.read_bytes()))
    parent_path=ROOT/child['parent_model'];parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    metadata_path=model_path.parent/'report.json';elimination_path=model_path.parent/'elimination.json'
    adhoc_validation=model_path.parent/'validation.json'
    parent_review_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json'
    metadata=json.loads(metadata_path.read_text());parent_review=json.loads(parent_review_path.read_text())
    assert metadata['all_passed']and parent_review['all_passed']
    for record in [metadata,parent_review]:
        for name,want in record.get('input_sha256',{}).items():assert sha(ROOT/name)==want,name
    for name,want in metadata['generated_sha256'].items():assert sha(model_path.parent/name)==want,name
    assert sha(parent_path)==child['parent_archive_sha256']
    assert child['input_names']==parent['input_names'][:-1]
    assert child['input_weights']==parent['input_weights'][:-1]
    assert parent['input_names'][-1]=='z33'and parent['input_weights'][-1]==6
    assert child['exact_eliminating_equation_index']==296
    assert len(child['nodes'])==12569 and len(child['all350_residual_nodes'])==350
    assert len(child['raw150_reconstruction_nodes'])==150
    for key in ['group_cut_raw_indices','positive_raw_indices']:
        assert child[key]==parent[key]
    assert child['parent_free_positive_columns']==parent['free_positive_columns']

    localization=R.from_expr(s.sympify(child['reduced_a_open']))
    parent_localization=R.from_expr(s.sympify(parent['reduced_a_open']))
    assert localization==parent_localization*(9*a-2)
    target=parent['all350_residual_nodes'][296]
    polynomial,ancestor_count=sparse_expansion(parent['nodes'],target,34)
    assert len(polynomial)==273 and sum(sum(m)>1 for m in polynomial)==269
    eliminated_mon=(0,)*33+(1,)
    assert all(mon==eliminated_mon for mon in polynomial if mon[33])
    coefficient=polynomial[eliminated_mon]
    expected=K(-8*(9*a-2)*(9*a+4))/K(3*(27*a*a-60*a+8))
    assert coefficient==K(expected)
    saved=json.loads(elimination_path.read_text())
    recorded={tuple(mon):K.from_expr(s.sympify(c))for mon,c in saved['eliminating_residual_QQ_a_polynomial']}
    assert polynomial==recorded
    solution={mon[:-1]:-c/coefficient for mon,c in polynomial.items()if mon!=eliminated_mon}
    assert len(solution)==272
    assert all(sum(d*w for d,w in zip(mon,child['input_weights'][1:]))==6 for mon in solution)
    denominator_supported(coefficient.numer,localization)
    denominator_supported(coefficient.denom,localization)
    for c in solution.values():denominator_supported(c.denom,localization)
    print('Independent sparse full273 terms and272-term solution PASS',flush=True)

    parent_graph=audit_graph(parent,localization);child_graph=audit_graph(child,localization)
    assert child_graph['weights'][child['eliminated_parent_input_node']]==6
    for key in ['all350_residual_nodes','raw150_reconstruction_nodes']:
        assert [parent_graph['weights'][i]for i in parent[key]]==[child_graph['weights'][i]for i in child[key]]
    normal=ArithmeticNormalForm(33)
    replacement={tuple(i for i,power in enumerate(mon)for _ in range(power)):c for mon,c in solution.items()}
    inputs=[normal.scalar(K(a))]+[normal.variable(i)for i in range(33)]
    parent_values=normal.evaluate(parent['nodes'],inputs+[replacement])
    child_values=normal.evaluate(child['nodes'],inputs)
    assert child_values[child['eliminated_parent_input_node']]==replacement
    for key in ['all350_residual_nodes','raw150_reconstruction_nodes']:
        assert all(parent_values[i]==child_values[j]for i,j in zip(parent[key],child[key]))
    print('Independent universal350+150 arithmetic composition PASS',flush=True)

    parent_jets=linear_jets(parent['nodes'],34)
    linear=parent_jets[target][1];assert linear[33]==coefficient
    linear_rows=0
    for i in parent['all350_residual_nodes']:
        constant,row=parent_jets[i];assert not constant
        if row:linear_rows+=1
        ratio=row.get(33,K.zero)/coefficient
        assert all(row.get(j,K.zero)==ratio*linear.get(j,K.zero)for j in range(34))
    assert linear_rows==2
    child_jets=linear_jets(child['nodes'],33)
    assert all(child_jets[i]==(K.zero,{})for i in child['all350_residual_nodes'])
    assert s.cancel(coefficient.as_expr().subs(R.symbols[0],1))==s.Rational(728,75)

    samples=[]
    for number,parameter in enumerate([Fraction(1),Fraction(2),Fraction(3,7)]):
        positives=[Fraction((-1)**j*((j+1)*(number+1)+1),j+number+3)for j in range(33)]
        assert all(positives)
        child_fraction=evaluate_fraction(child['nodes'],[parameter]+positives)
        eliminated=child_fraction[child['eliminated_parent_input_node']]
        parent_fraction=evaluate_fraction(parent['nodes'],[parameter]+positives+[eliminated])
        raw=[child_fraction[i]for i in child['raw150_reconstruction_nodes']]
        residual=[child_fraction[i]for i in child['all350_residual_nodes']]
        assert raw==[parent_fraction[i]for i in parent['raw150_reconstruction_nodes']]
        assert residual==[parent_fraction[i]for i in parent['all350_residual_nodes']]
        assert residual[296]==0
        physical=physical_equations(raw);assert residual==physical
        selected=[i for block in parent['implicit_blocks']for i in block['equation_rows']]
        assert len(selected)==74 and all(residual[i]==0 for i in selected)
        nonzero=sum(bool(value)for value in residual);assert nonzero>0
        samples.append(dict(a=str(parameter),positive_inputs=list(map(str,positives)),
                            eliminated_z33=str(eliminated),all150_parent_reconstruction_equal=True,
                            physical_divergence25_zero=True,all350_parent_and_physical_Jacobi_equal=True,
                            parent74_selected_equations_and_f296_zero=True,nonzero_residual_outputs=nonzero,
                            raw150_sha256=digest(raw),residual350_sha256=digest(residual)))
    print('Three independent denseQQ physical350+divergence25 tests PASS',flush=True)
    note=ROOT/'research/notes/Q_NONRESONANT_FINITE_MODEL_20261103_v63_REVIEW.md'
    paths=[Path(__file__),note,model_path,metadata_path,elimination_path,adhoc_validation,
           parent_path,parent_review_path,ROOT/'code/q_nonresonant_finite_model_20261103_v63.py',
           ROOT/'research/notes/Q_NONRESONANT_FINITE_MODEL_20261103_v63.md']
    result=dict(id='Q_NONRESONANT_FINITE_MODEL_20261103_V63_REVIEW',all_passed=True,
                no_producer_or_DAG_implementation_imports=True,
                independently_expanded_actual_F296_monomials=273,independently_checked_nonlinear_terms=269,
                exact_273_coefficient_archive_match=True,independent_dependency_ancestor_count=ancestor_count,
                full272_term_solution_matches_actual_child_eliminated_node=True,
                exact_coefficient=str(s.factor(coefficient.as_expr())),coefficient_at_a_one='728/75',
                all_actual_inverse_dependencies_ground_a_only=True,
                all_actual_inverse_numerators_and_denominators_supported_on_Dprime=True,
                child_inverse_nodes_checked=child_graph['inverse_count'],
                child_distinct_inverse_functions=child_graph['distinct_inverse_functions'],
                all_node_topology_and_weight_homogeneity_verified=True,
                universal350_residual_and150_reconstruction_composition_equalities=True,
                composition_method='Independent exact arithmetic normal form overQQ(a), no positive-degree truncation or numeric inference.',
                normalized_additive_atoms=len(normal.atoms)-33,
                exact_quotient_inverse_isomorphism_accepted_on_Dprime=True,
                inherited_parent_smooth_factors_affine17_projective16_accepted=True,
                parent_free_residual_linear_rank=1,parent_free_nonzero_linear_rows=linear_rows,
                full_positive_linear_rank_derived_from_accepted74_plus1=75,
                all_child_positive_linear_jets_zero=True,Q_a_one_transverse_tangent_dimension=34,
                Q_a_one_full_projective_tangent_dimension=50,
                dense_QQ_physical_tests=samples,
                producer_adhoc_samples_not_used_as_universal_proof=True,
                full_branch_primary_or_component_classification_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
                seconds=time.monotonic()-started)
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS exact nonresonant quotient/composition/physical review',result['seconds'],'seconds',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'reports/Q_NONRESONANT_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();review(args.model,args.output)
