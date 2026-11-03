#!/usr/bin/env python3
"""Full S61 rational coordinate circuits for the eleven additional components.

The generic point of the original family is retained as four affine variables;
this replaces the MPS-specific use of its diagonal torus. Scalar and Sp4 are
explicit. Prior component theorems identify the dense PGL saturation.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s

from additional_rank4_families_20261103_v62 import families,ROOT
from mps_s61_arithmetic_circuits_20261103_v62 import construct
from s61_model import S61,X,PAIRS,pushforward

def evaluate_algebraic(dag,point):
    K=s.QQ.algebraic_field(s.sqrt(3));values=[]
    for node in dag.nodes:
        op=node[0]
        if op=='q':value=K.from_sympy(s.Rational(node[1]))
        elif op=='v':value=K.from_sympy(s.sympify(point[dag.names[node[1]]]))
        elif op=='+':value=sum((values[i] for i in node[1]),K.zero)
        elif op=='*':
            value=K.one
            for i in node[1]:value*=values[i]
        elif op=='inv':value=K.one/values[node[1]]
        values.append(value)
    return [K.to_sympy(v) for v in values]

def verify_sample(dag,outputs,metadata,Q,parameters,algebraic):
    point={str(p):s.Integer((i%5)+2) for i,p in enumerate(parameters)}
    # The coefficient-field sample has the exact defining algebraic relation.
    # At its zero point the value is J, making an independent K-check bounded.
    if algebraic:
        point['r']=s.sqrt(3)
        for name in ['base1','base2','base3','base4']:point[name]=0
        point['scale']=1
    point.update({f'g{i}':v for i,v in enumerate([2,1,1,1,1,2,3,2,1,4])})
    values=evaluate_algebraic(dag,point) if algebraic else [s.Rational(v.numerator,v.denominator) for v in dag.evaluate(point)]
    matrix=s.Matrix([[values[j] for j in row] for row in metadata['linear_normalizing_matrix_nodes']])
    source=Q.subs({p:point[str(p)] for p in parameters});old=S61();normalized=old.normalize(pushforward(source,matrix))
    assert all(s.cancel(values[j]-v)==0 for j,v in zip(outputs,normalized['parameters']))
    checks=old.verify_point(normalized)
    return {k:str(v) for k,v in point.items()},checks

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--families',nargs='+');args=parser.parse_args();args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False)
    records=[]
    for record in families():
        name=record['component_id']
        if args.families and name not in args.families:continue
        start=time.monotonic();point=s.symbols('base1:5');scale=s.Symbol('scale');r=s.Symbol('r')
        source_parameters=[*record['parameters'],*([r] if record['field_relation'] is not None else []),scale,*point]
        translation=s.eye(5)
        for i,p in enumerate(point,1):translation[i,0]=1-p
        Q=pushforward(scale*record['P'],translation)
        dag,outputs,metadata=construct(record,Q,source_parameters)
        sample,checks=verify_sample(dag,outputs,metadata,Q,source_parameters,record['field_relation'] is not None)
        payload=dict(component_id=name,input_names=dag.names,nodes=dag.nodes,output_u1_through_u61=outputs,
            original_family_translation='Q=(translation(1-base))_*(scale*P); the universal fixed-one normalizer then evaluates P at(1:base1:base2:base3:base4)',
            normalization=metadata,field_relation=str(record['field_relation']) if record['field_relation'] is not None else None,
            coefficient_field='QQ[r]/(r^2-3)' if record['field_relation'] is not None else 'QQ',
            residual_sp4_included=True,overall_scalar_retained=True,full_source_dimension=len(source_parameters)+10-(record['field_relation'] is not None),
            expected_s61_dimension=record['projective_dimension']-9,
            generic_identity_proof_dependency='Exact generic family Jacobi; finite rational projective pushforwards, Euler quotient, perspective gauge, and existing S61 reduction theorem. Exact sample is implementation verification only.',
            component_identification_dependency=dict(path=record['source'],anchor=record['source_anchor']),
            exact_sample=sample,sample_checks=checks,seconds=time.monotonic()-start)
        encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode();path=args.output/(name+'.json.gz')
        with path.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as archive:archive.write(encoded)
        records.append(dict(component_id=name,path=str(path.relative_to(ROOT)),archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            circuit_sha256=hashlib.sha256(encoded).hexdigest(),nodes=len(dag.nodes),seconds=payload['seconds'],exact_sample_passed=True))
        print(name,'FULL S61 CIRCUIT PASS',len(dag.nodes),payload['seconds'],flush=True)
    sources=[Path(__file__),ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',ROOT/'code/s61_model.py']
    (args.output/'report.json').write_text(json.dumps(dict(id='ADDITIONAL_RANK4_S61_CIRCUITS_20261103_v62',all_passed=True,records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        new_components_claimed=0,global_exhaustion_claimed=False),indent=2)+'\n')
