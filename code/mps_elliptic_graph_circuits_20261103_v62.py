#!/usr/bin/env python3
"""Place the existing, proved T39/T40 global quadratic graphs in the catalogue."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import sympy as s

from mps_s61_arithmetic_circuits_20261103_v62 import DAG,ROOT
from s61_model import S61,coefficients


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);records=[];inputs=[Path(__file__)]
    for component,stem in [('T39','Q51_UNIVERSAL_SECTION_20261101'),('T40','Q52_UNIVERSAL_SECTION_20261101')]:
        source=ROOT/'reports'/f'{stem}.json';data=json.loads(source.read_text());assert data['all_passed']
        inputs.append(source);one=data['weight1_u_indices_zero_based'];two=data['weight2_u_indices_zero_based']
        names=[f'u{i+1}' for i in one];dag=DAG(names);outputs=[dag.zero]*61
        for i,name in zip(one,names):outputs[i]=dag.inputs[name]
        for i,value in zip(two,data['quadratic_section_q']):outputs[i]=dag.expression(s.sympify(value))
        point={name:j+1 for j,name in enumerate(names)};evaluated=dag.evaluate(point)
        values=[s.Rational(evaluated[i].numerator,evaluated[i].denominator) for i in outputs]
        model=S61();P=model.reconstruct(values)
        checks=model.verify_point(dict(parameters=values,coefficients=coefficients(P)))
        payload=dict(component_id=component,map_kind='existing_global_quadratic_graph',
                     input_names=names,nodes=dag.nodes,output_u1_through_u61=outputs,
                     source_report=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                     coefficient_field='QQ('+','.join(names)+')',open_condition='none',
                     full_s61_component_covered=True,residual_sp4_needed=False,
                     residual_sp4_explanation='The16 first-jet coordinates already parametrize the entire known graph; no residual-group saturation is omitted.',
                     general_134_identity_dependency='Existing universal graph report and independent review; this repackaging does not replace those proofs.',
                     exact_rational_sample=point,sample_checks=checks)
        encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode();path=args.output/f'{component}.json.gz'
        with path.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as archive:archive.write(encoded)
        records.append(dict(component_id=component,nodes=len(dag.nodes),
                            circuit_sha256=hashlib.sha256(encoded).hexdigest(),archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            exact_sample_checked=True))
        print(component,'GLOBAL GRAPH CIRCUIT PASS',len(dag.nodes),flush=True)
    inputs.extend([ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',ROOT/'code/s61_model.py',
                   ROOT/'reports/Q51_UNIVERSAL_SECTION_REVIEW_20261101.json',ROOT/'reports/Q52_UNIVERSAL_SECTION_REVIEW_20261101.json'])
    (args.output/'report.json').write_text(json.dumps(dict(id='MPS_ELLIPTIC_GRAPH_CIRCUITS_20261103_v62',all_passed=True,
        records=records,new_components_claimed=0,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}),indent=2)+'\n')
