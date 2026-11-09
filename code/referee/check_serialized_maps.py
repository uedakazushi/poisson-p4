#!/usr/bin/env python3
"""Reexecute the reviewed circuit constructor and compare all serialized nodes.

This is serialization reproducibility using the original constructor, not an
independent proof of its algebra. The referee report separately reviews the
universal change-of-coordinates proof and tests raw Jacobi with another code.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'code'))
import sympy as s
from mps_s61_arithmetic_circuits_20261103_v62 import construct
from mps_s61_catalogue_20261103_v62 import raw_family, DIAGRAMS
from additional_rank4_families_20261103_v62 import families
from s61_model import X, pushforward


def compare(name, P, parameters, path):
    dag, outputs, metadata = construct({}, P, parameters)
    encoded = gzip.decompress(path.read_bytes())
    data = json.loads(encoded)
    assert json.loads(json.dumps(dag.nodes)) == data['nodes'], name
    assert outputs == data['output_u1_through_u61'], name
    assert dag.names == data['input_names'], name
    assert metadata == data['normalization'], name
    print(name, 'all', len(dag.nodes), 'nodes identical', flush=True)
    return {'id': name, 'nodes':len(dag.nodes), 'all_nodes_outputs_inputs_normalization_identical':True,
            'archive':str(path.relative_to(ROOT)), 'archive_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def main(output):
    rows = []
    diagrams = {r['id']:r for r in json.loads(DIAGRAMS.read_text())['diagrams']}
    for n in range(1,39):
        name = f'T{n:02d}'
        p,t,P,_,_ = raw_family(diagrams[name])
        completed = ROOT/f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/{name}.json'
        if completed.exists():
            data = json.loads(completed.read_text()); P = s.zeros(5)
            for key,expression in data['homogeneous_bivector'].items():
                i,j = map(int,key.split(',')); P[i,j] = s.sympify(expression); P[j,i] = -P[i,j]
        rows.append(compare(name,P,(*p,*t),ROOT/f'reports/MPS_S61_FULL_CIRCUITS_20261103_v62/{name}.json.gz'))
    for record in families():
        name = record['component_id']; point = s.symbols('base1:5'); scale = s.Symbol('scale')
        parameters = [*record['parameters'],*([s.Symbol('r')] if record['field_relation'] is not None else []),scale,*point]
        translation = s.eye(5)
        for i,p in enumerate(point,1): translation[i,0] = 1-p
        P = pushforward(scale*record['P'],translation)
        rows.append(compare(name,P,parameters,ROOT/f'reports/ADDITIONAL_RANK4_S61_CIRCUITS_20261103_v62/{name}.json.gz'))
    output.write_text(json.dumps({'scope':__doc__,'families':rows,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
