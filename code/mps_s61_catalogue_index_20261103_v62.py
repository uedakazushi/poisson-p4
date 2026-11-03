#!/usr/bin/env python3
"""Freeze the forty known MPS component maps and their exact provenance."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def reference(path):return dict(path=str(path.relative_to(ROOT)),sha256=digest(path))

if __name__=='__main__':
    diagrams=json.loads((ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json').read_text())['diagrams']
    records=[]
    for source in diagrams:
        name=source['id'];number=int(name[1:]);graph=number>=39
        folder='MPS_ELLIPTIC_GRAPH_CIRCUITS_20261103_v62' if graph else 'MPS_S61_FULL_CIRCUITS_20261103_v62'
        archive=ROOT/'reports'/folder/f'{name}.json.gz';encoded=gzip.decompress(archive.read_bytes());body=json.loads(encoded)
        record=dict(component_id=name,map=reference(archive),decoded_sha256=hashlib.sha256(encoded).hexdigest(),
                    coordinate_order='u1,...,u61, using code/s61_model.py',input_names=body['input_names'],
                    map_kind='global_quadratic_graph' if graph else 'rational_period_edge_Sp4_arithmetic_circuit',
                    expected_s61_dimension=16 if graph else source['projective_component_dimension']-9,
                    source_diagram=source,full_s61_component_dominance=True,defining_prime_ideal_constructed=False)
        if graph:
            report=ROOT/body['source_report'];record['prior_global_graph_proof']=reference(report)
            record['opens']='none: the graph is polynomial on the entire sixteen-dimensional first-jet space'
        else:
            record.update(period_dimension=source['period_dimension'],edge_count=source['edge_count'],
                          differential_certificate=reference(ROOT/'reports/MPS_S61_DOMINANCE_20261103_v62'/f'{name}.json'),
                          opens=dict(period_pfaffian=body['source_period_pfaffian_open'],
                                     coefficient_denominators=body['source_coefficient_denominators_required_nonzero'],
                                     value_pivot_node=body['normalization']['value_pivot_node'],
                                     value_pfaffian_node=body['normalization']['value_pfaffian_node'],
                                     sp4_A_determinant_node=body['normalization']['sp4_big_cell_determinant_node']),
                          generic_identity_dependency='Generic source Jacobi, rational coordinate pushforwards, Euler quotient, and the proved S61 reduction theorem',
                          component_identification_dependency='MPS Theorem6.16, Corollary6.23, Example6.24 and the exact decorated-edge crosswalk')
            completion=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62'/f'{name}.json'
            if completion.exists():record['completed_source_tensor']=reference(completion)
        records.append(record)
    assert len(records)==40 and len({r['component_id'] for r in records})==40
    output=ROOT/'code/data/mps_s61_coordinate_catalogue_20261103_v62.json'
    output.write_text(json.dumps(dict(id='MPS40_S61_COORDINATE_CATALOGUE_20261103_v62',records=records,
        provenance_index_builder=reference(Path(__file__)),proof_note='research/notes/MPS_S61_CATALOGUE_20261103_v62.md',
        global_exhaustion_claimed=False,known_component_lower_bound_unchanged=57),indent=2)+'\n')
    print(str(output.relative_to(ROOT)),digest(output),'40 maps')
