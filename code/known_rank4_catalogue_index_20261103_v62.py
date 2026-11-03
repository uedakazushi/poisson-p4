#!/usr/bin/env python3
"""Combine all51 already-known generic-rank4 S61 coordinate maps."""
import gzip
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def ref(path):return dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
if __name__=='__main__':
    mps=ROOT/'code/data/mps_s61_coordinate_catalogue_20261103_v62.json'
    records=json.loads(mps.read_text())['records']
    report=ROOT/'reports/ADDITIONAL_RANK4_S61_CIRCUITS_20261103_v62/report.json'
    for item in json.loads(report.read_text())['records']:
        name=item['component_id'];path=ROOT/item['path'];encoded=gzip.decompress(path.read_bytes());body=json.loads(encoded)
        source=ROOT/'reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62'/f'{name}.json'
        records.append(dict(component_id=name,map=ref(path),decoded_sha256=hashlib.sha256(encoded).hexdigest(),
            coordinate_order='u1,...,u61, using code/s61_model.py',input_names=body['input_names'],
            map_kind='rational_family_point_scalar_Sp4_arithmetic_circuit',coefficient_field=body['coefficient_field'],
            field_relation=body['field_relation'],expected_s61_dimension=body['expected_s61_dimension'],
            full_source_dimension=body['full_source_dimension'],full_s61_component_dominance=True,
            generic_source_and_image_minor=ref(source),component_identification_dependency=body['component_identification_dependency'],
            opens=dict(value_pivot_node=body['normalization']['value_pivot_node'],
                       value_pfaffian_node=body['normalization']['value_pfaffian_node'],
                       sp4_A_determinant_node=body['normalization']['sp4_big_cell_determinant_node'],overall_scalar='scale != 0'),
            defining_prime_ideal_constructed=False))
    assert len(records)==51 and len({r['component_id'] for r in records})==51
    body=dict(id='KNOWN51_RANK4_S61_COORDINATE_CATALOGUE_20261103_v62',records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),mps,report]},
        proof_notes=['research/notes/MPS_S61_CATALOGUE_20261103_v62.md','research/notes/KNOWN51_RANK4_COORDINATE_MAPS_20261103_v62.md'],
        defining_prime_ideals_completed=False,global_exhaustion_claimed=False,known_total_component_lower_bound_unchanged=57)
    output=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json';output.write_text(json.dumps(body,indent=2)+'\n')
    print(output.relative_to(ROOT),hashlib.sha256(output.read_bytes()).hexdigest(),'51 maps')
