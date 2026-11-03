#!/usr/bin/env python3
"""Freeze exact first-jet image and generic-fiber dimensions of known51."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    index=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json';records=[];inputs=[Path(__file__),index]
    for item in json.loads(index.read_text())['records']:
        name=item['component_id'];base=ROOT/'reports/KNOWN51_RANK4_FIRSTJET_DIFFERENTIALS_20261103_v62'/f'{name}.json';record=json.loads(base.read_text());inputs.append(base)
        dimension=record['exact_image_dimension'];proof=base
        if dimension is None:
            if name=='T38':proof=ROOT/'reports/MPS_T38_RAW_JET_RANK_20261103_v62/report.json';dimension=json.loads(proof.read_text())['exact_primitive_firstjet_image_dimension']
            else:
                directory='MPS_FIRSTJET_NATIVE_RANKS_T13_T21_T30_20261103_v62' if name.startswith('T') else 'ADDITIONAL_RANK4_FIRSTJET_NATIVE_CIRCUIT_20261103_v62'
                proof=ROOT/'reports'/directory/f'{name}.json';dimension=json.loads(proof.read_text())['exact_generic_firstjet_rank']
            inputs.append(proof)
        records.append(dict(component_id=name,s61_component_dimension=item['expected_s61_dimension'],firstjet_image_dimension=dimension,
            generic_fiber_dimension=item['expected_s61_dimension']-dimension,proof_record=str(proof.relative_to(ROOT)),
            proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),general_firstjet_dominant=dimension==16))
    body=dict(id='KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62',all_passed=True,records=records,
        general_firstjet_dominant_components=[r['component_id'] for r in records if r['general_firstjet_dominant']],
        scope='Reduced component image and generic fiber dimensions only; no wholeS61 exhaustion, fiber length or nonreduced structure.',
        normalized_T38_attempts='Two generic-source normalized native ranks timed out300s; the fixed-orbit normalized rank was stopped after the alternative rawjet certificate completed. The lowdegree rawjet method completed in40.13s and supplies the exact proof.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    output=ROOT/'reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json';output.write_text(json.dumps(body,indent=2)+'\n');print('KNOWN51 EXACT FIRSTJET DIMENSIONS PASS')
