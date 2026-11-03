#!/usr/bin/env python3
"""Freeze completed native certificates separately from a timed-out T38 batch."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    source=ROOT/'reports/KNOWN_RANK4_FIRSTJET_NATIVE_CIRCUIT_20261103_v62'
    target=ROOT/'reports/MPS_FIRSTJET_NATIVE_RANKS_T13_T21_T30_20261103_v62';target.mkdir(exist_ok=False)
    records=[]
    for name in ['T13','T21','T30']:
        for suffix in ['.json','.in.gz','-native.log']:
            (target/(name+suffix)).write_bytes((source/(name+suffix)).read_bytes())
        body=json.loads((target/(name+'.json')).read_text());records.append(dict(component_id=name,
            rank=body['exact_generic_firstjet_rank'],sha256=hashlib.sha256((target/(name+'.json')).read_bytes()).hexdigest()))
    inputs=[Path(__file__),ROOT/'code/known_rank4_firstjet_native_circuit_20261103_v62.py',ROOT/'code/firstjet_circuit_function_rank_20261103_v62.cpp',
        ROOT/'legacy/v02/data/symplectic_slice_certificate.json',ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json']
    (target/'report.json').write_text(json.dumps(dict(id='MPS_FIRSTJET_NATIVE_RANKS_T13_T21_T30_20261103_v62',all_passed=True,records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        excluded_unfinished_computation='T38 in the original batch timed out after300seconds; no exact rank result for T38 is asserted here.'),indent=2)+'\n')
