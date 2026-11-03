#!/usr/bin/env python3
"""Additional honest source specializations for deficient first-jet witnesses.

Only a maximum lower bound is reported; failure to improve is not an upper.
"""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import random
from known_rank4_firstjet_differentials_20261103_v62 import differential,ROOT

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);index=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    base=ROOT/'reports/KNOWN51_RANK4_FIRSTJET_DIFFERENTIALS_20261103_v62';records=[];rng=random.Random(20261103)
    for item in json.loads(index.read_text())['records']:
        name=item['component_id'];prior=json.loads((base/(name+'.json')).read_text())
        if prior['firstjet_image_dimension_lower_bound']==prior['separate_geometric_upper_bound']:continue
        body=json.loads(gzip.decompress((ROOT/item['map']['path']).read_bytes()));best=prior['firstjet_image_dimension_lower_bound'];trials=[]
        for attempt in range(4):
            sample={n:('sqrt(3)' if n=='r' else str(rng.randint(1,127))) for n in body['input_names']}
            modified=copy.copy(body);modified['exact_rational_sample']=sample
            try:active,indices,matrix,minor=differential(modified)
            except AssertionError:trials.append(dict(sample=sample,chart_defined_mod32003=False));continue
            best=max(best,minor['rank']);trials.append(dict(sample=sample,chart_defined_mod32003=True,minor=minor,
                differential_matrix_mod32003=[[int(v) for v in row] for row in matrix.tolist()]))
            if best==prior['separate_geometric_upper_bound']:break
        payload=dict(component_id=name,prior_lower_bound=prior['firstjet_image_dimension_lower_bound'],new_lower_bound=best,
            geometric_upper_bound=prior['separate_geometric_upper_bound'],trials=trials,rank_upper_inferred_from_sampling=False)
        out=args.output/(name+'.json');out.write_text(json.dumps(payload,indent=2)+'\n')
        records.append(dict(component_id=name,lower=best,upper=prior['separate_geometric_upper_bound'],sha256=hashlib.sha256(out.read_bytes()).hexdigest()))
        print(name,prior['firstjet_image_dimension_lower_bound'],'->',best,'<=',prior['separate_geometric_upper_bound'],flush=True)
    (args.output/'report.json').write_text(json.dumps(dict(id='KNOWN51_RANK4_FIRSTJET_PROBES_20261103_v62',records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),index,base/'report.json',ROOT/'code/known_rank4_firstjet_differentials_20261103_v62.py']},
        no_sampled_upper_bound=True),indent=2)+'\n')
