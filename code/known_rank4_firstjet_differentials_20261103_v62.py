#!/usr/bin/env python3
"""Exact residue differential lower witnesses from the51 frozen rational maps.

Rank upper bounds come separately from component dimensions or the proved CA/CB
nine-dimensional incidence. A deficient sampled rank never proves an upper.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sympy as s
from mps_s61_dominance_20261103_v62 import rank_minor,PRIME
ROOT=Path(__file__).resolve().parents[1]

def differential(body):
    names=body['input_names'];active=[name for name in names if name!='r'];n=len(active)
    sample=body.get('exact_rational_sample',body.get('exact_sample'));values=[];gradients=[]
    def residue(value):
        if str(value)=='sqrt(3)':return int(s.sqrt_mod(3,PRIME))
        q=Fraction(str(value));return q.numerator*pow(q.denominator,-1,PRIME)%PRIME
    for node in body['nodes']:
        op=node[0];gradient=[0]*n
        if op=='q':value=residue(node[1])
        elif op=='v':
            name=names[node[1]];value=residue(sample[name])
            if name in active:gradient[active.index(name)]=1
        elif op=='+':
            value=sum(values[i] for i in node[1])%PRIME
            for i in node[1]:gradient=[(a+b)%PRIME for a,b in zip(gradient,gradients[i])]
        elif op=='*':
            ids=node[1];prefix=[1]
            for i in ids:prefix.append(prefix[-1]*values[i]%PRIME)
            value=prefix[-1];suffix=1
            for j in reversed(range(len(ids))):
                factor=prefix[j]*suffix%PRIME
                gradient=[(a+factor*b)%PRIME for a,b in zip(gradient,gradients[ids[j]])]
                suffix=suffix*values[ids[j]]%PRIME
        elif op=='inv':
            assert values[node[1]]%PRIME,'sample is outside the modular chart'
            value=pow(values[node[1]],-1,PRIME);gradient=[-value*value*v%PRIME for v in gradients[node[1]]]
        values.append(value);gradients.append(gradient)
    cert=json.loads((ROOT/'legacy/v02/data/symplectic_slice_certificate.json').read_text())
    indices=[i for i,item in enumerate(cert['free_variable_basis']) if item['weight']==1]
    matrix=s.Matrix([gradients[body['output_u1_through_u61'][i]] for i in indices])
    return active,indices,matrix,rank_minor(matrix)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);index=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    records=[]
    for item in json.loads(index.read_text())['records']:
        path=ROOT/item['map']['path'];body=json.loads(gzip.decompress(path.read_bytes()));active,indices,matrix,minor=differential(body)
        name=item['component_id'];upper=min(16,item['expected_s61_dimension']);dependency='S61 component dimension and dominance of the stored coordinate map'
        if name in ['R4_CA','R4_CB']:upper=9;dependency='research/notes/AMBIENT_WEIGHT1_IMAGE_20261030.md: closed line incidence of dimension9'
        assert minor['rank']<=upper
        record=dict(component_id=name,firstjet_indices_zero_based=indices,input_differential_columns=active,
            differential_matrix_mod32003=[[int(v) for v in row] for row in matrix.tolist()],nonzero_minor=minor,firstjet_image_dimension_lower_bound=minor['rank'],
            separate_geometric_upper_bound=upper,upper_dependency=dependency,
            exact_image_dimension=upper if upper==minor['rank'] else None,can_dominate_general_firstjet=(upper==16),
            circuit_archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            field_embedding_r_mod32003=int(s.sqrt_mod(3,PRIME)) if body.get('field_relation') else None)
        out=args.output/(name+'.json');out.write_text(json.dumps(record,indent=2)+'\n');records.append(dict(component_id=name,
            lower=minor['rank'],upper=upper,exact=upper==minor['rank'],sha256=hashlib.sha256(out.read_bytes()).hexdigest()))
        print(name,minor['rank'],upper,flush=True)
    (args.output/'report.json').write_text(json.dumps(dict(id='KNOWN51_RANK4_FIRSTJET_DIFFERENTIALS_20261103_v62',all_passed=True,
        records=records,general_firstjet_dominant_known_components=[r['component_id'] for r in records if r['lower']==16],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),index,ROOT/'code/mps_s61_dominance_20261103_v62.py',ROOT/'legacy/v02/data/symplectic_slice_certificate.json']},
        global_exhaustion_claimed=False),indent=2)+'\n')
