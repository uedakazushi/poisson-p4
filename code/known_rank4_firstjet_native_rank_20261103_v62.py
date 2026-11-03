#!/usr/bin/env python3
"""Export exact first-jet function matrices and use native FLINT rational rank."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from known_rank4_firstjet_function_field_20261103_v62 import exact_matrix,ROOT

def write_matrix(path,K,matrix):
    with path.open('w') as stream:
        stream.write(f'{len(K.gens)} {len(matrix)} {len(matrix[0])}\n')
        for row in matrix:
            for value in row:
                polys=[value.numer,value.denom];denominators=[int(s.Rational(c).q) for p in polys for c in p.values()]
                common=s.ilcm(*denominators) if len(denominators)>1 else denominators[0]
                for poly in polys:
                    stream.write(str(len(poly))+'\n')
                    for powers,coefficient in sorted(poly.items()):
                        stream.write(str(int(common*s.Rational(coefficient)))+' '+' '.join(map(str,powers))+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--families',nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--binary',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);index=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    source={r['component_id']:r for r in json.loads(index.read_text())['records']};records=[]
    for name in args.families:
        started=time.monotonic();archive=ROOT/source[name]['map']['path'];body=json.loads(gzip.decompress(archive.read_bytes()))
        assert not body.get('field_relation'),'native pilot presently QQ only'
        K,active,variables,matrix,nodes=exact_matrix(body);input_path=args.output/(name+'.in');write_matrix(input_path,K,matrix)
        print(name,'NATIVE MATRIX READY',nodes,time.monotonic()-started,flush=True)
        result=subprocess.run([str(args.binary.resolve()),str(input_path.resolve())],capture_output=True,text=True,timeout=300,check=True)
        log=args.output/(name+'-native.log');log.write_text(result.stdout+result.stderr)
        line=next(v for v in result.stdout.splitlines() if v.startswith('EXACT_FUNCTION_FIELD_RANK '));rank=int(line.split()[1])
        encoded=input_path.read_bytes();compressed=args.output/(name+'.in.gz')
        with compressed.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as out:out.write(encoded)
        payload=dict(component_id=name,exact_generic_firstjet_rank=rank,coefficient_field=str(K),function_variables=variables,
            differential_column_order=active,retained_sp4_lie_columns=True,residual_sp4_point='identity',needed_circuit_nodes=nodes,
            circuit_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(encoded).hexdigest(),
            compressed_input_sha256=hashlib.sha256(compressed.read_bytes()).hexdigest(),native_log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
            seconds=time.monotonic()-started)
        path=args.output/(name+'.json');path.write_text(json.dumps(payload,indent=2)+'\n');input_path.unlink()
        records.append(dict(component_id=name,rank=rank,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=payload['seconds']))
        print(name,'EXACT NATIVE FUNCTION RANK PASS',rank,payload['seconds'],flush=True)
    inputs=[Path(__file__),ROOT/'code/known_rank4_firstjet_function_field_20261103_v62.py',ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp',index]
    (args.output/'report.json').write_text(json.dumps(dict(id='KNOWN_RANK4_FIRSTJET_NATIVE_RANK_20261103_v62',records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        proof_dependency='Sp4 equivariance of primitive first-jet quotient; exact field elimination, no sampled upper'),indent=2)+'\n')
