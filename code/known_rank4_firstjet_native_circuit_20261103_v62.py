#!/usr/bin/env python3
"""Freeze native circuit AD/rank certificates without Python polynomial swelling."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time
ROOT=Path(__file__).resolve().parents[1]

def export(path,body):
    names=body['input_names'];groups=['g'+str(i) for i in range(10)];group=[1,0,0,1,0,0,0,0,0,0]
    variables=[n for n in names if n!='r' and n not in groups];active=[n for n in names if n!='r']
    cert=json.loads((ROOT/'legacy/v02/data/symplectic_slice_certificate.json').read_text())
    outputs=[body['output_u1_through_u61'][i] for i,b in enumerate(cert['free_variable_basis']) if b['weight']==1]
    with path.open('w') as out:
        out.write(f'{len(variables)} {len(names)} {len(active)} {len(body["nodes"])} 16\n')
        for name in names:
            kind,index=('r',0) if name=='r' else ('c',group[groups.index(name)]) if name in groups else ('v',variables.index(name))
            out.write(f'{kind} {index} {active.index(name) if name in active else -1}\n')
        for node in body['nodes']:
            op=node[0]
            if op in ['+','*']:out.write(f'{op} {len(node[1])} '+' '.join(map(str,node[1]))+'\n')
            else:out.write(f'{"i" if op=="inv" else op} {node[1]}\n')
        out.write(' '.join(map(str,outputs))+'\n')
    return variables,active

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--families',nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--binary',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);index=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    source={r['component_id']:r for r in json.loads(index.read_text())['records']};records=[]
    for name in args.families:
        start=time.monotonic();archive=ROOT/source[name]['map']['path'];body=json.loads(gzip.decompress(archive.read_bytes()))
        path=args.output/(name+'.in');variables,active=export(path,body);print(name,'NATIVE CIRCUIT READY',flush=True)
        result=subprocess.run([str(args.binary.resolve()),str(path.resolve())],capture_output=True,text=True,timeout=300,check=True)
        log=args.output/(name+'-native.log');log.write_text(result.stdout+result.stderr)
        rank=int(next(line for line in result.stdout.splitlines() if line.startswith('EXACT_CIRCUIT_FUNCTION_FIELD_RANK ')).split()[1])
        encoded=path.read_bytes();compressed=args.output/(name+'.in.gz')
        with compressed.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as output:output.write(encoded)
        payload=dict(component_id=name,exact_generic_firstjet_rank=rank,function_variables=variables,differential_column_order=active,
            coefficient_field='QQ(function_variables)[r]/(r^2-3)' if body.get('field_relation') else 'QQ(function_variables)',
            retained_sp4_lie_columns=True,residual_sp4_point='identity',
            circuit_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(encoded).hexdigest(),
            compressed_input_sha256=hashlib.sha256(compressed.read_bytes()).hexdigest(),native_log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
            seconds=time.monotonic()-start)
        report=args.output/(name+'.json');report.write_text(json.dumps(payload,indent=2)+'\n');path.unlink()
        records.append(dict(component_id=name,rank=rank,sha256=hashlib.sha256(report.read_bytes()).hexdigest(),seconds=payload['seconds']))
        print(name,'EXACT CIRCUIT RANK PASS',rank,payload['seconds'],flush=True)
    inputs=[Path(__file__),ROOT/'code/firstjet_circuit_function_rank_20261103_v62.cpp',ROOT/'legacy/v02/data/symplectic_slice_certificate.json',index]
    (args.output/'report.json').write_text(json.dumps(dict(id='KNOWN_RANK4_FIRSTJET_NATIVE_CIRCUIT_20261103_v62',records=records,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        proof_dependency='Primitive first-jet quotient is Sp4 equivariant; identity with all ten Lie columns; exact field AD and elimination.'),indent=2)+'\n')
