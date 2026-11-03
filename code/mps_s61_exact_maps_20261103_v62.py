#!/usr/bin/env python3
"""Freeze/verify general-parameter MPS maps in S61 with native exact arithmetic."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time

import sympy as s

from mps_s61_catalogue_20261103_v62 import raw_family,serialize,normalize_family,DIAGRAMS,ROOT
from mps_smoothing_completion_20261103_v62 import complete
from s61_rational_family_adapter_20261103_v62 import normalized_field
from s61_binary_zero_pair_stratum_removal_20261102 import compile_span
from s61_model import X,J,coefficients,jacobi_projective


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def write_native(path,model,K,parameters,normal):
    weights=model.reduced['parameter_weights'];D=K.one.denom
    for value,weight in zip(parameters,weights):
        if weight==1:D=D.lcm(value.denom)
    numerators=[]
    for value,weight in zip(parameters,weights):
        q,r=divmod(D**weight,value.denom);assert not r
        numerators.append(value.numer*q)
    with path.open('w') as stream:
        print(len(K.symbols),file=stream);print(*weights,file=stream)
        def poly(p):
            print(len(p),file=stream)
            for mon,c in sorted(p.items()):print(c,*mon,file=stream)
        def encoded(row):
            print(len(row),file=stream)
            for mon,num,den in row:print(len(mon),s.Rational(num,den),*mon,file=stream)
        poly(D)
        for p in numerators:poly(p)
        print(*coefficients(X[0]**2*J),file=stream)
        print(*model.certificate['raw_free_indices'],file=stream)
        print(len(model.data['eliminated_raw_variables']),file=stream)
        for item in model.data['eliminated_raw_variables']:
            print(item['raw_index'],file=stream);encoded(item['polynomial'])
        print(len(model.reduced['equations']),file=stream)
        for row in model.reduced['equations']:encoded(row)
        for value in normal:poly(value.numer);poly(value.denom)
    return dict(weight_one_common_denominator=str(D.as_expr()),
                weight_one_denominator_degree=max(map(sum,D.keys()),default=0),
                parameter_denominators_weight_compatible=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--families',nargs='+',default=['T03','T04','T05'])
    parser.add_argument('--native-timeout',type=int,default=180);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    source=ROOT/'code/s61_rational_map_verify_20261103_v62.cpp';executable=args.output/'verify-native'
    compile_span(source,executable)
    inputs=[Path(__file__),source,ROOT/'code/s61_rational_family_adapter_20261103_v62.py',
            ROOT/'code/mps_s61_catalogue_20261103_v62.py',ROOT/'code/mps_smoothing_completion_20261103_v62.py',
            ROOT/'code/s61_model.py',DIAGRAMS]
    inputs.extend(ROOT/'legacy/v02/data'/name for name in ['symplectic_slice_certificate.json',
                  'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    input_hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    diagrams={r['id']:r for r in json.loads(DIAGRAMS.read_text())['diagrams']};records=[]
    for name in args.families:
        started=time.monotonic();record=diagrams[name]
        periods,smoothings,P,pf,_=raw_family(record)
        if any(jacobi_projective(P)):
            periods,smoothings,P,pf,steps,_=complete(record)
        else:steps=[]
        variables=(*periods,*smoothings)
        model,K,values,normal,metadata=normalized_field(P,variables)
        print(name,'NORMALIZED',time.monotonic()-started,flush=True)
        native=args.output/f'{name}.in';common=write_native(native,model,K,values,normal)
        print(name,'EXPORTED',time.monotonic()-started,flush=True)
        log=args.output/f'{name}-native.log'
        with log.open('w') as stream:
            subprocess.run([str(executable.resolve()),str(native.resolve())],stdout=stream,stderr=subprocess.STDOUT,
                           check=True,timeout=args.native_timeout)
        assert 'PASS EQUATIONS 134 RECONSTRUCTION 150' in log.read_text()
        # The old expression normalizer is fast after rational specialization.
        sample=None
        for height in range(2,15):
            point={v:s.Integer(height)**(i+1) for i,v in enumerate(variables)}
            specialized=P.subs(point)
            if any(v.has(s.zoo,s.nan,s.oo,-s.oo) for v in specialized):continue
            try:_,old,_=normalize_family(specialized)
            except (StopIteration,ZeroDivisionError,ValueError):continue
            obtained=[s.cancel(K.to_sympy(v).subs(point)) for v in values]
            assert obtained==old['parameters'];sample={str(k):str(v) for k,v in point.items()};break
        assert sample is not None
        input_digest=sha(native)
        compressed=native.with_suffix('.in.gz')
        with compressed.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as archive:
                archive.write(native.read_bytes())
        native.unlink()
        denominators=sorted({str(s.factor(s.denom(c))) for i in range(5) for j in range(i+1,5)
                             for c in s.Poly(P[i,j],*X).coeffs() if s.denom(c)!=1})
        result=dict(component_id=name,all_passed=True,parameter_order=list(map(str,variables)),
                    period_pfaffian_open=str(pf),homogeneous_bivector=serialize(P),
                    coefficient_denominators_required_nonzero=denominators,
                    correction_steps=steps,normalization=metadata,common_denominator=common,
                    s61_coordinates_u1_through_u61=[str(v) for v in values],
                    exact_reduced_equations_checked=134,exact_reconstruction_coefficients_checked=150,
                    exact_projective_Jacobi_source_checked=True,
                    old_normalizer_exact_rational_cross_check=sample,
                    native_input_sha256=input_digest,native_input_archive=compressed.name,
                    native_input_archive_sha256=sha(compressed),native_log_sha256=sha(log),
                    input_sha256=input_hashes,
                    seconds=time.monotonic()-started,
                    scope='General-parameter coefficient chart in S61 on displayed pivot/Pfaffian opens. Residual Sp4 action supplies the full known component; density/label dependencies are recorded separately.')
        (args.output/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n')
        records.append(dict(component_id=name,seconds=result['seconds']))
        print(name,'FULL MAP PASS',result['seconds'],flush=True)
    executable.unlink()
    (args.output/'report.json').write_text(json.dumps(dict(id='MPS_S61_EXACT_MAPS_20261103_v62',all_passed=True,
         records=records,input_sha256=input_hashes,
         new_components_claimed=0,global_upper_bound_proved=False),indent=2)+'\n')
