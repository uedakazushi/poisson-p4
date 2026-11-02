#!/usr/bin/env python3
"""Verify three original unframed 2814-column quartic minors from the135 equations."""
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
import json
import subprocess
import tempfile
import time

from s61_binary_hyperplane_cubic_matrix_20261102 import compile_native
from s61_binary_hyperplane_root_model_20261102 import (
    ROOT,CODE,P,S61,digest,hashes,old_cut_rows,original_point,output_argument,
    residue,runtime_environment,source_paths,write_new,
)

CERT = CODE/'certificates/s61_binary_hyperplane_original_root_minors_20261102.json'
NATIVE = CODE/'s61_binary_hyperplane_root_minor_20261102.cpp'
NOTE = ROOT/'research/notes/S61_BINARY_HYPERPLANE_ORIGINAL_ROOT_QUARTICS_20261102.md'


def products(weights):
    one,two,three,four = [[j for j,w in enumerate(weights) if w==k] for k in (1,2,3,4)]
    powers = lambda r: list(combinations_with_replacement(one,r))
    mon1 = [(j,) for j in one]
    mon2 = [(j,) for j in two]+powers(2)
    mon3 = ([(j,) for j in three]+[tuple(sorted((i,j))) for i in one for j in two]+powers(3))
    mon4 = ([(j,) for j in four]+[tuple(sorted((i,j))) for i in one for j in three]
            +list(combinations_with_replacement(two,2))
            +[tuple(sorted(m+(j,))) for m in powers(2) for j in two]+powers(4))
    assert len(mon4)==len(set(mon4))==2814
    return [ [()],mon1,mon2,mon3 ], {m:i for i,m in enumerate(sorted(mon4))}


def verify():
    start=time.time()
    cert=json.loads(CERT.read_text())
    assert cert['prime']==P and cert['cuts']==[5,19,23] and cert['e']==2 and cert['d']==3
    model=S61()
    result={}
    generated={}
    with tempfile.TemporaryDirectory(prefix='poisson-original-root4-',dir='/tmp') as tmp:
        wd=Path(tmp)
        exe=wd/'minor'
        compile_native(NATIVE,exe)
        for name,item in cert['points'].items():
            a=Fraction(item['a'])
            assert {'a':a,'F':500*a+342,'D0':50*a+36}[name]==0
            _,_,old,weights,rw,full=original_point(model,a)
            assert old_cut_rows(old)==cert['cuts']
            assert [j for j,w in enumerate(weights) if w==0]==[19,23]
            shifts,columns=products(weights)
            equations=[{m:residue(c) for m,c in q.items()
                        if 19 not in m and 23 not in m and residue(c)}
                       for q in full]
            equations.append({(5,):1})
            labels=[(i,shift) for i,w in enumerate(rw+[1]) if w<=4 for shift in shifts[4-w]]
            assert len(labels)==5089 and sum(i==134 for i,_ in labels)==12
            selected=item['rows']
            assert len(selected)==len(set(selected))==2814
            assert all(0<=k<len(labels) and labels[k][0]!=134 for k in selected)
            input_path=wd/(name+'.in')
            output=wd/(name+'.out')
            with input_path.open('w') as stream:
                stream.write(f'2814 {P}\n')
                for k in selected:
                    i,shift=labels[k]
                    entries=[(columns[tuple(sorted(m+shift))],v) for m,v in equations[i].items()]
                    stream.write(str(len(entries))+''.join(f' {j} {v}' for j,v in entries)+'\n')
            subprocess.run([str(exe),str(input_path),str(output)],check=True,
                           capture_output=True,text=True,env=runtime_environment(),timeout=75)
            determinant=int(output.read_text())
            assert determinant==item['determinant'] and determinant%P
            result[name]={'a':str(a),'a_modp':residue(a),'original_134_cut_rows':5077,
                          'all_135_cut_rows':5089,'columns':2814,'minor_det_modp':determinant,
                          'rank_134':2814,'rank_135':2814,'weight4_quotient_dimension':0}
            generated[name]={'minor_input':digest(input_path),'native_output':digest(output)}
    static=source_paths()+[Path(__file__),CERT,NATIVE,NOTE,CODE/'s61_binary_hyperplane_cubic_matrix_20261102.py']
    return {'all_passed':True,'id':'S61_BINARY_HYPERPLANE_ORIGINAL_ROOT_QUARTICS_20261102',
            'scope':'Three specified rational points at e=2,d=3 with the original common cuts.',
            'original_equations_reconstructed':135,'cuts':[5,19,23],'points':result,
            'input_sha256':hashes(static),'generated_sha256':generated,'seconds':time.time()-start}


if __name__=='__main__':
    path=output_argument(__doc__)
    report=verify()
    write_new(path,report)
    print('PASS',report['points'],report['seconds'])
