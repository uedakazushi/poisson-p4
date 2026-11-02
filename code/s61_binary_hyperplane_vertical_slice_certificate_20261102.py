#!/usr/bin/env python3
"""Exact degree of a vertical slice minor and its two modular row-module controls."""
from pathlib import Path
import argparse
import json
import subprocess
import tempfile
import time

import numpy as np
import sympy as sp

from s61_binary_deeper_uniform_quadratic_20261101 import determinant_mod
from s61_binary_hyperplane_Bfree_horizontal_20261102 import compile_native
from s61_binary_hyperplane_root_model_20261102 import (
    ROOT,CODE,P,digest,hashes,runtime_environment,write_new,
)
from s61_binary_hyperplane_vertical_slice_data_20261102 import (
    at_two,coefficient_line,obtain_core,prepare,
)

CERT=CODE/'certificates/s61_binary_hyperplane_vertical_slice_20261102.json'
NOTE=ROOT/'research/notes/S61_BINARY_HYPERPLANE_VERTICAL_SLICE_CERTIFICATE_20261102.md'
NATIVES={
    'kernel':CODE/'s61_binary_hyperplane_leading_kernel_20261102.cpp',
    'peel':CODE/'s61_binary_hyperplane_vertical_rectangular_peel_20261102.cpp',
    'determinant':CODE/'s61_binary_hyperplane_vertical_minor_det_20261102.cpp',
    'module':CODE/'s61_binary_hyperplane_vertical_rowmodule_20261102.cpp',
}


def make_kernel_input(leading,output,selection):
    rows,cols,free=[selection[k] for k in ('rows','columns','free_columns')]
    assert len(rows)==len(cols)==len(set(rows))==len(set(cols))
    assert sorted(cols+free)==list(range(238))
    with leading.open() as stream:
        assert stream.readline().split()==['321','238']
        lines=stream.readlines()
    assert len(lines)==321*238
    values=np.array([at_two(coefficient_line(line)) for line in lines],
                    dtype=np.int64).reshape(321,238)
    minor=int(determinant_mod(values[np.ix_(rows,cols)],P))
    assert minor==selection['minor_modp'] and minor
    with output.open('w') as stream:
        stream.write(f'321 238 {len(cols)} {len(free)}\n')
        for indices in (rows,cols,free):
            stream.write(' '.join(map(str,indices))+'\n')
        stream.writelines(lines)
    return values


def verify_kernel_output(path,free,leading_at_two):
    values=np.zeros((238,len(free)),dtype=np.int64)
    anchor_values=[]
    terms=0
    max_degree=0
    with path.open() as stream:
        for k,anchor in enumerate(free):
            assert stream.readline().split()==['KERNEL',str(k)]
            denominator=coefficient_line(stream.readline())
            anchor_values.append(at_two(denominator))
            assert anchor_values[-1]
            for j in range(238):
                poly=coefficient_line(stream.readline())
                if j in free:
                    assert poly==(denominator if j==anchor else {})
                values[j,k]=at_two(poly)
                terms+=len(poly)
                max_degree=max(max_degree,max(poly,default=0))
        assert not stream.read().strip()
    assert not np.any(leading_at_two@values%P)
    return {'count':len(free),'numerator_terms':terms,'maximum_e_degree':max_degree,
            'anchor_values_at_e2':anchor_values,'full_QQ_identities_checked_by_native':321*len(free),
            'independent_modular_product_zero':True,'free_block_diagonal_verified':True}


def direct_minor(path):
    data={}
    evaluations=[]
    for line in path.read_text().splitlines():
        z=line.split()
        if z[0]=='DETERMINANT_DEGREE':
            data.update(degree=int(z[1]),leading_coefficient=int(z[3]))
        elif z[0]=='DETERMINANT_COEFFICIENTS':
            assert len(z)==int(z[1])+2
            data['coefficients']=list(map(int,z[2:]))
        elif z[0]=='EVALUATION':
            evaluations.append([int(z[1]),int(z[2])])
        elif z[0]!='TOTAL_SECONDS':
            raise ValueError(line)
    assert [v[0] for v in evaluations]==[0,3,17]
    data['independent_evaluations']=evaluations
    return data


def polynomial(coordinates):
    return sp.Poly.from_list(list(reversed(coordinates)),sp.Symbol('a'),modulus=P)


def module_output(path,expected):
    a=sp.Symbol('a')
    factors={k:sp.Poly(v,a,modulus=P).monic() for k,v in
             {'a':a,'F':500*a+342,'D0':50*a+36,'C':100*a+52}.items()}
    orders={}
    data={}
    for line in path.read_text().splitlines():
        z=line.split()
        if z[0]=='ROWMODULE':
            data.update(rank=int(z[1]),degree=int(z[3]),row_operations=int(z[5]))
        elif z[0]=='DETERMINANT_COEFFICIENTS':
            assert len(z)==int(z[1])+2
            det=polynomial(list(map(int,z[2:]))).monic()
        elif z[0]=='RESIDUAL_SIZE':
            data.update(residual_size=int(z[1]),factor_count=int(z[3]))
        elif z[0]=='FACTOR':
            assert len(z)==int(z[2])+3
            f=polynomial(list(map(int,z[3:]))).monic()
            matching=[k for k,val in factors.items() if val==f]
            assert len(matching)==1 and matching[0] not in orders
            orders[matching[0]]=int(z[1])
        else:
            raise ValueError(line)
    assert orders==expected and data['rank']==238 and data['degree']==287
    product=sp.Poly(1,a,modulus=P)
    for name,count in orders.items():
        product*=factors[name]**count
    assert product==det and det.degree()==287
    data.update(known_orders=orders,monic_coefficients=[int(det.nth(i))%P for i in range(288)])
    return data


def verify(wd,core_path=None):
    start=time.time()
    cert=json.loads(CERT.read_text())
    assert cert['prime']==P and cert['fixed_d_for_degree']==3 and cert['fixed_e_modp']==2
    core,static,source=obtain_core(wd,core_path,cert['core_sha256'])
    full,leading,field_matrices,prep=prepare(core,wd,cert)
    env=runtime_environment()
    times={}
    logs={}
    generated={q.name:digest(q) for q in [full,leading,*field_matrices.values()]}
    for name,src in NATIVES.items():
        compile_native(src,wd/name)
    def run(name,args,label):
        t=time.time()
        cp=subprocess.run([str(wd/name),*map(str,args)],check=True,capture_output=True,
                          text=True,env=env,timeout=125)
        times[label]=time.time()-t
        logs[label]=cp.stderr.splitlines()[-1] if cp.stderr else ''
        print(label,round(times[label],2),flush=True)
    audits=[]
    for stage,selection in enumerate(cert['kernel_pivots']):
        input_path=wd/f'kernel_input{stage}'
        kernel=wd/f'kernel_output{stage}'
        finite=make_kernel_input(leading,input_path,selection)
        run('kernel',[input_path,kernel],f'kernel_{stage}')
        audit=verify_kernel_output(kernel,selection['free_columns'],finite)
        audit['lower_minor_modp']=selection['minor_modp']
        audits.append(audit)
        generated[input_path.name]=digest(input_path)
        generated[kernel.name]=digest(kernel)
        if stage<2:
            next_full=wd/f'full{stage+1}'
            next_leading=wd/f'leading{stage+1}'
            run('peel',[full,input_path,kernel,next_full,next_leading,2-stage],f'peel_{stage}')
            full,leading=next_full,next_leading
            generated[full.name]=digest(full)
            generated[leading.name]=digest(leading)
    coranks=[q['count'] for q in audits]
    assert coranks==[76,33,10]
    rows=cert['selected_original_rows']
    assert len(rows)==len(set(rows))==238 and all(0<=j<321 for j in rows)
    row_sum=sum(cert['row_bounds'][j] for j in rows)
    col_sum=sum(cert['column_bounds'])
    upper=row_sum+col_sum
    assert upper==cert['selected_degree_upper_before_peels']==694
    selected=wd/'selected'
    selected.write_text('238\n'+' '.join(map(str,rows))+'\n')
    det_path=wd/'determinant_output'
    run('determinant',[field_matrices[3],selected,det_path],'direct_original_minor')
    det=direct_minor(det_path)
    assert upper-sum(coranks)==det['degree']==cert['selected_modular_degree']==575
    assert det['leading_coefficient']==cert['selected_modular_lc']==16008
    generated[det_path.name]=digest(det_path)
    modules={}
    for d in (3,4):
        output=wd/f'module_d{d}'
        run('module',[field_matrices[d],output],f'full_module_d{d}')
        modules[str(d)]=module_output(output,cert['rowmodule_known_orders'])
        generated[output.name]=digest(output)
    assert modules['3']['monic_coefficients']==modules['4']['monic_coefficients']
    static.update(hashes([Path(__file__),CERT,NOTE,*NATIVES.values(),
        CODE/'s61_binary_hyperplane_vertical_slice_data_20261102.py',
        CODE/'s61_binary_hyperplane_root_model_20261102.py',
        CODE/'s61_binary_deeper_cubic_field_20261101.hpp',
        CODE/'s61_binary_deeper_uniform_quadratic_20261101.py',
        CODE/'s61_binary_hyperplane_cubic_matrix_20261102.py']))
    return {'all_passed':True,'id':'S61_BINARY_HYPERPLANE_VERTICAL_SLICE_CERTIFICATE_20261102',
            'scope':'One original 238-row minor at d=3 over QQ(e), and full modular polynomial row modules at e=2,d=3,4.',
            'source':source,'prime':P,'preparation':prep,'kernel_dimensions':coranks,
            'leading_ranks_QQ_e':[238-k for k in coranks],'kernel_audits':audits,
            'delta_valuation_lower_all_maximal_minors':sum(coranks),'selected_original_rows':rows,
            'selected_row_weight_sum':row_sum,'column_weight_sum':col_sum,'selected_degree_upper':upper,
            'selected_minor_exact_a_degree_QQ_e':det['degree'],'direct_modular_minor':det,
            'full_modular_row_modules':modules,'input_sha256':static,'generated_sha256':generated,
            'native_stage_seconds':times,'native_log_tails':logs,'seconds':time.time()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--core-matrix',type=Path,help='Use the SHA-verified original BFREE_HORIZONTAL matrix.')
    args=parser.parse_args()
    if args.output.exists():
        parser.error(f'output already exists: {args.output}')
    with tempfile.TemporaryDirectory(prefix='poisson-vertical-slice-',dir='/tmp') as tmp:
        value=verify(Path(tmp),args.core_matrix)
    write_new(args.output,value)
    print('PASS',value['kernel_dimensions'],value['selected_minor_exact_a_degree_QQ_e'],value['seconds'])
