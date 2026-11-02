#!/usr/bin/env python3
"""Prove a finite 217-generator model and exclude an entire line of residual directions."""
import argparse
import itertools
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path
import numpy as np
import sympy as sp
from s61_affine_finite_normal_module_20261102 import slice_model,ROOT,P,compile_native,sha,pivot_minor

STEM='S61_AFFINE_BORDER_SUPPORT_20261102'


def source_entries(quadrics,retained):
    point=np.zeros(46,dtype=np.int64)
    point[0]=1
    point[26:]=retained
    return np.column_stack([((quadrics+quadrics.transpose(0,2,1))@point)[:,1:26],
                            np.einsum('i,rij,j->r',point,quadrics,point)])%P


def reduced_rows(wd,quadrics,retained_inputs,native):
    mon=list(itertools.combinations_with_replacement(range(1,26),2))
    leading=np.column_stack([quadrics[:,i,j]+(quadrics[:,j,i] if i!=j else 0) for i,j in mon])%P
    with (wd/'input').open('w') as f:
        f.write(str(len(retained_inputs))+'\n'+' '.join(map(str,leading.ravel()))+'\n')
        for i,v in enumerate(retained_inputs):
            f.write(str(i)+'\n'+' '.join(map(str,source_entries(quadrics,v).ravel()))+'\n')
    compile_native(native,wd/'border')
    env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
    subprocess.run([wd/'border',wd/'input',wd/'output'],env=env,check=True,timeout=150)
    samples=[];ranks=[]
    with (wd/'output').open()as f:
        tag,rank2,rank3,minor2,minor3=f.readline().split()
        assert tag=='fixed' and list(map(int,[rank2,rank3]))==[134,2925]
        p2=list(map(int,f.readline().split()));p3=list(map(int,f.readline().split()))
        for i in range(len(retained_inputs)):
            label,rank=map(int,f.readline().split());assert label==i
            rows=np.array([np.fromstring(f.readline(),dtype=np.int64,sep=' ') for _ in range(425)])
            assert rows.shape==(425,217)
            samples.append(rows);ranks.append(rank)
    return np.array(samples),ranks,dict(leading_quadratic_minor=int(minor2),leading_quadratic_pivot_columns=p2,
        leading_cubic_minor=int(minor3),leading_cubic_pivot_rows=p3)


def certify():
    start=time.monotonic()
    _,quadrics,_,paths=slice_model()
    direction=np.random.default_rng(1625433).integers(-4,5,size=20)
    native=ROOT/'code/s61_affine_border_reduction_20261102.cpp'
    rowmodule=ROOT/'code/s61_binary_hyperplane_vertical_rowmodule_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-border-support-') as tmp:
        wd=Path(tmp)
        values,ranks,pivots=reduced_rows(wd,quadrics,[i*direction for i in range(6)],native)
        assert ranks==[212,217,217,217,217,217]
        vandermonde=np.array([[pow(i,d,P) for d in range(4)] for i in range(6)],dtype=np.int64)
        inverse=np.array(sp.Matrix(vandermonde[:4].tolist()).inv_mod(P)).astype(np.int64)
        coeffs=(inverse@values[:4].reshape(4,-1)%P).reshape(4,425,217)
        assert np.array_equal((vandermonde@coeffs.reshape(4,-1)%P).reshape(6,425,217),values)
        # Degrees follow from the constant Schur reductions and homogeneous source quadrics.
        limits=[1]*191+[2]*25+[3]
        highest=np.zeros((425,217),dtype=np.int64)
        for i,bound in enumerate(limits):
            assert not np.any(coeffs[bound+1:,:,i])
            highest[:,i]=coeffs[bound,:,i]
        infinity_rows,cols,infinity_det=pivot_minor(highest)
        assert cols==list(range(217)) and infinity_det
        with (wd/'univariate').open('w') as f:
            f.write('32003 425 217 0\n0 0\n')
            for row in coeffs.transpose(1,2,0):
                parts=[]
                for i,poly in enumerate(row):
                    support=np.flatnonzero(poly)
                    if len(support):
                        size=int(support[-1])+1
                        parts.extend([str(i),str(size)]+list(map(str,poly[:size])))
                f.write(str(np.count_nonzero(np.any(row,axis=1)))+' '+' '.join(parts)+'\n')
        compile_native(rowmodule,wd/'rows')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        subprocess.run([wd/'rows',wd/'univariate',wd/'module'],env=env,check=True,
                       timeout=40,capture_output=True)
        module=(wd/'module').read_text().splitlines()
        assert 'ROWMODULE 217 DEGREE 5 ' in module[0]
        det=list(map(int,module[1].split()[2:]))
        assert len(det)==6 and not any(det[:5]) and det[5]
        assert 'FACTOR 5 2 0 1' in module
    paths.extend([Path(__file__).resolve(),native,rowmodule,ROOT/'code/s61_affine_finite_normal_module_20261102.py',
                  ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py',
                  ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',ROOT/'research/notes'/(STEM+'.md')])
    return dict(id=STEM,all_passed=True,prime=P,affine_variable_count=25,retained_base_variables=20,
                reduced_affine_algebra_generator_counts=[1,25,191],finite_module_bound=217,
                polynomial_obstruction_matrix_shape=[425,217],obstruction_column_degrees=limits,
                pivots=pivots,retained_sample_direction=direction.tolist(),sample_ranks=ranks,
                all_lambda_rowmodule_determinant_modp=det,lambda_infinity_217_minor_rows=infinity_rows,
                lambda_infinity_217_minor_mod_prime=infinity_det,projective_section_dimension=26,
                projective_section_reduced_length=5,normal_direction_image_dimension_upper=18,
                whole_residual_support_empty_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    out=parser.parse_args().output
    assert not out.exists();result=certify()
    with out.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(f"PASS 425×217, line determinant {result['all_lambda_rowmodule_determinant_modp']}, {result['seconds']:.3f} seconds")
