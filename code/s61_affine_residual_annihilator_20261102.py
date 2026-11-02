#!/usr/bin/env python3
"""Construct a nonzero-origin determinant annihilating the entire residual algebra."""
import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path
import numpy as np
import sympy as sp
from itertools import combinations_with_replacement
from s61_affine_border_support_20261102 import ROOT,P,compile_native,sha,pivot_minor,slice_model,source_entries
from s61_affine_normal_projection_20261102 import residue,load_fiber


def certify():
    start=time.monotonic()
    section,Q,complement,paths=slice_model()
    fives=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    j=json.loads(fives.read_text())
    prior=ROOT/'reports/S61_AFFINE_BORDER_SUPPORT_20261102.json'
    border=json.loads(prior.read_text())
    assert border['all_passed']
    for a,h in border['input_sha256'].items():assert sha(ROOT/a)==h
    F,ids,center,*_=load_fiber();T=sp.Symbol('t')
    points=[]
    root=int(j['toric_quadratic_conjugate_witness']['t_mod_prime'])
    for value in (0,root,-root):
        points.append([residue(sp.sympify(j['toric_S61_coordinates'][ids[i]]).subs(T,value)-center[i]) for i in range(3)]+[0]*22)
    for key in ('q51_S61_coordinates','q52_S61_coordinates'):
        points.append([residue(sp.sympify(j[key][ids[i]])-center[i]) for i in range(3)]+[0]*22)
    mon=list(combinations_with_replacement(range(25),2))
    H=np.column_stack([Q[:,i+1,k+1]+(Q[:,k+1,i+1] if i!=k else 0) for i,k in mon])%P
    native=ROOT/'code/s61_affine_origin_differential_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-312-annihilator-')as tmp:
        wd=Path(tmp)
        with (wd/'probes').open('w')as f:
            for point in points:
                V=np.array([1]+point+[0]*20,dtype=np.int64)
                assert not np.any(np.einsum('i,rij,j->r',V,Q,V)%P)
                derivative=((Q+Q.transpose(0,2,1))@V)[:,26:]%P
                f.write(' '.join(map(str,point))+'\n'+' '.join(map(str,derivative.ravel()))+'\n')
        (wd/'input').write_text('1\n'+' '.join(map(str,H.ravel()))+'\n0\n'+' '.join(map(str,source_entries(Q,[0]*20).ravel()))+'\n')
        compile_native(native,wd/'native')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        subprocess.run([wd/'native',wd/'input',wd/'output',wd/'probes'],env=env,check=True,timeout=100)
        with (wd/'output').open()as f:
            tag,h2,h3,m2,m3=f.readline().split()
            assert (tag,int(h2),int(h3),int(m2),int(m3))==('fixed',134,2925,18251,22132)
            p2=list(map(int,f.readline().split()))
            assert p2==border['pivots']['leading_quadratic_pivot_columns']
            assert list(map(int,f.readline().split()))==border['pivots']['leading_cubic_pivot_rows']
            d_eval=np.array([[np.fromstring(f.readline(),dtype=np.int64,sep=' ') for _ in range(425)] for _ in points])
            assert f.readline().strip()=='0 212'
            M0=np.array([np.fromstring(f.readline(),dtype=np.int64,sep=' ') for i in range(425)])
    pivot_rows,pivot_cols,pivot0=pivot_minor(M0)
    assert len(pivot_cols)==212
    free=[i for i in range(217) if i not in pivot_cols]
    X2free=[mon[i] for i in range(325) if i not in p2]
    evaluations=np.array([[z[i]*z[k]%P for i,k in X2free]+z+[1] for z in points],dtype=np.int64).T
    assert not np.any(M0@evaluations%P)
    transform=np.array(sp.Matrix(evaluations[free].tolist()).inv_mod(P)).astype(np.int64)
    normalized=evaluations@transform%P
    assert np.array_equal(normalized[free],np.eye(5,dtype=np.int64))
    # d_eval[s,r,i] is the row derivative in y_i evaluated at known point s.
    deriv=(d_eval.transpose(1,2,0)@transform)%P
    augmented=np.column_stack([M0[:,pivot_cols],deriv.reshape(425,-1)])
    selected,cols,det=pivot_minor(augmented)
    assert len(selected)==312 and cols==list(range(312)) and det
    # The weaker row degrees below already bound the determinant in all y.
    degrees=border['obstruction_column_degrees']
    upper=sum(degrees[q] for q in pivot_cols)+100*2
    st='S61_AFFINE_RESIDUAL_ANNIHILATOR_20261102'
    paths+=[fives,prior,native,Path(__file__).resolve(),ROOT/'code/s61_affine_finite_normal_module_20261102.py',
            ROOT/'code/s61_affine_border_support_20261102.py',ROOT/'code/s61_affine_normal_projection_20261102.py',
            ROOT/'research/notes'/(st+'.md')]
    return dict(id=st,all_passed=True,prime=P,retained_variables=20,matrix_at_origin_rank=212,
                constant_pivot_rows=pivot_rows,constant_pivot_columns=pivot_cols,constant_pivot_minor_modp=pivot0,
                constant_free_columns=free,augmented_shape=[425,312],linear_residual_directions=100,
                determinant_rows=selected,annihilator_at_origin_modp=det,degree_upper=upper,
                determinant_definition='det([M(y)[:,pivot columns], H_1(y),..., H_20(y)][selected rows,:]); M(y)K=sum y_i H_i(y), with polynomial telescoping along the ordered y coordinates, K normalized on five free columns.',
                annihilates_complete_residual_affine_algebra=True,known_five_restriction_kernel=True,
                nonzero_at_origin=True,entire_residual_empty_claimed=False,input_sha256={str(a.relative_to(ROOT)):sha(a) for a in paths},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output;assert not out.exists()
    cert=certify()
    with out.open('x')as f:json.dump(cert,f,indent=2);f.write('\n')
    print(f"PASS D(0)={cert['annihilator_at_origin_modp']} mod 32003, degree <= {cert['degree_upper']}, {cert['seconds']:.3f} seconds")
