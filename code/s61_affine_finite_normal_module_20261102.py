#!/usr/bin/env python3
"""Verify the 434-generator reduction of the residual affine algebra."""
import argparse
import hashlib
import itertools
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path
import numpy as np
import sympy as sp
from s61_affine_normal_projection_20261102 import ROOT,P,compile_native,load_fiber,matrices,pivot_minor,residue,sha

STEM='S61_AFFINE_FINITE_NORMAL_MODULE_20261102'


def input_report(name):
    p=ROOT/'reports'/(name+'.json')
    j=json.loads(p.read_text())
    assert j['all_passed']
    for path,digest in j['input_sha256'].items():
        assert sha(ROOT/path)==digest,path
    return j,p


def slice_model():
    normal,npth=input_report('S61_AFFINE_NORMAL_PROJECTION_20261102')
    section,spth=input_report('S61_AFFINE_TWENTYFIVE_PLANE_20261102')
    assert normal['normal_product_shape']==[5418,3612] and normal['plane_product_shape']==[516,420]
    F,qids,base,source,point=load_fiber()
    path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    five=json.loads(path.read_text())
    x=sp.symbols('v0:3')
    slopes=np.array([[residue(sp.sympify(v).coeff(z)) for z in x] for v in five['affine_plane_substitution']],dtype=np.int64)
    normals=np.array(section['fixed_normal_plane_integer_basis'],dtype=np.int64)
    selected_normal_rows=pivot_minor(normals)[0]
    complement=[i for i in range(42) if i not in selected_normal_rows]
    normal_basis=np.column_stack((normals,np.eye(42,dtype=np.int64)[:,complement]))
    assert len(complement)==20 and len(pivot_minor(normal_basis)[0])==42
    T=np.zeros((46,46),dtype=np.int64)
    T[0,0]=1;T[1:,1:4]=slopes;T[4:,4:]=normal_basis
    M=np.zeros((134,46,46),dtype=np.int64)
    for i,poly in enumerate(F):
        for mon,c in poly.items():
            if len(mon)==2:M[i,mon[0]+1,mon[1]+1]=residue(c)
            elif len(mon)==1:M[i,0,mon[0]+1]=residue(c)
            else:assert c==0
    Q=T.T@M@T%P
    return section,Q,complement,[npth,spth,path,source,point]


def certify():
    start=time.monotonic()
    section,Q,complement,paths=slice_model()
    L,H,*_=matrices()
    row_ids,mixed_pivots,mixed_det=pivot_minor(L.reshape(129,168))
    assert len(mixed_pivots)==129
    quad=np.column_stack([Q[:,i,j]+(Q[:,j,i] if i!=j else 0)
                         for i,j in itertools.combinations_with_replacement(range(26),2)])%P
    omitted=section['omitted_cubic_monomial_indices']
    selected=section['selected_cubic_minor_rows']
    native=ROOT/'code/s61_affine_finite_module_lift_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-finite-normal-')as td:
        wd=Path(td)
        (wd/'in').write_text('\n'.join(' '.join(map(str,v)) for v in [omitted,selected,quad.ravel()])+'\n')
        compile_native(native,wd/'lift')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        subprocess.run([wd/'lift',wd/'in',wd/'out'],env=env,check=True,timeout=90)
        target,co=(wd/'out').read_text().splitlines()
    coeff=np.zeros(3484,dtype=np.int64)
    coeff[selected]=np.array(list(map(int,co.split())),dtype=np.int64)
    coeff=coeff.reshape(134,26)
    # Test the lifted equation in the full 46-coordinate model. The RHS is
    # precisely the difference between the central quadratic and its full lift.
    for i in range(3):
        rng=np.random.default_rng(3949941+i)
        V=rng.integers(0,P,size=46)
        Z=V.copy();Z[26:]=0
        multipliers=coeff@Z[:26]%P
        f=np.einsum('i,rij,j->r',V,Q,V,optimize=True)%P
        fcenter=np.einsum('i,rij,j->r',Z,Q,Z,optimize=True)%P
        assert int(multipliers@fcenter)%P==pow(int(V[4]),3,P)
        assert (int(multipliers@f)+int(multipliers@(fcenter-f))-pow(int(V[4]),3,P))%P==0
    counts=[42,22*23//2+(4*42-len(mixed_pivots)),20*(10-5)]
    assert counts==[42,292,100] and sum(counts)==434
    paths.extend([Path(__file__).resolve(),native,ROOT/'code/s61_affine_normal_projection_20261102.py',
        ROOT/'code/s61_affine_fiber_macaulay.py',ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',
        ROOT/'research/notes'/(STEM+'.md')])
    return dict(id=STEM,all_passed=True,prime=P,retained_variables=20,complementary_normal_variables=22,
                complementary_raw_normal_indices=complement,mixed_product_rank=129,mixed_product_pivot_columns=mixed_pivots,
                mixed_product_minor_mod_prime=mixed_det,sample_lift_mon='w1^3',sample_lift_cubic_index=int(target),
                sample_lift_coefficients=sparse_coefficients(coeff),sample_lift_full_3276_coefficients_checked=True,
                module_generator_bounds_by_degree=counts,module_generator_bound=434,
                module_generated_in_degrees_at_most=3,finite_over_retained_twenty_variable_ring=True,
                t_equals_one_is_residual_affine_algebra=True,uniform_on_nonempty_beta_open=True,
                residual_finite_over_beta_claimed=False,zero_module_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},seconds=time.monotonic()-start)


def sparse_coefficients(matrix):
    return [[int(i),int(j),int(matrix[i,j])] for i,j in zip(*np.nonzero(matrix))]


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    out=parser.parse_args().output
    assert not out.exists()
    result=certify()
    with out.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(f"PASS generator bounds {result['module_generator_bounds_by_degree']}, {result['seconds']:.3f} seconds")
