#!/usr/bin/env python3
"""Certify a projective 25-plane section by its full cubic ideal."""
import argparse
import json
import os
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement
from pathlib import Path
import numpy as np
import sympy as sp
from s61_affine_normal_projection_20261102 import ROOT, P, compile_native, load_fiber, residue, pivot_minor, sha

STEM='S61_AFFINE_TWENTYFIVE_PLANE_20261102'


def certify():
    started=time.monotonic()
    equations, q_ids, center, original, point_record=load_fiber()
    five_path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    five=json.loads(five_path.read_text())
    assert five['affine_span_dimension']==3 and five['intersection_reduced']
    # Homogeneous coordinates t,x1,x2,x3,n1,...,n22.
    n=26
    coordinates=np.zeros((46,n),dtype=np.int64)
    coordinates[0,0]=1
    x=sp.symbols('v0:3')
    coordinates[1:,1:4]=[[residue(sp.sympify(v).coeff(y)) for y in x] for v in five['affine_plane_substitution']]
    normals=np.random.default_rng(1496271).integers(-4,5,size=(42,22))
    assert len(pivot_minor(normals)[0])==22
    coordinates[4:,4:]=normals % P
    squares=np.zeros((134,46,46),dtype=np.int64)
    for i,poly in enumerate(equations):
        for mon,coeff in poly.items():
            if len(mon)==2:
                squares[i,mon[0]+1,mon[1]+1]=residue(coeff)
            elif len(mon)==1:
                squares[i,0,mon[0]+1]=residue(coeff)
            else:
                assert coeff==0
    Q=coordinates.T @ squares @ coordinates % P
    quadratic_monomials=list(combinations_with_replacement(range(n),2))
    polynomial_rows=np.column_stack([Q[:,i,j]+(Q[:,j,i] if i!=j else 0) for i,j in quadratic_monomials]) % P
    # The five geometrical points split at this prime. Evaluate only pure plane monomials.
    t=sp.Symbol('t');cubic=sp.sympify(five['toric_cubic'])
    conjugate=int(five['toric_quadratic_conjugate_witness']['t_mod_prime'])
    roots=[0,conjugate,-conjugate]
    assert all(residue(cubic.subs(t,r))==0 for r in roots) and len(set(r%P for r in roots))==3
    points=[]
    for r in roots:
        points.append([1]+[residue(sp.sympify(five['toric_S61_coordinates'][q_ids[i]]).subs(t,r)-center[i]) for i in range(3)])
    for name in ('q51_S61_coordinates','q52_S61_coordinates'):
        points.append([1]+[residue(sp.sympify(five[name][q_ids[i]])-center[i]) for i in range(3)])
    assert len({tuple(row) for row in points})==5
    for row in points:
        z=np.array(row+[0]*22,dtype=np.int64)
        assert not np.any(polynomial_rows@np.array([z[i]*z[j]%P for i,j in quadratic_monomials])%P)
    all_cubics=list(combinations_with_replacement(range(n),3))
    basis=list(combinations_with_replacement(range(4),3))
    evaluation=np.array([[int(row[i])*int(row[j])*int(row[k])%P for i,j,k in basis] for row in points],dtype=np.int64)
    rr,cc,evdet=pivot_minor(evaluation)
    assert len(rr)==5
    omitted=[all_cubics.index(basis[c]) for c in cc]
    native=ROOT/'code/s61_affine_twentyfive_plane_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-25-plane-') as tmp:
        wd=Path(tmp)
        (wd/'input').write_text(' '.join(map(str,omitted))+'\n'+' '.join(map(str,polynomial_rows.ravel()))+'\n')
        compile_native(native,wd/'rank')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        subprocess.run([wd/'rank',wd/'input',wd/'out'],check=True,env=env,timeout=90)
        head,detail=(wd/'out').read_text().splitlines()
    rows,cols,rank,det=map(int,head.split())
    selected=list(map(int,detail.split()))
    assert (rows,cols,rank)==(3484,3276,3271) and det and len(selected)==len(set(selected))==rank
    paths=[original,point_record,five_path,Path(__file__).resolve(),native,ROOT/'code/s61_affine_normal_projection_20261102.py',
           ROOT/'code/s61_affine_fiber_macaulay.py',ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',
           ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json',ROOT/'reports/S61_FIVE_COVER_PRODUCT_20261101.json',
           ROOT/'research/notes'/(STEM+'.md')]
    return dict(id=STEM,all_passed=True,sample_prime=P,fixed_normal_plane_integer_basis=normals.tolist(),
                plane_section_projective_dimension=25,known_plane_projective_dimension=3,source_equations=134,
                degree=3,macaulay_shape=[rows,cols],macaulay_exact_rank=rank,
                omitted_cubic_monomial_indices=omitted,omitted_cubic_monomials=[basis[i] for i in cc],
                evaluation_minor_mod_prime=evdet,selected_cubic_minor_rows=selected,selected_cubic_minor_mod_prime=det,
                section_reduced_length=5,coordinate_normal_direction_plane_dimension=21,
                generic_beta_and_Grassmannian_open=True,full_fiber_exhaustion_proved=False,
                input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    out=parser.parse_args().output
    assert not out.exists()
    data=certify()
    with out.open('x') as f:
        json.dump(data,f,indent=2);f.write('\n')
    print(f"PASS: rank {data['macaulay_exact_rank']}, actual minor {data['selected_cubic_minor_mod_prime']}, {data['seconds']:.3f} seconds")
