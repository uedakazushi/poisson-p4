#!/usr/bin/env python3
"""Certify a cyclic univariate presentation of the full affine algebra over an open chart."""
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

from s61_affine_finite_normal_module_20261102 import ROOT, P, compile_native, load_fiber, residue, sha, slice_model

STEM = 'S61_AFFINE_FIXED_PAIR_UNIVARIATE_20261102'


def certify():
    start = time.monotonic()
    section, quadrics, axes, inputs = slice_model()
    normal_basis = np.column_stack([section['fixed_normal_plane_integer_basis'],np.eye(42,dtype=np.int64)[:,axes]])
    known_path = ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    known = json.loads(known_path.read_text())
    _, variable_ids, base, _, _ = load_fiber()
    idx = [variable_ids.index(i) for i in [5,41]]
    plane = sp.symbols('v0:3')
    slopes = sp.Matrix([[sp.sympify(known['affine_plane_substitution'][i]).coeff(z) for z in plane] for i in idx])
    assert slopes[:,:2].det()==-sp.Rational(237707685,145350016)
    coord = np.array([[residue(base[i])]+[residue(v) for v in slopes[j,:]]+list(normal_basis[i-3])
                      for j,i in enumerate(idx)],dtype=np.int64)%P
    inverse = np.array(sp.Matrix(coord[:,1:3].tolist()).inv_mod(P)).astype(np.int64)
    free = [0]+list(range(3,46))
    substitution = np.zeros((46,44),dtype=np.int64)
    substitution[free] = np.eye(44,dtype=np.int64)
    def fix_pair(pair):
        specified = coord.copy()
        specified[:,0] = (specified[:,0]-pair)%P
        substitution[1:3] = -inverse@specified[:,free]%P
        assert not np.any(specified@substitution%P)
        return substitution.T@quadrics@substitution%P

    point_pair = [residue(base[i]) for i in idx]
    directions = np.random.default_rng(63971831)
    linear_form = directions.integers(1,600,size=23)
    retained = list(map(int,directions.integers(-10,11,size=20)))
    samples = [(1,[0,0],[0]*20),(1,point_pair,[0]*20),(1,point_pair,retained),(0,[0,0],[0]*20)]
    monomials = list(itertools.combinations_with_replacement(range(1,24),2))
    def highest(matrix):
        return np.column_stack([matrix[:,i,j]+(matrix[:,j,i] if i!=j else 0) for i,j in monomials])%P

    leading = highest(fix_pair([0,0]))
    native = ROOT/'code/s61_affine_krylov_field_fiber_20261102.cpp'
    results=[]
    with tempfile.TemporaryDirectory(prefix='s61-cyclic-affine-') as tmp:
        wd=Path(tmp)
        with (wd/'in').open('w') as f:
            f.write(f'{len(samples)}\n'+' '.join(map(str,leading.ravel()))+'\n'+' '.join(map(str,linear_form))+'\n')
            for j,(t,pair,y) in enumerate(samples):
                matrix=fix_pair(pair)
                assert np.array_equal(leading,highest(matrix))
                origin=np.array([t]+[0]*23+y,dtype=np.int64)
                linear=((matrix+matrix.transpose(0,2,1))@origin)[:,1:24]
                constant=np.einsum('i,rij,j->r',origin,matrix,origin,optimize=True)
                if j==1:
                    assert not np.any(constant%P)
                low=np.column_stack([linear,constant])%P
                f.write(f'{j}\n'+' '.join(map(str,low.ravel()))+'\n')
        compile_native(native,wd/'cyclic')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        subprocess.run([wd/'cyclic',wd/'in',wd/'out'],check=True,env=env,timeout=120)
        lines=iter((wd/'out').read_text().splitlines())
        assert next(lines)=='pivots 15981 3631'
        for j,(t,pair,y) in enumerate(samples):
            a,b,label,rank,dt,det=next(lines).split()
            assert int(a)==int(b)==j and label=='cyclic' and dt=='det'
            assert (int(rank),int(det))==[(166,30375),(166,1827),(166,3231),(3,0)][j]
            item=dict(t=t,pair_mod_prime=pair,retained_y=y,cyclic_rank=int(rank),cyclic_minor=int(det))
            if int(rank)==166:
                gcd_line,step_line=next(lines).split('steps')
                entry=gcd_line.split()
                assert entry[0]=='univariate' and entry[2]=='gcd'
                gcd=list(map(int,entry[3:]))
                assert gcd==[[1],[0,1],[1]][j] and list(map(int,step_line.split()))==[0]
                characteristic=list(map(int,next(lines).split()[1:]))
                assert len(characteristic)==167 and characteristic[-1]==1
                item.update(monic_gcd_mod_prime=gcd,characteristic_mod_prime=characteristic,
                            first_original_row_already_gives_final_gcd=True)
            results.append(item)
        assert list(lines)==[]
    dependency=ROOT/'reports/S61_AFFINE_FIXED_PAIR_LOCALIZATION_20261102.json'
    report=json.loads(dependency.read_text())
    assert report['all_passed'] and report['generator_bound']==166
    inputs.extend([dependency,known_path,native,Path(__file__).resolve(),ROOT/'research/notes'/(STEM+'.md'),
                   ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
                   ROOT/'code/s61_affine_finite_normal_module_20261102.py',
                   ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py',
                   ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py'])
    return dict(id=STEM,all_passed=True,prime=P,affine_original_equations=134,retained_parameters=20,
                projected_coordinates=['u6','u42'],toric_sample_exact_coordinate_values=[str(base[i]) for i in idx],
                eliminated_variables=23,univariate_degree_upper=166,generator_counts=[1,23,142],
                linear_form_coefficients=list(map(int,linear_form)),sample_checks=results,
                actual_affine_algebra_reconstructed_on_cyclic_open=True,full_fiber_empty_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    output=parser.parse_args().output
    assert not output.exists()
    result=certify()
    with output.open('x')as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(f"PASS three cyclic algebras and the homogeneous noncyclic control; {result['seconds']:.3f}s")
