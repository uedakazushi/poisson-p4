#!/usr/bin/env python3
"""Finite polynomial reduction over twenty normals after prescribing u6 and u42."""
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

STEM='S61_AFFINE_FIXED_PAIR_LOCALIZATION_20261102'


def certify():
    start=time.monotonic()
    section,Q,standard_axes,paths=slice_model()
    basis=np.column_stack([np.array(section['fixed_normal_plane_integer_basis']),np.eye(42,dtype=np.int64)[:,standard_axes]])
    five_path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    projection=ROOT/'reports/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.json'
    five=json.loads(five_path.read_text())
    _,ids,base,_,_=load_fiber()
    names=sp.symbols('v0:3')
    idx=[ids.index(i) for i in [5,41]]
    slopes=sp.Matrix([[sp.sympify(five['affine_plane_substitution'][i]).coeff(z) for z in names] for i in idx])
    pivot=sp.factor(slopes[:,:2].det());assert pivot==-sp.Rational(237707685,145350016) and residue(pivot)
    T=np.array([[residue(base[i])]+[residue(v) for v in slopes[j,:]]+list(basis[i-3]) for j,i in enumerate(idx)],dtype=np.int64)%P
    free=[0]+list(range(3,46));inverse=np.array(sp.Matrix(T[:,1:3].tolist()).inv_mod(P)).astype(np.int64)
    substitution=np.zeros((46,44),dtype=np.int64);substitution[free]=np.eye(44,dtype=np.int64)
    substitution[1:3]=-inverse@T[:,free]%P;assert not np.any(T@substitution%P)
    transformed=substitution.T@Q@substitution%P
    mon=list(itertools.combinations_with_replacement(range(1,24),2))
    quadratic=np.column_stack([transformed[:,a,b]+(transformed[:,b,a] if a!=b else 0) for a,b in mon])%P
    directions=np.random.default_rng(993217).integers(-5,6,size=(2,20))
    terms=[(a,b) for a in range(4) for b in range(4-a)]
    samples=terms+[(5,2),(3,3)]
    native=ROOT/'code/s61_affine_fixed_pair_border_20261102.cpp'
    mac=ROOT/'code/s61_affine_fixed_pair_parameter_matrix_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-166-')as tmp:
        wd=Path(tmp)
        with (wd/'input').open('w')as f:
            f.write(str(len(samples))+'\n'+' '.join(map(str,quadratic.ravel()))+'\n')
            for i,(a,b) in enumerate(samples):
                z=np.array([1]+[0]*23+list(a*directions[0]+b*directions[1]),dtype=np.int64)
                linear=((transformed+transformed.transpose(0,2,1))@z)[:,1:24]
                constant=np.einsum('i,rij,j->r',z,transformed,z,optimize=True)
                lower=np.column_stack([linear,constant])%P
                f.write(str(i)+'\n'+' '.join(map(str,lower.ravel()))+'\n')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        compile_native(native,wd/'border')
        subprocess.run([wd/'border',wd/'input',wd/'out'],check=True,env=env,timeout=100)
        matrices=[]
        with (wd/'out').open()as f:
            header=f.readline().split()
            assert header==['fixed','134','2300','15981','3631']
            pivot2=list(map(int,f.readline().split()));pivot3=list(map(int,f.readline().split()))
            for i in range(len(samples)):
                assert list(map(int,f.readline().split()))==[i,166]
                matrices.append(np.array([np.fromstring(f.readline(),dtype=np.int64,sep=' ') for _ in range(782)]))
        matrices=np.array(matrices)
        V=np.array([[pow(a,x,P)*pow(b,y,P)%P for x,y in terms] for a,b in samples],dtype=np.int64)
        inverse=np.array(sp.Matrix(V[:10].tolist()).inv_mod(P)).astype(np.int64)
        coeff=(inverse@matrices[:10].reshape(10,-1)%P).reshape(10,782,166)
        assert np.array_equal((V@coeff.reshape(10,-1)%P).reshape(matrices.shape),matrices)
        with (wd/'poly').open('w')as f:
            for index,(i,j) in enumerate(terms):
                f.write(f'{i} {j} '+' '.join(map(str,coeff[index].ravel()))+'\n')
        compile_native(mac,wd/'mac')
        subprocess.run([wd/'mac',wd/'poly',wd/'macout','4'],check=True,env=env,timeout=45)
        det_line,selected,summary=(wd/'macout').read_text().splitlines()
        assert list(map(int,summary.split()))==[4,2346,1097,1097,0]
        det=int(det_line);rows=list(map(int,selected.split()));assert det and len(set(rows))==1097
    paths.extend([five_path,projection,ROOT/'reports/S61_AFFINE_BORDER_SUPPORT_20261102.json',Path(__file__).resolve(),
                  native,mac,ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_finite_normal_module_20261102.py',
                  ROOT/'code/s61_affine_fiber_macaulay.py',ROOT/'research/notes'/(STEM+'.md')])
    return dict(id=STEM,all_passed=True,prime=P,fixed_coordinates=['u6','u42'],fixed_values_in_check=[0,0],affine_frame_minor_QQ=str(pivot),
                unknown_variables=23,retained_base_variables=20,leading_quadratic_rows=134,leading_quadratic_columns=276,
                leading_quadratic_pivot_columns=pivot2,leading_quadratic_minor_mod_prime=15981,
                leading_cubic_rows=3082,leading_cubic_columns=2300,leading_cubic_pivot_rows=pivot3,
                leading_cubic_minor_mod_prime=3631,generator_counts=[1,23,142],generator_bound=166,
                relation_matrix_shape=[782,166],retained_directions=directions.tolist(),parameter_degree=4,
                chosen_direction_plane_empty_projective_solution=True,parameter_relation_shape=[2346,1097],
                parameter_minor_mod_prime=det,parameter_minor_rows=rows,all_twenty_parameters_empty_proved=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},seconds=time.monotonic()-start)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;assert not out.exists()
    value=certify()
    with out.open('x')as f:json.dump(value,f,indent=2);f.write('\n')
    print(f"PASS: 166 generators and degree-four minor {value['parameter_minor_mod_prime']}; {value['seconds']:.3f} seconds")
