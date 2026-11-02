#!/usr/bin/env python3
"""Recover the full finite algebra at a normal direction and separate its nonzero-t part."""
import argparse
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np
import sympy as sp

from s61_affine_finite_normal_module_20261102 import ROOT, P, compile_native, load_fiber, residue, sha, slice_model
from s61_model import S61, X, J

STEM = 'S61_AFFINE_HOMOGENIZER_FINITE_ACTION_20261102'


def certify():
    started=time.monotonic()
    section,quadrics,axes,inputs=slice_model()
    raw_polys,ids,base,_,_=load_fiber()
    normal_basis=sp.Matrix(section['fixed_normal_plane_integer_basis']).row_join(sp.eye(42)[:,axes])
    assert residue(normal_basis.det())!=0
    known_path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    known=json.loads(known_path.read_text())
    variables=sp.symbols('v0:3')
    slopes=sp.Matrix([[sp.sympify(expression).coeff(z) for z in variables]
                      for expression in known['affine_plane_substitution']])
    coord_rows=[ids.index(i) for i in [5,41]]
    coord=sp.Matrix([[base[i],*slopes[i,:],*normal_basis[i-3,:]] for i in coord_rows])
    assert coord[:,1:3].det()==-sp.Rational(237707685,145350016)
    full_inverse=coord[:,1:3].inv()
    environment=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
    border=ROOT/'reports/S61_AFFINE_BORDER_SUPPORT_20261102.json'
    shift=np.array(json.loads(border.read_text())['retained_sample_direction'],dtype=np.int64)

    # Compute the genuine flag as a rational source point before reducing it.
    model=S61()
    form=X[0]**2*J
    form[1,3]+=X[4]**2
    form[3,1]-=X[4]**2
    normalized=model.normalize(form)['parameters']
    q=sp.Matrix([normalized[i] for i in ids])
    assert q==sp.eye(45)[:,8] and all(q[i]==0 for i in coord_rows)
    assert all(sum(c*q[mon[0]]*q[mon[1]] for mon,c in poly.items() if len(mon)==2)==0 for poly in raw_polys)
    offsets=q-slopes*q[:3,:]
    assert offsets[:3,:]==sp.zeros(3,1)
    flag_full=sp.Matrix([0,*q[:3,:],*(normal_basis.inv()*offsets[3:,:])])
    reduce_matrix=lambda x:np.array(x.applyfunc(residue).tolist(),dtype=np.int64)
    flag_mod=reduce_matrix(flag_full).ravel()
    flag_y=flag_mod[26:]
    free=[0]+list(range(3,46))
    parameters=[sp.zeros(2,1),sp.Matrix([base[i] for i in coord_rows])]
    native=ROOT/'code/s61_affine_finite_direction_action_20261102.cpp'
    pairs=[]
    with tempfile.TemporaryDirectory(prefix='s61-191-action-') as temporary:
        wd=Path(temporary)
        compile_native(native,wd/'action')
        for case,pair in enumerate(parameters):
            # Substitute the two fixed coordinates, then the translation y=Y-vt.
            specified=coord.copy()
            specified[:,0]-=pair
            substitution=sp.zeros(46,44)
            for j,row in enumerate(free):
                substitution[row,j]=1
            substitution[1:3,:]=-full_inverse*specified[:,free]
            translation=sp.eye(44)
            translation[24:,0]=-sp.Matrix(shift.tolist())
            substitution*=translation
            assert specified*substitution==sp.zeros(2,44)
            sub=reduce_matrix(substitution)
            raw=sub.T@quadrics@sub%P
            mon=list(itertools.combinations_with_replacement(range(24),2))
            highest=np.column_stack([raw[:,i,j]+(raw[:,j,i] if i!=j else 0) for i,j in mon])%P
            values=[np.zeros(20,dtype=np.int64),shift,flag_y]
            with (wd/'in').open('w') as stream:
                stream.write(f'{len(values)}\n'+' '.join(map(str,highest.ravel()))+'\n')
                for i,y in enumerate(values):
                    zero=np.array([0]*24+list(y),dtype=np.int64)
                    linear=((raw+raw.transpose(0,2,1))@zero)[:,:24]%P
                    constant=np.einsum('i,rij,j->r',zero,raw,zero,optimize=True)%P
                    if i==2 or (case,i)==(1,1):
                        x=flag_mod[free][:24] if i==2 else np.array([1]+[0]*23,dtype=np.int64)
                        original=flag_mod if i==2 else np.array([1]+[0]*45,dtype=np.int64)
                        actual=sub@np.concatenate([x,y])%P
                        assert np.array_equal(actual,original)
                        assert not np.any((highest@np.array([x[k]*x[l]%P for k,l in mon])+linear@x+constant)%P)
                    block=np.column_stack([linear,constant])
                    stream.write(f'{i}\n'+' '.join(map(str,block.ravel()))+'\n')
            subprocess.run([wd/'action',wd/'in',wd/'out'],check=True,env=environment,timeout=150)
            lines=iter((wd/'out').read_text().splitlines())
            tag,d2,d3=next(lines).split()
            assert tag=='pivots' and [int(d2),int(d3)]==[[8028,17223],[6685,2189]][case]
            cols=list(map(int,next(lines).split()));rows=list(map(int,next(lines).split()))
            assert len(set(cols))==134 and len(set(rows))==2600
            checks=[]
            for i in range(3):
                tokens=next(lines).split()
                assert tokens[:3]==['exact',str(i),'trace'] and tokens[-2]=='dimension'
                trace=list(map(int,tokens[3:-2]))
                j,t,tag,s,dtag,length=next(lines).split()
                assert [j,tag,dtag]==['t','stable','dim']
                length,t,s=map(int,[length,t,s])
                expected=([(191,25,0),(0,0,0),(1,0,0)] if case==0 else [(191,25,0),(1,1,1),(1,0,0)])
                assert (length,t,s)==expected[i] and trace[-1]==length
                checks.append(dict(direction=['zero','v','flag'][i],commutator_and_stable_kernel_dimensions=trace,
                                   total_length=length,homogenizer_rank=t,nonzero_homogenizer_length=s))
            assert next(lines)=='toy trace 1 0' and list(lines)==[]
            pairs.append(dict(value=list(map(str,pair)),quadratic_minor=int(d2),cubic_minor=int(d3),
                              quadratic_pivots=cols,cubic_pivots=rows,checks=checks))
    inputs.extend([known_path,border,native,Path(__file__).resolve(),ROOT/'research/notes'/(STEM+'.md'),
                   ROOT/'reports/S61_AFFINE_FIXED_PAIR_LOCALIZATION_20261102.json',
                   ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',
                   ROOT/'code/s61_affine_finite_normal_module_20261102.py',
                   ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py',
                   ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py'])
    return dict(id=STEM,all_passed=True,prime=P,retained_direction_dimension=19,
                normal_shift=list(map(int,shift)),exact_flag_direction=list(map(str,flag_full[26:,:])),
                leading_matrix_shapes=[[134,300],[3216,2600]],generator_counts=[1,24,166],
                homogenizer_candidate_size=191,selected_pair_checks=pairs,
                arbitrary_field_fiber_and_t_action_exact=True,all_directions_excluded_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    path=parser.parse_args().output
    assert not path.exists()
    report=certify()
    with path.open('x')as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print(f"PASS exact finite direction algebras and t-actions; {report['seconds']:.3f}s")
