#!/usr/bin/env python3
"""Certify projection from the known three-plane by two polarized Macaulay ranks."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
import sympy as sp

from s61_affine_fiber_macaulay import load_fiber
from s61_binary_hyperplane_cubic_matrix_20261102 import compile_native

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_AFFINE_NORMAL_PROJECTION_20261102'
P = 32003
NORMALS = 42
PLANES = 4


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residue(value):
    value = Fraction(value)
    assert value.denominator % P
    return value.numerator * pow(value.denominator, -1, P) % P


def pivot_minor(matrix):
    """Return the row and column choices and the signed determinant of that minor."""
    reduced = matrix.copy() % P
    order = list(range(len(reduced)))
    rows, columns = [], []
    pivots_product = 1
    for j in range(matrix.shape[1]):
        r = len(rows)
        if r == len(reduced):
            break
        choices = np.flatnonzero(reduced[r:,j])
        if not len(choices):
            continue
        swap = r + int(choices[0])
        reduced[[r,swap]] = reduced[[swap,r]]
        order[r],order[swap] = order[swap],order[r]
        rows.append(order[r])
        columns.append(j)
        pivot = int(reduced[r,j])
        pivots_product = pivots_product * pivot % P
        reduced[r] = reduced[r] * pow(pivot,-1,P) % P
        reduced[r+1:,j:] = (reduced[r+1:,j:] - reduced[r+1:,j,None] * reduced[r,j:]) % P
    return rows, columns, pivots_product


def matrices():
    F, q_columns, center, original, point = load_fiber()
    known_path = ROOT / 'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    known = json.loads(known_path.read_text())
    assert len(F) == 134 and known['q_to_u_indices_zero_based'] == q_columns
    assert known['affine_plane_free_q_indices'] == [0,1,2]
    plane_vars = sp.symbols('v0:3')
    expressions = [sp.sympify(value) for value in known['affine_plane_substitution']]
    # Translate the source rational affine chart to the rational toric point.
    # Coordinates z1,z2,z3 measure the displacements of the first three q's.
    center_values = dict(zip(plane_vars,center[:3],strict=True))
    assert all(expr.subs(center_values) == val for expr,val in zip(expressions,center,strict=True))
    embedding = np.array([[residue(expr.coeff(x)) for x in plane_vars] for expr in expressions],dtype=np.int64)
    assert np.array_equal(embedding[:3],np.eye(3,dtype=np.int64))
    linear = np.zeros((134,45),dtype=np.int64)
    quadratic = np.zeros((134,45,45),dtype=np.int64)
    for i,eq in enumerate(F):
        for mon,coef in eq.items():
            if len(mon) == 1:
                linear[i,mon[0]] = residue(coef)
            elif len(mon) == 2:
                quadratic[i,mon[0],mon[1]] = residue(coef)
            else:
                assert coef == 0
    mixed = (quadratic + quadratic.transpose(0,2,1)) % P
    pure = embedding.T @ quadratic @ embedding % P
    restriction = np.column_stack([linear @ embedding % P]+[
        pure[:,i,j]+(pure[:,j,i] if i!=j else 0) for i,j in combinations_with_replacement(range(3),2)]) % P
    rows, cols, det = pivot_minor(restriction)
    assert rows == [0,1,2,3,4] and cols == [0,1,2,3,4]
    inverse = np.array(sp.Matrix(restriction[:5,cols].tolist()).inv_mod(P)).astype(np.int64)
    change = restriction[:,cols] @ inverse % P
    assert np.array_equal(restriction,change @ restriction[:5] % P)
    assert np.array_equal(change,np.array([[residue(c) for c in r] for r in known['all134_restricted_span_coefficients']]))
    # z0 is the homogenizing t. The chosen normal complement is q indices 3..44.
    cross = np.concatenate([linear[:,:,None], mixed @ embedding % P],axis=2)[:,3:,:]
    residual_linear = (cross[5:]-(change[5:] @ cross[:5].reshape(5,-1)).reshape(129,42,4)) % P
    nq = quadratic[:,3:,3:]
    residual_quadratic = (nq[5:]-(change[5:] @ nq[:5].reshape(5,-1)).reshape(129,42,42)) % P
    return residual_linear,residual_quadratic,det,[original,point,known_path]


def certify():
    started = time.monotonic()
    cross, quad, plane_det, paths = matrices()
    pairs = list(combinations_with_replacement(range(PLANES),2))
    plane_product = np.zeros((129*PLANES, NORMALS*len(pairs)),dtype=np.int64)
    for j in range(PLANES):
        for k in range(PLANES):
            col = pairs.index(tuple(sorted([j,k]))) * NORMALS
            plane_product[j*129:(j+1)*129,col:col+NORMALS] = cross[:,:,k]
    small_rows,small_cols,small_det = pivot_minor(plane_product)
    assert len(small_rows) == 420 and small_cols == list(range(420))
    # A 15-plane through the known three-plane can already be exhausted at degree two.
    normal_slice=np.random.default_rng(1496271).integers(-4,5,size=(42,12))
    mixed_slice=normal_slice.T @ cross % P
    quadratic_slice=normal_slice.T @ quad @ normal_slice % P
    section_rows=np.column_stack([mixed_slice.reshape(129,-1)]+[
        quadratic_slice[:,i,j]+(quadratic_slice[:,j,i] if i!=j else 0)
        for i,j in combinations_with_replacement(range(12),2)]) % P
    section_pivots,section_cols,section_det=pivot_minor(section_rows)
    assert len(section_pivots)==126 and section_cols==list(range(126)) and section_det==2962
    native = ROOT / 'code/s61_affine_normal_projection_20261102.cpp'
    with tempfile.TemporaryDirectory(prefix='s61-affine-projection-') as tmp:
        wd=Path(tmp)
        (wd/'input').write_text(' '.join(map(str,cross.ravel()))+'\n')
        compile_native(native,wd/'rank')
        env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
        subprocess.run([wd/'rank',wd/'input',wd/'output'],check=True,env=env,timeout=120)
        raw=(wd/'output').read_text().splitlines()
    nr,nc,det=map(int,raw[0].split())
    assert (nr,nc)==(5418,3612) and det
    large_rows=list(map(int,raw[1].split()))
    assert len(large_rows)==len(set(large_rows))==3612
    cof=np.zeros((129*NORMALS,PLANES),dtype=np.int64)
    cof[large_rows]=np.array([list(map(int,row.split())) for row in raw[2:]],dtype=np.int64)
    cof=cof.reshape(NORMALS,129,PLANES)
    # Independently reconstruct normal cubic numerators on the chart n0=1.
    # The relation is z_j n0^2 + sum c_{li} n_l Q_i(n) = 0.
    numerators=-(cof.transpose(0,2,1)@quad.reshape(129,-1)).reshape(42,4,42,42) % P
    triples=list(combinations_with_replacement(range(NORMALS),3))
    coeff=np.zeros((PLANES,len(triples)),dtype=np.int64)
    from itertools import permutations
    for m,mon in enumerate(triples):
        for i,j,k in set(permutations(mon)):
            coeff[:,m]+=numerators[i,:,j,k]
    coeff%=P
    # All 4 coordinates of arbitrary (z,n) check the specialized cubic identity.
    for l in (1,2,3):
        rng=np.random.default_rng(6112026+l)
        z=rng.integers(0,P,size=4);n=rng.integers(0,P,size=42)
        vals=(np.einsum('i,rij,j->r',n,cross,z,optimize=True)+np.einsum('i,rij,j->r',n,quad,n,optimize=True))%P
        L=(n@cof.reshape(42,-1)%P).reshape(129,4)
        rhs=coeff@np.array([int(n[i])*int(n[j])*int(n[k])%P for i,j,k in triples])%P
        assert np.array_equal((z*(int(n[0])**2)-rhs)%P,vals@L%P)
    paths+=[native,Path(__file__).resolve(),ROOT/'code/s61_affine_fiber_macaulay.py',
            ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',
            ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json',
            ROOT/'reports/S61_FIVE_COVER_PRODUCT_20261101.json',
            ROOT/'research/notes'/(STEM+'.md')]
    return dict(id=STEM,all_passed=True,prime=P,equations_used=134,independent_plane_quadrics=5,
                center='rational toric point for coefficients (2,3,5,7,11,13)',restriction_pivot_minor=plane_det,
                normal_equations=129,plane_product_shape=[516,420],plane_product_minor_rows=small_rows,
                plane_product_minor_mod_prime=small_det,normal_product_shape=[nr,nc],normal_product_minor_rows=large_rows,
                normal_product_minor_mod_prime=det,chart_n0_modp_cubic_terms=[int(np.count_nonzero(row)) for row in coeff],
                chart_n0_modp_cubic_digest=hashlib.sha256(coeff.astype('<u2').tobytes()).hexdigest(),
                exact_QQ_generic_cubic_reconstruction=True,projective_normal_coordinates=42,
                normalized_direction_chart_count=42,highest_projected_generator_degree=4,
                fifteen_plane_normal_basis=normal_slice.tolist(),fifteen_plane_sample_minor_rows=section_pivots,
                fifteen_plane_sample_minor=section_det,fifteen_plane_projective_scheme_length=5,
                all_generic_affine_points_exhausted=False,input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output
    assert not out.exists()
    cert=certify()
    with out.open('x') as f:
        json.dump(cert,f,indent=2);f.write('\n')
    print(f"PASS ranks 420 and 3612 in {cert['seconds']:.3f} seconds")
