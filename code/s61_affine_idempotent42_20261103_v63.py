#!/usr/bin/env python3
"""Normalize one row of the linear normal matrix: a global degree42 idempotent."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from s61_affine_normal_projection_20261102 import matrices
from s61_affine_polarization_splitter_20261103_v63 import determinant_including_zero

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    packet=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63';old=json.loads((packet/'report.json').read_text())
    assert old['all_passed']
    for path,value in old['input_sha256'].items():assert sha(ROOT/path)==value
    for path,value in old['evidence_sha256'].items():assert sha(packet/path)==value
    data=ROOT/'code/data/s61_affine_idempotent42_20261103_v63.json';spec=json.loads(data.read_text())
    selpath=ROOT/'code/data/s61_affine_residual42_20261103_v63.json';selected=json.loads(selpath.read_text())['selected_residual_rows_zero_based']
    slot=spec['replaced_selected_row_slot'];columns=spec['replacement_residual_row_indices'];assert len(columns)==5
    cross,quad,*_=matrices()
    knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(knownpath.read_text())['known_five_mod_32003'];matrices_at_points=[];cofactor_evaluations=[]
    for point in known:
        A=np.einsum('rnk,k->rn',cross,np.array([1,*point[1:]],dtype=np.int64))%P
        square=A[selected,:];d=determinant_including_zero(square);assert d
        inverse=np.array(sp.Matrix(square.tolist()).inv_mod(P)).astype(np.int64)
        cofactor_evaluations.append((d*(A@inverse[:,slot]%P)%P)[columns].tolist())
        matrices_at_points.append(A)
    C=np.array(cofactor_evaluations,dtype=np.int64);minor=determinant_including_zero(C)
    assert C.tolist()==spec['knownfive_replacement_determinant_matrix']and minor==spec['normalization_minor_mod_prime']==91
    weights=np.array(sp.Matrix(C.tolist()).inv_mod(P)).astype(np.int64)@np.ones(5,dtype=np.int64)%P
    assert weights.tolist()==spec['replacement_weights_mod_prime']
    values=[]
    for A in matrices_at_points:
        square=A[selected,:].copy();square[slot]=weights@A[columns,:]%P
        values.append(determinant_including_zero(square))
    assert values==spec['normalized_known_five_determinants']==[1]*5
    rng=np.random.default_rng(20261103);z=rng.integers(0,P,size=4);z[0]=0;n=rng.integers(0,P,size=42)
    A=(np.einsum('rnk,k->rn',cross,z)+np.einsum('rjk,k->rj',quad,n))%P
    square=A[selected,:].copy();square[slot]=weights@A[columns,:]%P
    leading=determinant_including_zero(square)
    assert leading==spec['normalized_leading_degree42_determinant']==17048
    (args.output/'idempotent42-certificate.json').write_text(json.dumps(dict(
        selected_normal_rows=selected,replaced_slot=slot,replacement_rows=columns,
        coefficient_matrix_at_knownfive=C.tolist(),coefficient_matrix_determinant_mod_prime=minor,
        replacement_weights_at_frozen_prime=weights.tolist(),known_five_determinants=values,
        leading_degree42_determinant_mod_prime=leading,
        generic_beta_recipe='For each replacement row letc_r∈E be its determinant with41fixedrows; solve sum(w_r*c_r)=1_E in the actual5dimensional known algebra.',
        weights_are_beta_rational_not_frozen_integers=True,
        generic_weights_expanded=False,only_beta_5x5_basis_minor_inverted=True,
        polynomial_E42='det(linear42matrix with normalizedreplacementrow)',
        exact_full_algebra_identities=['E42^2=E42','E42=e','E42*n_j=0','I_res=I+(E42)'],
        exact_q_degree=42,no_q_or_cyclic_or_H_localization=True,
        unknown_residual_emptiness_claimed=False),indent=2)+'\n')
    inputs=[Path(__file__),data,selpath,packet/'report.json',knownpath,
        ROOT/'code/s61_affine_normal_projection_20261102.py',
        ROOT/'research/notes/S61_AFFINE_IDEMPOTENT42_20261103_v63.md']
    report=dict(id='S61_AFFINE_IDEMPOTENT42_20261103_v63',all_passed=True,
        normalized_degree42_determinant_is_global_idempotent=True,
        full_unknown_residual_exactly_I_plus_E42=True,
        generic_beta_weights_finite_5x5_recipe_not_expanded=True,
        knownfive_allnormalized_to_one=True,unit_or_full_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS normalizedlinear42det isglobalidempotent;unknownresidualretained;',report['seconds'])


if __name__=='__main__':main()
