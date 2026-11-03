#!/usr/bin/env python3
"""Exact five determinant residual equations and a degree43 scalar separator.

All generic coefficients are finite five-dimensional algebra recipes.  The
finite-field calculations certify beta-unit opens and nonzero leading terms;
they are not an emptiness argument.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time

import numpy as np
import sympy as sp

from s61_affine_fiber_macaulay import load_fiber
from s61_affine_polarization_splitter_20261103_v63 import residue,determinant_including_zero

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def rref(matrix):
    matrix=np.array(matrix,dtype=np.int64).copy()%P;r=0;piv=[]
    for column in range(matrix.shape[1]):
        candidates=np.flatnonzero(matrix[r:,column])
        if not len(candidates):continue
        k=r+int(candidates[0]);matrix[[k,r]]=matrix[[r,k]]
        matrix[r]=matrix[r]*pow(int(matrix[r,column]),-1,P)%P
        for i in range(matrix.shape[0]):
            if i!=r:matrix[i]=(matrix[i]-matrix[i,column]*matrix[r])%P
        piv.append(column);r+=1
        if r==len(matrix):break
    return matrix,piv


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    packet=ROOT/'reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63'
    dependency=json.loads((packet/'report.json').read_text());assert dependency['all_passed']
    for path,value in dependency['input_sha256'].items():assert sha(ROOT/path)==value
    for path,value in dependency['evidence_sha256'].items():assert sha(packet/path)==value
    sep=json.loads((packet/'known-five-separator.json').read_text());rows=sep['selected_original_rows_zero_based']
    F,qids,base,*_=load_fiber()
    fivepath=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(fivepath.read_text())
    vs=sp.symbols('v0:3');slopes=np.array([[residue(sp.sympify(expr).coeff(v))for v in vs]for expr in five['affine_plane_substitution']],dtype=np.int64)
    knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(knownpath.read_text())['known_five_mod_32003']
    ps=[slopes@np.array(point[1:],dtype=np.int64)%P for point in known]
    affine=np.vstack([np.ones(5,dtype=np.int64),np.array(ps,dtype=np.int64).T]);reduced,piv=rref(affine)
    assert len(piv)==4
    free=next(j for j in range(5)if j not in piv);weights=np.zeros(5,dtype=np.int64);weights[free]=1
    for i,j in enumerate(piv):weights[j]=-reduced[i,free]%P
    assert not np.any(affine@weights%P)
    hs=sep['separating_h_values'];V=np.array([[pow(h,j,P)for j in range(5)]for h in hs],dtype=np.int64)
    VI=np.array(sp.Matrix(V.tolist()).inv_mod(P)).astype(np.int64)
    functional=weights@V%P;assert not functional[0]
    normal=next(j for j in range(4,0,-1)if functional[j]);scale=pow(int(functional[normal]),-1,P)
    functional=functional*scale%P;weights=weights*scale%P
    assert normal==4 and int(functional[normal])==1 and np.array_equal(functional@VI%P,weights)
    assert np.all(weights!=0)
    absolute=np.array([(point+np.array(list(map(residue,base)),dtype=np.int64))%P for point in ps],dtype=np.int64)
    qcoords=VI@absolute%P
    constraints=np.vstack([np.array([1,0,0,0,0]),qcoords.T])
    assert not np.any(constraints@functional%P)
    # Four original coordinate functions fix the generic functional uniquely.
    # The universal3plane proves the OTHER42 constraints too, over QQ(beta).
    _,basisrows=rref(constraints[:,:4].T);basisrows=basisrows[:4];assert len(basisrows)==4
    minor=determinant_including_zero(constraints[basisrows,:4]);assert minor
    solved=-np.array(sp.Matrix(constraints[basisrows,:4].tolist()).inv_mod(P)).astype(np.int64)@constraints[basisrows,4]%P
    assert np.array_equal(solved,functional[:4])
    J0=np.zeros((45,45),dtype=np.int64)
    for ii,r in enumerate(rows):
        for mon,c in F[r].items():
            if len(mon)==1:J0[ii,mon[0]]=(J0[ii,mon[0]]+residue(c))%P
    def Hessian_linear(q):
        matrix=np.zeros((45,45),dtype=np.int64)
        for ii,r in enumerate(rows):
            for mon,c in F[r].items():
                if len(mon)==2:
                    a,b=mon;co=residue(c)
                    matrix[ii,a]=(matrix[ii,a]+co*q[b])%P;matrix[ii,b]=(matrix[ii,b]+co*q[a])%P
        return matrix
    rng=random.Random(20261103);direction=[rng.randrange(P)for _ in range(45)]
    leading=Hessian_linear(direction)*pow(2,-1,P)%P;lead_det=determinant_including_zero(leading);assert lead_det
    leading_inverse=np.array(sp.Matrix(leading.tolist()).inv_mod(P)).astype(np.int64);top=[]
    for point in ps:
        constant=(J0+Hessian_linear(point)*pow(2,-1,P))%P;K=leading_inverse@constant%P
        tr=int(np.trace(K))%P;tr2=int(np.trace(K@K%P))%P
        top.append([lead_det,lead_det*tr%P,lead_det*(tr*tr-tr2)*pow(2,-1,P)%P])
    ell_top=[int(weights@np.array([v[j]for v in top],dtype=np.int64))%P for j in range(3)]
    assert ell_top==[0,0,15015]
    diagonal=sep['known_five_jacobian_determinants']
    mu=[sum(pow(h,j,P)*pow(d,-1,P)for h,d in zip(hs,diagonal))%P for j in range(5)]
    Dknown=VI@np.diag(diagonal)%P
    assert np.array_equal(np.array(mu,dtype=np.int64)@Dknown%P,np.ones(5,dtype=np.int64))
    rvalues=[int(c)*int(d)%P for c,d in zip(weights,diagonal)];rnorm=1
    for value in rvalues:assert value;rnorm=rnorm*value%P
    assert rnorm==17591 and lead_det*mu[0]%P==10048
    evidence=dict(prime=P,known_five_affine_coordinate_span_rank=4,
        generic_affine_span_upper_dependency='Universal affine3plane; exact all16beta proof, not this modular rank alone.',
        functional_on_h_coefficient_basis=functional.tolist(),split_functional_weights=weights.tolist(),
        generic_functional_recipe='lambda4=1; solve lambda0..3 using the selected4functions, then all46constraints follow from rank4.',
        selected_function_rows_zero_based=basisrows,selected_constraint4minor_mod_prime=minor,
        constraint_functions='row0=1;row(j+1)=originalq_j',
        frozen_q_coordinates_in_h_basis=qcoords.tolist(),D_coefficient_matrix_at_known_five=Dknown.tolist(),
        trace_d_inverse_h_powers=mu,degree43_separator_values_at_known_five=rvalues,
        degree43_separator_norm_mod_prime=rnorm,q_direction=direction,
        leading_H_over2_determinant=lead_det,
        individual_D_coefficients_at_q_degrees45_44_43=top,
        ell_D_coefficients_at_q_degrees45_44_43=ell_top,
        raw_trace_idempotent_degree45_coefficient=10048,
        raw_trace_averaging_does_not_lower_degree=True,
        determinant_coefficient_system='WriteD_p(q)=D0(q)+hD1(q)+...+h^4D4(q) inE[q].',
        exact_residual_ideal_equalities='I+(e)=I+(D0,D1,D2,D3,D4)=I+(R43)',
        adapted_five_generator_degree_bounds=[45,44,44,44,43],
        R43_definition='R43=ell(D_p)=sum(lambda_j*D_j)',
        R43_generic_degree_exact=43,R43_unit_on_known_five=True,
        residual_emptiness_or_full_length5_claimed=False,generic_polynomials_expanded=False)
    (args.output/'residual43-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),packet/'report.json',packet/'known-five-separator.json',fivepath,knownpath,
        ROOT/'code/s61_affine_polarization_splitter_20261103_v63.py',
        ROOT/'research/notes/TORIC_UNIVERSAL_THREE_PLANE_20261101.md',
        ROOT/'research/notes/S61_FIVE_COVER_PRODUCT_20261101.md',
        ROOT/'research/notes/S61_AFFINE_RESIDUAL43_20261103_v63.md']
    report=dict(id='S61_AFFINE_RESIDUAL43_20261103_v63',all_passed=True,
        exact_global_residual_five_determinant_coefficient_ideal=True,
        exact_global_single_degree43_residual_equation=True,
        unknown_nonreduced_support_retained=True,no_q_or_cyclic_or_H_localization=True,
        beta_only_unit_opens_contain_frozen_center=True,
        finite_recipe_not_expanded_polynomial=True,full_generic_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exactresidualI+(R43)=I+(5Dcoeff)=I+(e);degree43;unknownretained;',report['seconds'])


if __name__=='__main__':main()
