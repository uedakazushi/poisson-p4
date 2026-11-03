#!/usr/bin/env python3
"""A linear42-square determinant defines the entire unknown affine residual.

This retains the unresolved residual.  It does not infer emptiness from the
finite-field calculations.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from s61_affine_fiber_macaulay import load_fiber
from s61_affine_normal_projection_20261102 import matrices
from s61_affine_polarization_splitter_20261103_v63 import determinant_including_zero,residue

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def add(poly,mon,coefficient):
    value=poly.get(mon,Fraction(0))+coefficient
    if value:poly[mon]=value
    else:poly.pop(mon,None)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    specpath=ROOT/'code/data/s61_affine_residual42_20261103_v63.json';spec=json.loads(specpath.read_text())
    selected=spec['selected_residual_rows_zero_based'];assert len(selected)==len(set(selected))==42
    F,qids,center,original,pointpath=load_fiber()
    fivepath=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(fivepath.read_text())
    zs=sp.symbols('v0:3');expressions=list(map(sp.sympify,five['affine_plane_substitution']))
    embedding=[[Fraction(expression.coeff(z))for z in zs]for expression in expressions]
    assert embedding[:3]==[[Fraction(int(i==j))for j in range(3)]for i in range(3)]
    linearimages=[]
    for i,row in enumerate(embedding):
        image={(j,):coefficient for j,coefficient in enumerate(row)if coefficient}
        if i>=3:image[i,]=Fraction(1)
        linearimages.append(image)
    transformed=[]
    for polynomial in F:
        out={}
        for mon,c in polynomial.items():
            if not mon:add(out,(),c)
            elif len(mon)==1:
                for term,value in linearimages[mon[0]].items():add(out,term,c*value)
            else:
                assert len(mon)==2
                for a,ca in linearimages[mon[0]].items():
                    for b,cb in linearimages[mon[1]].items():add(out,tuple(sorted(a+b)),c*ca*cb)
        transformed.append(out)
    change=[[Fraction(c)for c in row]for row in five['all134_restricted_span_coefficients']]
    factors=[];normal_residuals=[]
    for i in range(5,134):
        residual=dict(transformed[i])
        for j,c in enumerate(change[i]):
            if c:
                for mon,value in transformed[j].items():add(residual,mon,-c*value)
        A=[{}for _ in range(42)]
        for mon,c in residual.items():
            normal=next((v for v in mon if v>=3),None);assert normal is not None
            remaining=list(mon);remaining.remove(normal);assert len(remaining)<=1
            add(A[normal-3],tuple(remaining),c)
        replay={}
        for j,column in enumerate(A):
            for mon,c in column.items():add(replay,tuple(sorted((j+3,)+mon)),c)
        assert replay==residual
        assert all(c.denominator%P for column in A for c in column.values())
        factors.append(A);normal_residuals.append(residual)
    # Separate existing finite-field engine checks every factor coefficient.
    cross,quad,plane_minor,_=matrices();assert plane_minor==spec['plane_restriction_minor']
    for r,A in enumerate(factors):
        for j,column in enumerate(A):
            assert residue(column.get((),0))==int(cross[r,j,0])
            for z in range(3):assert residue(column.get((z,),0))==int(cross[r,j,z+1])
            for k in range(42):assert residue(column.get((k+3,),0))==int(quad[r,j,k])
    knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(knownpath.read_text())['known_five_mod_32003'];values=[]
    for p in known:
        zz=np.array([1,*p[1:]],dtype=np.int64);A=np.einsum('rnk,k->rn',cross,zz)%P
        values.append(determinant_including_zero(A[selected,:]))
    assert values==spec['known_five_determinants']==[1680,23338,12913,13967,26226]
    norm=1
    for value in values:assert value;norm=norm*value%P
    rng=np.random.default_rng(20261103);zz=rng.integers(0,P,size=4);zz[0]=0;nn=rng.integers(0,P,size=42)
    leading=(np.einsum('rnk,k->rn',cross,zz)+np.einsum('rjk,k->rj',quad,nn))%P
    leading_det=determinant_including_zero(leading[selected,:]);assert leading_det==spec['top_degree42_direction_determinant']==31743
    encode=lambda poly:[[list(mon),str(c)]for mon,c in sorted(poly.items())]
    sparse=dict(coordinate_rule='deltaq_i=sum_j embedding_ij*z_j + n_(i-3) fori>=3; first3deltaq=z',
        embedding=[[str(c)for c in row]for row in embedding],
        selected_residual_rows_zero_based=selected,
        pure_plane_subtraction_coefficients=[[str(c)for c in row]for row in change[5:]],
        all129_linear_factor_rows=[[[[list(mon),str(c)]for mon,c in sorted(column.items())]for column in row]for row in factors])
    (args.output/'linear129x42-QQ.json.gz').write_bytes(gzip.compress(json.dumps(sparse,separators=(',',':')).encode(),mtime=0))
    (args.output/'residual42-certificate.json').write_text(json.dumps(dict(
        exact_QQ_factorization_rows=129,exact_QQ_all134_inputs_accounted_for=True,
        identity='F_i-sum_j change_ij F_j = sum_k A_ik(q)*n_k, i5..133,j0..4',
        selected_linear_matrix_shape=[42,42],selected_residual_rows_zero_based=selected,
        known_five_determinants_mod_prime=values,normal_determinant_norm_mod_prime=norm,
        leading_degree42_direction_z=zz.tolist(),leading_degree42_direction_n=nn.tolist(),
        leading_degree42_determinant_mod_prime=leading_det,
        generic_recipe='Use universal3plane frame and5independent restrictionquadrics; subtracttheirspan; factor byfirstnormal variable; selectsaved42rows.',
        generic_coefficients_finite_inverse_recipe=True,generic_expanded_polynomial_saved=False,
        exact_ideal_equality='I_res=I+(R42)=I+(e), R42=det A_selected',
        adjugate_identity='R42*n = adj(A_selected)*F_selected_residual inI^42',
        generic_q_degree_exact=42,no_q_or_cyclic_or_H_localization=True,
        unknown_support_and_nonreduced_structure_retained=True,
        full_generic_exhaustion_claimed=False),indent=2)+'\n')
    dependency=ROOT/'reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63/report.json'
    old=json.loads(dependency.read_text());assert old['all_passed']
    for path,value in old['input_sha256'].items():assert sha(ROOT/path)==value
    for path,value in old['evidence_sha256'].items():assert sha(dependency.parent/path)==value
    inputs=[Path(__file__),specpath,original,pointpath,fivepath,knownpath,dependency,
        ROOT/'code/s61_affine_normal_projection_20261102.py',
        ROOT/'research/notes/TORIC_UNIVERSAL_THREE_PLANE_20261101.md',
        ROOT/'research/notes/S61_FIVE_COVER_PRODUCT_20261101.md',
        ROOT/'research/notes/S61_AFFINE_RESIDUAL42_20261103_v63.md']
    report=dict(id='S61_AFFINE_RESIDUAL42_20261103_v63',all_passed=True,
        all129_exact_QQ_linear_factorizations_replayed=True,
        full_generic_linear_matrix_from_universal_plane_recipe=True,
        one_degree42_equation_defines_full_unknown_residual_exactly=True,
        coefficient_degree_one_matrix_size42=True,beta_only_unit_open_contains_frozen_center=True,
        no_unknown_q_support_removed=True,unit_or_exhaustion_computed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS129QQlinear-factor;exactresidualI+(R42);degree42;unitOPEN;',report['seconds'])


if __name__=='__main__':main()
