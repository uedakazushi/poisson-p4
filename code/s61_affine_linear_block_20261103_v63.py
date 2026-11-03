#!/usr/bin/env python3
"""Whole affine35-variable model and exact Schur differential on the admissible slice.

Only a beta-dependent determinant is inverted.  No finite-over-beta conclusion
or Gröbner calculation is made.
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
from s61_nonsquare_cubic_row_basis_review_20261103_v62 import determinant_mod

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def residue(value):
    value=Fraction(value);assert value.denominator%P
    return value.numerator*pow(value.denominator,-1,P)%P


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    specpath=ROOT/'code/data/s61_affine_linear_block_20261103_v63.json'
    spec=json.loads(specpath.read_text());X=spec['selected_columns'];RR=spec['allowed_rows']
    Y=spec['remaining_q_columns'];TT=spec['remaining_equation_rows'];TR=spec['critical_selected_remaining_rows']
    assert len(X)==len(RR)==10 and len(Y)==35 and len(TT)==124 and len(TR)==35
    assert sorted(X+Y)==list(range(45))and sorted(RR+TT)==list(range(134))
    packet=ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62'
    old=json.loads((packet/'report.json').read_text());assert old['all_passed']
    for path,value in old['input_sha256'].items():assert sha(ROOT/path)==value
    for path,value in old['evidence_sha256'].items():assert sha(packet/path)==value
    encoded=json.loads(gzip.decompress((packet/'restricted-original134.json.gz').read_bytes()))
    assert len(encoded)==134
    coefficient_zero=[]
    for row in encoded:
        out={}
        for mon,coef in row:
            assert len(mon)<=2
            for exponent,c in coef:
                cc=Fraction(c);assert cc.denominator%P
            out[tuple(mon)]=next((Fraction(c)for e,c in coef if not any(e)),Fraction(0))
        coefficient_zero.append(out)
    block=[];remainders=[]
    for r in RR:
        row=dict((tuple(mon),coef)for mon,coef in encoded[r]);arow=[]
        for mon in row:
            assert not(len(mon)==2 and any(c in X for c in mon))
        for c in X:arow.append(row.pop((c,),[]))
        assert all(not any(c in X for c in mon)for mon in row)
        block.append(arow);remainders.append([[list(mon),coef]for mon,coef in row.items()])
    A0=sp.Matrix([[sp.Rational(coefficient_zero[r].get((c,),0))for c in X]for r in RR])
    exact_det=sp.factor(A0.det());assert exact_det!=0
    AA=np.array([[residue(v)for v in row]for row in A0.tolist()],dtype=np.int64)
    assert AA.tolist()==spec['constant_block_mod_p']
    detA=determinant_mod(AA.tolist(),P);assert detA==spec['detA_mod_p']==28259
    assert residue(exact_det)==detA
    AI=np.array(sp.Matrix(AA.tolist()).inv_mod(P)).astype(np.int64)
    # Independent original134 derivatives at each saved known point.  The
    # top block remains AA because its ten columns contain no quadratic term.
    F,qids,base,*_=load_fiber()
    fivepath=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    five=json.loads(fivepath.read_text());vs=sp.symbols('v0:3')
    slopes=np.array([[residue(sp.sympify(expr).coeff(v))for v in vs]for expr in five['affine_plane_substitution']],dtype=np.int64)
    knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    points=json.loads(knownpath.read_text())['known_five_mod_32003'];small=[];large=[]
    sel=[TT.index(r)for r in TR]
    for point in points:
        delta=slopes@np.array(point[1:],dtype=np.int64)%P;J=np.zeros((134,45),dtype=np.int64)
        for r,polynomial in enumerate(F):
            ev=0
            for mon,c in polynomial.items():
                co=residue(c);value=co
                for i in mon:value=value*int(delta[i])%P
                ev=(ev+value)%P
                if len(mon)==1:J[r,mon[0]]=(J[r,mon[0]]+co)%P
                elif len(mon)==2:
                    a,b=mon;J[r,a]=(J[r,a]+co*delta[b])%P;J[r,b]=(J[r,b]+co*delta[a])%P
            assert ev==0
        assert np.array_equal(J[np.ix_(RR,X)],AA)
        Schur=(J[np.ix_(TT,Y)]-J[np.ix_(TT,X)]@AI%P@J[np.ix_(RR,Y)])%P
        small.append(determinant_mod(Schur[sel,:].tolist(),P))
        large.append(determinant_mod(J[np.ix_(RR+TR,X+Y)].tolist(),P))
    assert small==spec['known_five_Schur_determinants']==[24510,4117,15608,23212,723]
    assert large==spec['known_five_original_block_determinants']==[detA*d%P for d in small]
    assert all(small)
    (args.output/'linear-block.json').write_text(json.dumps(dict(
        original_q_columns_zero_based=X,original_equation_rows_zero_based=RR,
        remaining_q_columns_zero_based=Y,remaining_equation_rows_zero_based=TT,
        beta_coefficient_matrix=block,quadratic_remainder_vector=remainders,
        exact_beta_center_determinant=str(exact_det),determinant_mod32003=detA,
        all_selected_equations_exactly_A_times_X_plus_B=True,
        all_block_entries_independent_of_all_q=True,
        all_B_entries_in_remaining35q_only=True,
        model_recipe='d=det(A); X=-adj(A)*B/d; keep every other original equation after this substitution.',
        base_ring='Z_(32003)[ss1,...,ss5,1/d], using the existing p-integral rational center',
        finite_over_beta_asserted=False),indent=2)+'\n')
    (args.output/'Schur-critical.json').write_text(json.dumps(dict(
        selected_remaining_rows_zero_based=TR,
        definition='S=J_TY-J_TX*A^(-1)*J_RY, followed by X=-A^(-1)*B',
        differential_presentation_shape=[124,35],
        chosen_critical_minor_shape=[35,35],known_five_determinants_mod32003=small,
        known_five_original45block_determinants_mod32003=large,
        all_positive_dimensional_support_contained_in_chosen_minor_zero=True,
        original_global_E_replaced_by_this_new_minor=False,
        exact_full_critical_ideal='all35x35minors of the124x35Schur matrix; Fitting_0(Omega)',
        original45det_equals_d_times_Schur35det=True,
        no_improved_degree_bound_or_unit_result_claimed=True),indent=2)+'\n')
    # A fully specified finite recipe, retaining all134 equations until the
    # ten solved equations are removed.  We do not execute a field GB.
    def coefficient_text(coeff):
        return '+'.join('('+str(Fraction(c))+')'+''.join('*ss'+str(i+1)+(('^'+str(e))if e!=1 else'')for i,e in enumerate(exponents)if e)for exponents,c in coeff)or'0'
    equations=[]
    for row in encoded:
        equations.append('+'.join('('+coefficient_text(c)+')'+''.join('*q'+str(i)for i in mon)for mon,c in row)or'0')
    text='A=QQ[ss1,ss2,ss3,ss4,ss5];K=frac A;R=K['+','.join(f'q{i}'for i in range(45))+'];\n'
    text+='FF=matrix{{'+','.join(equations)+'}};\n'
    text+='xx={'+','.join('q'+str(i)for i in X)+'};yy={'+','.join('q'+str(i)for i in Y)+'};\n'
    text+='rowsX={'+','.join(map(str,RR))+'};rowsT={'+','.join(map(str,TT))+'};\n'
    text+='JJ=transpose(jacobian FF);\nAA=submatrix(JJ,rowsX,{'+','.join(map(str,X))+'});\n'
    text+='BB=transpose(submatrix(FF,{0},rowsX))-AA*transpose(matrix{xx});\n'
    text+='xxsol=-inverse AA*BB;\nS=K['+','.join('y'+str(i)for i in range(35))+'];\n'
    images=[]
    for i in range(45):images.append('xxsol_('+str(X.index(i))+',0)'if i in X else'q'+str(i))
    text+='-- Substitute xxsol for xx first; all other original124 equations are retained.\n'
    text+='elim=map(R,R,{'+','.join(images)+'});\n'
    text+='toS=map(S,R,{'+','.join('0'if i in X else'y'+str(Y.index(i))for i in range(45))+'});\n'
    text+='G=toS elim submatrix(FF,{0},rowsT);\nI35=ideal G;\n'
    text+='-- Exact model definition only; no gb or unit result is asserted.\n'
    (args.output/'exact35-model.m2').write_text(text)
    inputs=[Path(__file__),specpath,packet/'report.json',packet/'restricted-original134.json.gz',
        ROOT/'code/s61_affine_fiber_macaulay.py',fivepath,knownpath,
        ROOT/'research/notes/S61_AFFINE_LINEAR_BLOCK_20261103_v63.md']
    report=dict(id='S61_AFFINE_LINEAR_BLOCK_20261103_v63',all_passed=True,
        whole_affine_fiber_on_admissible_slice_isomorphic_to35variable_model=True,
        beta_only_localization_contains_frozen_center=True,all134_original_equations_accounted_for=True,
        solved_original_equations=10,retained_original_equations=124,remaining_q_variables=35,
        no_q_dependent_localization=True,exact_p_integral_base_localization=True,
        finite_over_beta_proved=False,full_generic_dimension_zero_proved=False,
        critical_unit_or_GB_computed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS whole35variable model;betaonly10block;124x35Schur;notfiniteoverbeta;',report['seconds'])


if __name__=='__main__':main()
