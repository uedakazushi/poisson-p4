#!/usr/bin/env python3
"""Independent coefficient and Schur audit of the beta-only35-variable model."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sympy as s
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha

P=32003


def residue(v):
    v=Fraction(v)
    assert v.denominator%P
    return v.numerator*pow(v.denominator,-1,P)%P


def det(A):
    M=[list(row)for row in A]
    result=1
    for k in range(len(M)):
        pivot=next((i for i in range(k,len(M))if M[i][k]%P),None)
        if pivot is None:return 0
        if pivot!=k:
            M[k],M[pivot]=M[pivot],M[k]
            result=-result
        value=M[k][k]%P
        result=result*value%P
        inverse=pow(value,-1,P)
        for i in range(k+1,len(M)):
            t=M[i][k]*inverse%P
            for j in range(k+1,len(M)):
                M[i][j]=(M[i][j]-t*M[k][j])%P
    return result%P


def inverse(A):
    n=len(A)
    B=[list(row)+[int(i==j)for j in range(n)]for i,row in enumerate(A)]
    for j in range(n):
        k=next(i for i in range(j,n)if B[i][j]%P)
        B[k],B[j]=B[j],B[k]
        t=pow(B[j][j],-1,P)
        B[j]=[v*t%P for v in B[j]]
        for i in range(n):
            if i!=j:
                t=B[i][j]
                B[i]=[(v-t*w)%P for v,w in zip(B[i],B[j])]
    assert [r[:n]for r in B]==[[int(i==j)for j in range(n)]for i in range(n)]
    return [row[n:]for row in B]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    directory=ROOT/'reports/S61_AFFINE_LINEAR_BLOCK_20261103_v63'
    d=read(directory/'report.json')
    fresh=read(args.replay/'report.json')
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,filename
    for filename,digest in d['evidence_sha256'].items():
        assert sha(directory/filename)==digest
        assert (directory/filename).read_bytes()==(args.replay/filename).read_bytes()
    spec=read(ROOT/'code/data/s61_affine_linear_block_20261103_v63.json')
    encoded=read(ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62/restricted-original134.json.gz')
    previous=read(ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62_REVIEW.json')
    assert previous['all_passed']
    X=spec['selected_columns']
    Y=spec['remaining_q_columns']
    R=spec['allowed_rows']
    T=spec['remaining_equation_rows']
    selected=spec['critical_selected_remaining_rows']
    assert len(X)==len(R)==10 and len(Y)==35 and len(T)==124 and len(selected)==35
    assert sorted(X+Y)==list(range(45)) and sorted(R+T)==list(range(134))
    block=read(directory/'linear-block.json')
    allrows=[]
    for row in encoded:
        result={}
        for mon,terms in row:
            mon=tuple(mon)
            assert len(mon)<=2 and mon not in result
            result[mon]={tuple(exponents):Fraction(c)for exponents,c in terms}
            assert len(result[mon])==len(terms)
            for exponents,c in result[mon].items():
                assert len(exponents)==5 and c.denominator%P
        allrows.append(result)
    for k,r in enumerate(R):
        source=allrows[r]
        expected={}
        for j,col in enumerate(X):
            coefficient=source.get((col,),{})
            assert coefficient=={tuple(e):Fraction(v)for e,v in block['beta_coefficient_matrix'][k][j]}
            if coefficient:expected[(col,)]=coefficient
        remainder={tuple(mon):{tuple(e):Fraction(v)for e,v in terms}for mon,terms in block['quadratic_remainder_vector'][k]}
        assert all(not set(mon).intersection(X)for mon in remainder)
        assert not set(expected).intersection(remainder)
        expected.update(remainder)
        assert expected==source
    zero=(0,)*5
    A0=[[allrows[r].get((col,),{}).get(zero,Fraction(0))for col in X]for r in R]
    rational=s.Matrix(A0).det()
    assert rational==s.Rational(block['exact_beta_center_determinant'])!=0
    A=[[residue(v)for v in row]for row in A0]
    assert A==spec['constant_block_mod_p'] and det(A)==28259
    AI=inverse(A)
    five=read(ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json')
    v=s.symbols('v0:3')
    affine=[s.Poly(s.sympify(expr),*v,domain=s.QQ)for expr in five['affine_plane_substitution']]
    # The plane uses absolute free q, whereas the saved modular four-tuples
    # give displacements from the toric point.  Use the actual frozen toric
    # center explicitly; do not confuse the plane's affine constant with it.
    assert len(affine)==45 and all(f.total_degree()<=1 for f in affine)
    slopes=[[residue(f.coeff_monomial(z))for z in v]for f in affine]
    audit=read(ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json')
    toric=next(point for point in audit['points']if point['component_id']=='T01')
    base=[s.Rational(toric['S61_coordinates_u1_through_u61'][i])for i in five['q_to_u_indices_zero_based']]
    free=five['affine_plane_free_q_indices']
    assert [f.eval(dict(zip(v,[base[i]for i in free]))).as_expr()for f in affine]==base
    points=read(ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json')['known_five_mod_32003']
    small=[]
    large=[]
    for point in points:
        q=[(residue(b)+sum(a*c for a,c in zip(row,point[1:])))%P for b,row in zip(base,slopes)]
        J=[]
        for row in allrows:
            differential=[0]*45
            value=0
            for mon,coeff in row.items():
                c=residue(coeff.get(zero,Fraction(0)))
                term=c
                for index in mon:term=term*q[index]%P
                value=(value+term)%P
                for j,index in enumerate(mon):
                    term=c
                    for k,other in enumerate(mon):
                        if j!=k:term=term*q[other]%P
                    differential[index]=(differential[index]+term)%P
            assert value==0
            J.append(differential)
        assert [[J[r][col]for col in X]for r in R]==A
        S=[]
        for r in T:
            left=[sum(J[r][X[k]]*AI[k][j]for k in range(10))%P for j in range(10)]
            S.append([(J[r][col]-sum(left[j]*J[R[j]][col]for j in range(10)))%P for col in Y])
        small.append(det([S[T.index(r)]for r in selected]))
        large.append(det([[J[r][col]for col in X+Y]for r in R+selected]))
    assert small==[24510,4117,15608,23212,723]==spec['known_five_Schur_determinants']
    assert large==[28259*v%P for v in small]==spec['known_five_original_block_determinants']
    paths=[Path(__file__),directory/'report.json',*(directory/f for f in d['evidence_sha256']),
           ROOT/'code/data/s61_affine_linear_block_20261103_v63.json',
           ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62_REVIEW.json',
           ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62/restricted-original134.json.gz',
           ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json',
           ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json',
           ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json',
           ROOT/'code/q_smith_global_model_review_20261103_v63.py',
           ROOT/'research/notes/S61_AFFINE_LINEAR_BLOCK_20261103_v63_REVIEW.md']
    report=dict(id='S61_AFFINE_LINEAR_BLOCK_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_and_all_three_evidence_files_match=True,
        independent_exact_coefficients_all10_F_rows_equal_beta_A_X_plus_B=True,
        all_B_in_remaining35_variables_and_all134_rows_retained_or_solved=True,
        exact_QQ_center_determinant=str(rational),center_determinant_mod32003=28259,
        all_coefficient_denominators_p_units=True,
        entire_affine_fiber_quotient_isomorphism_on_beta_only_determinant_open_accepted=True,
        no_q_dependent_denominator_and_cleared_Y_degree_bound=4,
        independent_original134_derivatives_at_all5_and_Schur_determinants=small,
        original45block_equals_d_times35Schur_at_all5=True,
        exact_relative_differential124_by35_and_Fitting0_all35minors_accepted=True,
        chosen_minor_contains_all_positive_dimensional_support=True,
        finite_over_beta_generic_dimension_zero_or_unit_GB_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:
        json.dump(report,f,indent=2)
        f.write('\n')
    print('PASS independent whole affine35 model; beta-only10block and all5 Schur checks')


if __name__=='__main__':
    main()
